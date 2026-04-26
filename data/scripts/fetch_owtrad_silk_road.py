#!/usr/bin/env python3
"""
fetch_owtrad_silk_road.py
=========================
Download and parse OWTRAD GIS datasets to build a comprehensive Silk Road
route network for the Gerizim beru-grid study.

SOURCE
------
Old World Trade Routes (OWTRAD) GIS project by T. Matthew Ciolek.
Research School of Pacific and Asian Studies, Australian National University.

  URL:     http://www.ciolek.com/owtrad.html
  Licence: Open Publication License v1.0 (http://www.opencontent.org/opl.shtml)
  Cite:    Ciolek, T.M. (2004+). Old World Trade Routes (OWTRAD) Project.
           www.ciolek.com/owtrad.html  Canberra: Asia Pacific Research Online.

DATASETS DOWNLOADED
-------------------
  tmcZCAm1000  Central Asia 1–1400 CE          lon 25–135°E, lat 25–50°N
  tmcIRa0100   Iran / Persia 50 BCE–300 CE      lon 30–80°E,  lat 10–45°N
  tmcCNm1000   NW China 100–1400 CE             lon 74–94°E,  lat 36–44°N
  tmcKGa0100d  Kyrgyzstan 100 BCE–1400 CE       lon 65–82°E,  lat 36–45°N  (local)

OUTPUT
------
  data/store/silk_road/owtrad_routes.csv
    Columns: dataset, lon1, lat1, lon2, lat2, node1, node2, country1, country2

  data/store/silk_road/owtrad_nodes.csv
    Columns: dataset, name, lon, lat, country, node_id

USAGE
-----
    python3 data/scripts/fetch_owtrad_silk_road.py          # download & parse
    python3 data/scripts/fetch_owtrad_silk_road.py --local  # parse local files only
    python3 data/scripts/fetch_owtrad_silk_road.py --validate
"""

import argparse
import csv
import hashlib
import io
import re
import sys
import time
import zipfile
from pathlib import Path

try:
    import requests
    requests.packages.urllib3.disable_warnings()
except ImportError:
    sys.exit("pip install requests  (required)")

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT        = Path(__file__).resolve().parents[2]
CACHE_DIR   = ROOT / "data" / "store" / "silk_road" / "owtrad_cache"
OUTPUT_DIR  = ROOT / "data" / "store" / "silk_road"
ROUTES_CSV  = OUTPUT_DIR / "owtrad_routes.csv"
NODES_CSV   = OUTPUT_DIR / "owtrad_nodes.csv"
LOCAL_DIR   = ROOT / "data" / "owtrad-gis-tmcKGa0100d"   # existing local dataset

USER_AGENT  = "gerizim-beru-grid/1.0 (academic research; owtrad-silk-road)"
BASE_URL    = "http://www.ciolek.com/OWTRAD/DATA/interchange"

# ── Datasets to download ──────────────────────────────────────────────────────
# Ordered west-to-east by coverage; deduplication happens at parse time.
DATASETS = [
    {
        "id":      "tmcZCAm1000",
        "desc":    "Central Asia 1–1400 CE",
        "bbox":    (25, 135, 25, 50),  # lon_min, lon_max, lat_min, lat_max
        "url":     f"{BASE_URL}/owtrad-gis-tmcZCAm1000.zip",
        "cite":    "Ciolek 2004+, tmcZCAm1000",
    },
    {
        "id":      "tmcIRa0100",
        "desc":    "Iran/Persia 50 BCE–300 CE",
        "bbox":    (30, 80, 10, 45),
        "url":     f"{BASE_URL}/owtrad-gis-tmcIRa0100.zip",
        "cite":    "Ciolek 2004+, tmcIRa0100",
    },
    {
        "id":      "tmcCNm1000",
        "desc":    "NW China 100–1400 CE",
        "bbox":    (74, 94, 36, 44),
        "url":     f"{BASE_URL}/owtrad-gis-tmcCNm1000.zip",
        "cite":    "Ciolek 2004+, tmcCNm1000",
    },
]

ROUTE_FIELDNAMES = ["dataset", "lon1", "lat1", "lon2", "lat2",
                    "node1", "node2", "country1", "country2"]
NODE_FIELDNAMES  = ["dataset", "name", "lon", "lat", "country", "node_id"]


# ── Download ──────────────────────────────────────────────────────────────────

def fetch_zip(url: str, cache_path: Path, retries: int = 3) -> bytes:
    if cache_path.exists():
        print(f"    [cached] {cache_path.name}")
        return cache_path.read_bytes()
    for attempt in range(retries):
        try:
            print(f"    Downloading {url} ...", flush=True)
            r = requests.get(url, headers={"User-Agent": USER_AGENT},
                             timeout=120, verify=False, stream=True)
            r.raise_for_status()
            chunks = []
            total = 0
            for chunk in r.iter_content(65536):
                chunks.append(chunk)
                total += len(chunk)
                print(f"\r    {total/1024:.0f} KB", end="", flush=True)
            print()
            data = b"".join(chunks)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_bytes(data)
            return data
        except Exception as e:
            if attempt < retries - 1:
                print(f"    [retry {attempt+1}/{retries}: {e}]")
                time.sleep(10)
            else:
                raise
    return b""


# ── MIF/MID parsers ───────────────────────────────────────────────────────────

def _parse_mif_coords(mif_text: str) -> list:
    """Extract list of (lon1, lat1, lon2, lat2) from a routes .mif file."""
    segs = []
    # Line features: "Line lon1 lat1 lon2 lat2" (possibly with leading whitespace)
    pattern = re.compile(
        r'\bLine\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', re.I
    )
    for m in pattern.finditer(mif_text):
        try:
            lon1, lat1, lon2, lat2 = (float(m.group(i)) for i in range(1, 5))
            segs.append((lon1, lat1, lon2, lat2))
        except ValueError:
            pass
    return segs


def _parse_mid_routes(mid_text: str, dataset_id: str) -> list:
    """
    Parse a routes .mid file.
    Column order (from MIF header): NODE1, COUNTRY1, NODE2, COUNTRY2, DETAIL,
    USES, TYPE, ROLE, GOODS1-3, DIR, DIST, TRAVMODE, TRAVTIME, EARLYDATE,
    LATEDATE, DATAQLTY, SRC, LONG1, LAT1, COORDSRC1, LONG2, LAT2, COORDSRC2,
    NODEID1, NODEID2, PROBL, DATAID
    """
    rows = []
    for line in mid_text.splitlines():
        line = line.strip()
        if not line:
            continue
        # CSV-quoted fields
        parts = next(csv.reader([line]))
        if len(parts) < 25:
            continue
        try:
            node1    = parts[0].strip()
            country1 = parts[1].strip()
            node2    = parts[2].strip()
            country2 = parts[3].strip()
            lon1     = float(parts[19])
            lat1     = float(parts[20])
            lon2     = float(parts[22])
            lat2     = float(parts[23])
            rows.append({
                "dataset":  dataset_id,
                "lon1":     f"{lon1:.6f}",
                "lat1":     f"{lat1:.6f}",
                "lon2":     f"{lon2:.6f}",
                "lat2":     f"{lat2:.6f}",
                "node1":    node1,
                "node2":    node2,
                "country1": country1,
                "country2": country2,
            })
        except (ValueError, IndexError):
            continue
    return rows


def _parse_mid_nodes(mid_text: str, dataset_id: str) -> list:
    """
    Parse a nodes .mid file.
    Column order: NODE1, LONG1, LAT1, COUNTRY1, DETAIL, DATASTATUS, NODEID1, DATAID
    """
    rows = []
    for line in mid_text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = next(csv.reader([line]))
        if len(parts) < 7:
            continue
        try:
            name    = parts[0].strip()
            lon     = float(parts[1])
            lat     = float(parts[2])
            country = parts[3].strip()
            node_id = parts[6].strip()
            rows.append({
                "dataset": dataset_id,
                "name":    name,
                "lon":     f"{lon:.6f}",
                "lat":     f"{lat:.6f}",
                "country": country,
                "node_id": node_id,
            })
        except (ValueError, IndexError):
            continue
    return rows


# ── Read from zip ─────────────────────────────────────────────────────────────

def parse_zip(data: bytes, dataset_id: str) -> tuple:
    """Returns (route_rows, node_rows)."""
    zf  = zipfile.ZipFile(io.BytesIO(data))
    names = zf.namelist()

    # Find routes and nodes files (case-insensitive)
    routes_mid = next((n for n in names if "routes" in n.lower() and n.endswith(".mid")), None)
    nodes_mid  = next((n for n in names if "nodes"  in n.lower() and n.endswith(".mid")), None)

    route_rows = []
    node_rows  = []

    if routes_mid:
        mid_text = zf.read(routes_mid).decode("latin-1", errors="replace")
        route_rows = _parse_mid_routes(mid_text, dataset_id)

    if nodes_mid:
        mid_text = zf.read(nodes_mid).decode("latin-1", errors="replace")
        node_rows = _parse_mid_nodes(mid_text, dataset_id)

    return route_rows, node_rows


# ── Parse local tmcKGa0100d ───────────────────────────────────────────────────

def parse_local_kga() -> tuple:
    """Parse the pre-existing local Kyrgyzstan dataset."""
    ds = "tmcKGa0100d"
    route_rows = []
    node_rows  = []

    routes_mid = LOCAL_DIR / "tmcKGa0100d_routes.mid"
    nodes_mid  = LOCAL_DIR / "tmcKGa0100d_nodes.mid"

    if routes_mid.exists():
        route_rows = _parse_mid_routes(routes_mid.read_text("latin-1"), ds)
    if nodes_mid.exists():
        node_rows  = _parse_mid_nodes(nodes_mid.read_text("latin-1"),  ds)

    return route_rows, node_rows


# ── Deduplication ─────────────────────────────────────────────────────────────

def dedup_routes(rows: list) -> list:
    """Remove duplicate edges (same endpoint pair regardless of dataset)."""
    seen = set()
    out  = []
    for r in rows:
        key = (round(float(r["lon1"]), 4), round(float(r["lat1"]), 4),
               round(float(r["lon2"]), 4), round(float(r["lat2"]), 4))
        rev = (key[2], key[3], key[0], key[1])
        if key not in seen and rev not in seen:
            seen.add(key)
            out.append(r)
    return out


def dedup_nodes(rows: list) -> list:
    seen = set()
    out  = []
    for r in rows:
        key = (round(float(r["lon"]), 3), round(float(r["lat"]), 3))
        if key not in seen:
            seen.add(key)
            out.append(r)
    return out


# ── Write CSV ─────────────────────────────────────────────────────────────────

def write_routes(rows: list, fetch_date: str) -> None:
    with open(ROUTES_CSV, "w", newline="", encoding="utf-8") as f:
        f.write("# owtrad_routes.csv — OWTRAD Silk Road route network (edge list)\n")
        f.write("# DO NOT HAND-EDIT. Regenerate: python3 data/scripts/fetch_owtrad_silk_road.py\n")
        f.write(f"# Generated: {fetch_date}\n")
        f.write("# Source: Ciolek, T.M. (2004+). OWTRAD Project. www.ciolek.com/owtrad.html\n")
        f.write("# Datasets: tmcZCAm1000, tmcIRa0100, tmcCNm1000, tmcKGa0100d\n")
        f.write("# Licence: Open Publication License v1.0\n")
        f.write(f"# N = {len(rows)} edges\n")
        w = csv.DictWriter(f, fieldnames=ROUTE_FIELDNAMES)
        w.writeheader()
        w.writerows(rows)


def write_nodes(rows: list, fetch_date: str) -> None:
    with open(NODES_CSV, "w", newline="", encoding="utf-8") as f:
        f.write("# owtrad_nodes.csv — OWTRAD Silk Road network nodes\n")
        f.write("# DO NOT HAND-EDIT. Regenerate: python3 data/scripts/fetch_owtrad_silk_road.py\n")
        f.write(f"# Generated: {fetch_date}\n")
        f.write("# Source: Ciolek, T.M. (2004+). OWTRAD Project. www.ciolek.com/owtrad.html\n")
        f.write(f"# N = {len(rows)} nodes\n")
        w = csv.DictWriter(f, fieldnames=NODE_FIELDNAMES)
        w.writeheader()
        w.writerows(rows)


# ── Validate ──────────────────────────────────────────────────────────────────

def run_validate() -> int:
    for path in [ROUTES_CSV, NODES_CSV]:
        if not path.exists():
            print(f"  MISSING: {path.relative_to(ROOT)}")
            return 1
        rows = []
        with open(path) as f:
            non_comment = (l for l in f if not l.startswith("#"))
            rows = list(csv.DictReader(non_comment))
        sha = hashlib.sha256(path.read_bytes()).hexdigest()[:16]
        print(f"  {path.name}: {len(rows)} rows  SHA-256={sha}...")
    return 0


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    from datetime import datetime, timezone
    parser = argparse.ArgumentParser()
    parser.add_argument("--local",    action="store_true", help="Parse local files only (no download)")
    parser.add_argument("--validate", action="store_true")
    args = parser.parse_args()

    if args.validate:
        sys.exit(run_validate())

    fetch_date = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    print("=" * 72)
    print("  OWTRAD SILK ROAD ROUTE DOWNLOADER")
    print(f"  Run at: {fetch_date}")
    print("=" * 72)

    all_routes = []
    all_nodes  = []

    # 1. Local Kyrgyzstan dataset
    print("\n  Parsing local tmcKGa0100d (Kyrgyzstan) ...")
    r, n = parse_local_kga()
    print(f"    Routes: {len(r)}  Nodes: {len(n)}")
    all_routes += r
    all_nodes  += n

    # 2. Remote datasets
    if not args.local:
        for ds in DATASETS:
            print(f"\n  Dataset: {ds['id']} — {ds['desc']}")
            cache = CACHE_DIR / f"owtrad-gis-{ds['id']}.zip"
            data  = fetch_zip(ds["url"], cache)
            r, n  = parse_zip(data, ds["id"])
            print(f"    Routes: {len(r)}  Nodes: {len(n)}")
            all_routes += r
            all_nodes  += n

    # 3. Deduplicate
    all_routes = dedup_routes(all_routes)
    all_nodes  = dedup_nodes(all_nodes)
    print(f"\n  After deduplication: {len(all_routes)} route edges, {len(all_nodes)} nodes")

    # 4. Write
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_routes(all_routes, fetch_date)
    write_nodes(all_nodes,  fetch_date)

    rsha = hashlib.sha256(ROUTES_CSV.read_bytes()).hexdigest()
    nsha = hashlib.sha256(NODES_CSV.read_bytes()).hexdigest()
    print(f"\n  Written: {ROUTES_CSV.relative_to(ROOT)}  ({len(all_routes)} edges)  SHA-256={rsha[:16]}...")
    print(f"  Written: {NODES_CSV.relative_to(ROOT)}   ({len(all_nodes)} nodes)   SHA-256={nsha[:16]}...")

    # 5. Summary stats
    lons1 = [float(r["lon1"]) for r in all_routes]
    lons2 = [float(r["lon2"]) for r in all_routes]
    lats1 = [float(r["lat1"]) for r in all_routes]
    lats2 = [float(r["lat2"]) for r in all_routes]
    all_lons = lons1 + lons2
    all_lats = lats1 + lats2
    print(f"\n  Network extent: lon {min(all_lons):.1f}–{max(all_lons):.1f}°E, "
          f"lat {min(all_lats):.1f}–{max(all_lats):.1f}°N")

    datasets_used = sorted({r["dataset"] for r in all_routes})
    print(f"  Datasets: {', '.join(datasets_used)}")
    print("\n  Done.")


if __name__ == "__main__":
    main()