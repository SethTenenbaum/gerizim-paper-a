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
  Remote (downloaded from interchange/):
  Mediterranean & Turkey:
    tmcXMEm1400  Mediterranean ~1400 CE          lon  -10– 40°E, lat 25–55°N
    tmcTRm1200a  Turkey 1200 CE (part a)         lon   25– 45°E, lat 35–45°N
    tmcTRm1200b  Turkey 1200 CE (part b)         lon   25–120°E, lat 35–45°N
    tmcTRm1300   Turkey / Anatolia 1300 CE       lon   12– 48°E, lat 35–45°N
  Persia / Iran:
    tmcIRa0100   Iran/Persia 50 BCE–300 CE       lon   30– 80°E, lat 10–45°N
    tmcIRa0500   Iran/Persia 300–700 CE          lon   30– 80°E, lat 10–45°N
  Central Asia:
    tmcZCAm0200  Iran & China 200 BCE–500 CE     lon   25–120°E, lat 15–55°N
    tmcZCAm0600  Mediterranean–China 200 BCE–    lon    0–120°E, lat 15–55°N
                 1400 CE
    tmcZCAm0800  Central Asia 700–1000 CE        lon   25–120°E, lat 25–55°N
    tmcZCAm1000  Central Asia 1–1400 CE          lon   25–135°E, lat 25–50°N
    tmcZCAm1000a Central Asia 1–1400 CE (var a)  lon   25–135°E, lat 25–50°N
    tmcTMm1100   Turkmenistan 1100 CE            lon   50– 64°E, lat 35–45°N
    tmcTJm0400a  Tajikistan 400 CE (part a)      lon   65– 82°E, lat 36–45°N
    tmcTJm0400b  Tajikistan 400 CE (part b)      lon   65– 82°E, lat 36–45°N
    tmcKGa0100a  Kyrgyzstan (part a)             lon   65– 82°E, lat 36–45°N
    tmcKGa0100b  Kyrgyzstan (part b)             lon   65– 82°E, lat 36–45°N
    tmcKGa0100c  Kyrgyzstan (part c)             lon   65– 82°E, lat 36–45°N
    tmcKGa0100e  Kyrgyzstan (part e)             lon   65– 82°E, lat 36–45°N
  Middle East & India:
    tmcZMEm1300  Middle East & India 1300–       lon   25– 90°E, lat  5–45°N
                 1600 CE
    tmcINm1550   India 1550–1650 CE              lon   65– 90°E, lat  8–35°N
  China (all periods):
    tmcCNm1000   NW China 100–1400 CE            lon   74– 94°E, lat 36–44°N
    tmcCNm0620   China 620 CE (Tang routes)      lon   87–104°E, lat 30–45°N
    tmcCNm0680a  China 680 CE (part a)           lon   88–104°E, lat 30–45°N
    tmcCNm0680b  China 680 CE (part b)           lon   88–104°E, lat 30–45°N
    tmcCNm0680c  China 680 CE (part c)           lon   88–104°E, lat 30–45°N
    tmcCNm1500   China Ming courier 1368–1644 CE lon  102–121°E, lat 20–50°N
    tmcCNm1700   China Qing routes 1644–1800 CE  lon  115–121°E, lat 25–45°N
    tmcCNm1850   China routes 1800–1900 CE       lon   98–117°E, lat 25–45°N
    tmcCNm1920   China routes 1900–1949 CE       lon   80–117°E, lat 20–50°N
  Southeast Asia:
    tmcKHm1200   Khmer / Cambodia 1200–1300 CE   lon  100–110°E, lat 10–20°N

  Local (pre-existing):
  tmcKGa0100d  Kyrgyzstan 100 BCE–1400 CE        lon   65– 82°E, lat 36–45°N

NOTE: Not all datasets have a zip in the interchange directory. fetch_zip()
will skip with a warning on HTTP 404 rather than crashing.

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
from typing import Optional

try:
    import requests
    # Suppress InsecureRequestWarning from verify=False — OWTRAD server has
    # a self-signed / expired certificate; verify=False is intentional.
    try:
        import urllib3
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    except Exception:
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
# Not all datasets have a zip in the interchange directory — fetch_zip()
# returns None on 404 and those are silently skipped.
DATASETS = [
    # ── Mediterranean & Turkey ─────────────────────────────────────────────────
    {
        "id":   "tmcXMEm1400",
        "desc": "Mediterranean ~1400 CE",
        "bbox": (-10, 40, 25, 55),
        "url":  f"{BASE_URL}/owtrad-gis-tmcXMEm1400.zip",
        "cite": "Ciolek 2004+, tmcXMEm1400",
    },
    {
        "id":   "tmcTRm1200a",
        "desc": "Turkey 1200 CE (part a)",
        "bbox": (25, 45, 35, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcTRm1200a.zip",
        "cite": "Ciolek 2004+, tmcTRm1200a",
    },
    {
        "id":   "tmcTRm1200b",
        "desc": "Turkey 1200 CE (part b)",
        "bbox": (25, 120, 35, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcTRm1200b.zip",
        "cite": "Ciolek 2004+, tmcTRm1200b",
    },
    {
        "id":   "tmcTRm1300",
        "desc": "Turkey / Anatolia 1300 CE",
        "bbox": (12, 48, 35, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcTRm1300.zip",
        "cite": "Ciolek 2004+, tmcTRm1300",
    },
    # ── Persia / Iran ──────────────────────────────────────────────────────────
    {
        "id":   "tmcIRa0100",
        "desc": "Iran/Persia 50 BCE–300 CE",
        "bbox": (30, 80, 10, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcIRa0100.zip",
        "cite": "Ciolek 2004+, tmcIRa0100",
    },
    {
        "id":   "tmcIRa0500",
        "desc": "Iran/Persia 300–700 CE",
        "bbox": (30, 80, 10, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcIRa0500.zip",
        "cite": "Ciolek 2004+, tmcIRa0500",
    },
    # ── Central Asia (Zone CA) ─────────────────────────────────────────────────
    {
        "id":   "tmcZCAm0200",
        "desc": "Iran & China 200 BCE–500 CE",
        "bbox": (25, 120, 15, 55),
        "url":  f"{BASE_URL}/owtrad-gis-tmcZCAm0200.zip",
        "cite": "Ciolek 2004+, tmcZCAm0200",
    },
    {
        "id":   "tmcZCAm0600",
        "desc": "Mediterranean–China 200 BCE–1400 CE",
        "bbox": (0, 120, 15, 55),
        "url":  f"{BASE_URL}/owtrad-gis-tmcZCAm0600.zip",
        "cite": "Ciolek 2004+, tmcZCAm0600",
    },
    {
        "id":   "tmcZCAm0800",
        "desc": "Central Asia 700–1000 CE",
        "bbox": (25, 120, 25, 55),
        "url":  f"{BASE_URL}/owtrad-gis-tmcZCAm0800.zip",
        "cite": "Ciolek 2004+, tmcZCAm0800",
    },
    {
        "id":   "tmcZCAm1000",
        "desc": "Central Asia 1–1400 CE",
        "bbox": (25, 135, 25, 50),
        "url":  f"{BASE_URL}/owtrad-gis-tmcZCAm1000.zip",
        "cite": "Ciolek 2004+, tmcZCAm1000",
    },
    {
        "id":   "tmcZCAm1000a",
        "desc": "Central Asia 1–1400 CE (variant a)",
        "bbox": (25, 135, 25, 50),
        "url":  f"{BASE_URL}/owtrad-gis-tmcZCAm1000a.zip",
        "cite": "Ciolek 2004+, tmcZCAm1000a",
    },
    # ── Turkmenistan ───────────────────────────────────────────────────────────
    {
        "id":   "tmcTMm1100",
        "desc": "Turkmenistan 1100 CE",
        "bbox": (50, 64, 35, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcTMm1100.zip",
        "cite": "Ciolek 2004+, tmcTMm1100",
    },
    # ── Tajikistan ─────────────────────────────────────────────────────────────
    {
        "id":   "tmcTJm0400a",
        "desc": "Tajikistan 400 CE (part a)",
        "bbox": (65, 82, 36, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcTJm0400a.zip",
        "cite": "Ciolek 2004+, tmcTJm0400a",
    },
    {
        "id":   "tmcTJm0400b",
        "desc": "Tajikistan 400 CE (part b)",
        "bbox": (65, 82, 36, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcTJm0400b.zip",
        "cite": "Ciolek 2004+, tmcTJm0400b",
    },
    # ── Kyrgyzstan (additional parts) ─────────────────────────────────────────
    {
        "id":   "tmcKGa0100a",
        "desc": "Kyrgyzstan 100 BCE–1400 CE (part a)",
        "bbox": (65, 82, 36, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcKGa0100a.zip",
        "cite": "Ciolek 2004+, tmcKGa0100a",
    },
    {
        "id":   "tmcKGa0100b",
        "desc": "Kyrgyzstan 100 BCE–1400 CE (part b)",
        "bbox": (65, 82, 36, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcKGa0100b.zip",
        "cite": "Ciolek 2004+, tmcKGa0100b",
    },
    {
        "id":   "tmcKGa0100c",
        "desc": "Kyrgyzstan 100 BCE–1400 CE (part c)",
        "bbox": (65, 82, 36, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcKGa0100c.zip",
        "cite": "Ciolek 2004+, tmcKGa0100c",
    },
    {
        "id":   "tmcKGa0100e",
        "desc": "Kyrgyzstan 100 BCE–1400 CE (part e)",
        "bbox": (65, 82, 36, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcKGa0100e.zip",
        "cite": "Ciolek 2004+, tmcKGa0100e",
    },
    # ── Middle East & India ────────────────────────────────────────────────────
    {
        "id":   "tmcZMEm1300",
        "desc": "Middle East & India 1300–1600 CE",
        "bbox": (25, 90, 5, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcZMEm1300.zip",
        "cite": "Ciolek 2004+, tmcZMEm1300",
    },
    {
        "id":   "tmcINm1550",
        "desc": "India 1550–1650 CE",
        "bbox": (65, 90, 8, 35),
        "url":  f"{BASE_URL}/owtrad-gis-tmcINm1550.zip",
        "cite": "Ciolek 2004+, tmcINm1550",
    },
    # ── China (all periods) ────────────────────────────────────────────────────
    {
        "id":   "tmcCNm1000",
        "desc": "NW China 100–1400 CE",
        "bbox": (74, 94, 36, 44),
        "url":  f"{BASE_URL}/owtrad-gis-tmcCNm1000.zip",
        "cite": "Ciolek 2004+, tmcCNm1000",
    },
    {
        "id":   "tmcCNm0620",
        "desc": "China 620 CE (Tang routes)",
        "bbox": (87, 104, 30, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcCNm0620.zip",
        "cite": "Ciolek 2004+, tmcCNm0620",
    },
    {
        "id":   "tmcCNm0680a",
        "desc": "China 680 CE (part a)",
        "bbox": (88, 104, 30, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcCNm0680a.zip",
        "cite": "Ciolek 2004+, tmcCNm0680a",
    },
    {
        "id":   "tmcCNm0680b",
        "desc": "China 680 CE (part b)",
        "bbox": (88, 104, 30, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcCNm0680b.zip",
        "cite": "Ciolek 2004+, tmcCNm0680b",
    },
    {
        "id":   "tmcCNm0680c",
        "desc": "China 680 CE (part c)",
        "bbox": (88, 104, 30, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcCNm0680c.zip",
        "cite": "Ciolek 2004+, tmcCNm0680c",
    },
    {
        "id":   "tmcCNm1500",
        "desc": "China Ming-era courier routes 1368–1644 CE",
        "bbox": (102, 121, 20, 50),
        "url":  f"{BASE_URL}/owtrad-gis-tmcCNm1500.zip",
        "cite": "Ciolek 2004+, tmcCNm1500",
    },
    {
        "id":   "tmcCNm1700",
        "desc": "China Qing-era routes 1644–1800 CE",
        "bbox": (115, 121, 25, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcCNm1700.zip",
        "cite": "Ciolek 2004+, tmcCNm1700",
    },
    {
        "id":   "tmcCNm1850",
        "desc": "China routes 1800–1900 CE",
        "bbox": (98, 117, 25, 45),
        "url":  f"{BASE_URL}/owtrad-gis-tmcCNm1850.zip",
        "cite": "Ciolek 2004+, tmcCNm1850",
    },
    {
        "id":   "tmcCNm1920",
        "desc": "China routes 1900–1949 CE",
        "bbox": (80, 117, 20, 50),
        "url":  f"{BASE_URL}/owtrad-gis-tmcCNm1920.zip",
        "cite": "Ciolek 2004+, tmcCNm1920",
    },
    # ── Southeast Asia ─────────────────────────────────────────────────────────
    {
        "id":   "tmcKHm1200",
        "desc": "Khmer/Cambodia 1200–1300 CE",
        "bbox": (100, 110, 10, 20),
        "url":  f"{BASE_URL}/owtrad-gis-tmcKHm1200.zip",
        "cite": "Ciolek 2004+, tmcKHm1200",
    },
]

ROUTE_FIELDNAMES = ["dataset", "lon1", "lat1", "lon2", "lat2",
                    "node1", "node2", "country1", "country2"]
NODE_FIELDNAMES  = ["dataset", "name", "lon", "lat", "country", "node_id"]


# ── Download ──────────────────────────────────────────────────────────────────

def fetch_zip(url: str, cache_path: Path, retries: int = 3) -> Optional[bytes]:
    """Download a zip file, returning None if the server returns 404."""
    if cache_path.exists():
        print(f"    [cached] {cache_path.name}")
        return cache_path.read_bytes()
    for attempt in range(retries):
        try:
            print(f"    Downloading {url} ...", flush=True)
            r = requests.get(url, headers={"User-Agent": USER_AGENT},
                             timeout=120, verify=False, stream=True)
            if r.status_code == 404:
                print(f"    [skip] 404 — no interchange zip for this dataset")
                return None
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
                print(f"    [skip] download failed after {retries} attempts: {e}")
                return None
    return None


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
    dataset_ids = ", ".join(ds["id"] for ds in DATASETS) + ", tmcKGa0100d"
    with open(ROUTES_CSV, "w", newline="", encoding="utf-8") as f:
        f.write("# owtrad_routes.csv — OWTRAD Silk Road route network (edge list)\n")
        f.write("# DO NOT HAND-EDIT. Regenerate: python3 data/scripts/fetch_owtrad_silk_road.py\n")
        f.write(f"# Generated: {fetch_date}\n")
        f.write("# Source: Ciolek, T.M. (2004+). OWTRAD Project. www.ciolek.com/owtrad.html\n")
        f.write(f"# Datasets attempted: {dataset_ids}\n")
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
        skipped = []
        for ds in DATASETS:
            print(f"\n  Dataset: {ds['id']} — {ds['desc']}")
            cache = CACHE_DIR / f"owtrad-gis-{ds['id']}.zip"
            data  = fetch_zip(ds["url"], cache)
            if data is None:
                skipped.append(ds["id"])
                continue
            r, n  = parse_zip(data, ds["id"])
            print(f"    Routes: {len(r)}  Nodes: {len(n)}")
            all_routes += r
            all_nodes  += n
        if skipped:
            print(f"\n  [info] No interchange zip found for: {', '.join(skipped)}")

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