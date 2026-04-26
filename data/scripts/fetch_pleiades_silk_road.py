#!/usr/bin/env python3
"""
fetch_pleiades_silk_road.py
===========================
Download and filter the PLEIADES ancient-world gazetteer to produce a
reference set of ancient settlements in the Silk Road corridor.

SOURCE
------
PLEIADES is a peer-reviewed gazetteer of ancient places, published by the
Institute for the Study of the Ancient World (ISAW) at NYU.

  URL:     https://pleiades.stoa.org
  Data:    https://atlantides.org/downloads/pleiades/gis/pleiades_gis_data.zip
  Licence: CC BY 3.0  (https://creativecommons.org/licenses/by/3.0/)
  Cite:    Bagnall, R., Talbert, R., Bond, S., et al. (eds.).
           Pleiades: A community-built gazetteer and graph of ancient places.
           https://pleiades.stoa.org  DOI: 10.5281/zenodo.596503

FILTER CRITERIA (pre-specified, do not change post-hoc)
-------------------------------------------------------
  Geographic box : lon 28–122 °E, lat 18–56 °N
                   (Levant → East China; full overland Silk Road corridor)
  Feature types  : settlement, urban, settlement-modern, populated-place,
                   city, town, village (any key in SETTLEMENT_TYPES)
  Coordinates    : representative_longitude and representative_latitude present
  Precision      : any (rough-precision sites retained; flagged in output)

OUTPUT
------
  data/store/silk_road/pleiades_silk_road_cities.csv

  Columns: pleiades_id, title, lon, lat, place_types, location_precision

USAGE
-----
    python3 data/scripts/fetch_pleiades_silk_road.py           # download & filter
    python3 data/scripts/fetch_pleiades_silk_road.py --validate  # recheck stored file
"""

import argparse
import csv
import hashlib
import io
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("pip install requests  (required)")

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT       = Path(__file__).resolve().parents[2]
OUTPUT_DIR = ROOT / "data" / "store" / "silk_road"
OUTPUT_CSV = OUTPUT_DIR / "pleiades_silk_road_cities.csv"

# ── Source ───────────────────────────────��────────────────────────────────────
PLEIADES_ZIP_URL = "https://atlantides.org/downloads/pleiades/gis/pleiades_gis_data.zip"
USER_AGENT       = "gerizim-beru-grid/1.0 (academic research; pleiades-silk-road-filter)"

# ── Filter constants (pre-specified) ──────────────────────────────���──────────
LON_MIN, LON_MAX = 28.0, 122.0
LAT_MIN, LAT_MAX = 18.0,  56.0

SETTLEMENT_TYPES = {
    "settlement", "urban", "settlement-modern", "populated-place",
    "city", "town", "village", "municipium", "colonia", "polis",
    "emporion", "vicus", "pagus",
}

FIELDNAMES = ["pleiades_id", "title", "lon", "lat", "place_types", "location_precision"]


# ── Download ──────────────────────────────────────────────────────────────��───

def fetch_zip(retries: int = 3) -> bytes:
    for attempt in range(retries):
        try:
            print(f"  Downloading {PLEIADES_ZIP_URL} ...", flush=True)
            resp = requests.get(
                PLEIADES_ZIP_URL,
                headers={"User-Agent": USER_AGENT},
                timeout=180,
                stream=True,
            )
            resp.raise_for_status()
            chunks = []
            total = 0
            for chunk in resp.iter_content(65536):
                chunks.append(chunk)
                total += len(chunk)
                print(f"\r  Downloaded {total/1024/1024:.1f} MB", end="", flush=True)
            print()
            return b"".join(chunks)
        except requests.exceptions.Timeout:
            if attempt < retries - 1:
                print(f"  [timeout, retry {attempt+1}/{retries}]")
                time.sleep(15)
            else:
                raise
        except Exception as e:
            if attempt < retries - 1:
                print(f"  [error: {e}, retry {attempt+1}/{retries}]")
                time.sleep(10)
            else:
                raise
    return b""


# ── Parse & filter ────────────────────────────────────────────────────────────

def _read_csv_from_zip(zf, path):
    with zf.open(path) as f:
        text = io.TextIOWrapper(f, encoding="utf-8-sig")
        yield from csv.DictReader(text)


def build_type_index(zf) -> dict:
    """Build {place_id: [place_type, ...]} from places_place_types.csv."""
    idx = {}
    for row in _read_csv_from_zip(zf, "data/gis/places_place_types.csv"):
        pid = row["place_id"].strip()
        pt  = row["place_type"].strip().lower()
        idx.setdefault(pid, []).append(pt)
    return idx


def filter_places(zf, type_idx: dict) -> list:
    kept = []
    total = 0
    for row in _read_csv_from_zip(zf, "data/gis/places.csv"):
        total += 1
        try:
            lon = float(row["representative_longitude"])
            lat = float(row["representative_latitude"])
        except (ValueError, KeyError):
            continue
        if not (LON_MIN <= lon <= LON_MAX and LAT_MIN <= lat <= LAT_MAX):
            continue
        pid   = row["id"].strip()
        types = type_idx.get(pid, [])
        if not any(t in SETTLEMENT_TYPES for t in types):
            continue
        kept.append({
            "pleiades_id":       pid,
            "title":             row.get("title", ""),
            "lon":               f"{lon:.6f}",
            "lat":               f"{lat:.6f}",
            "place_types":       ";".join(types),
            "location_precision": row.get("location_precision", ""),
        })
    print(f"  Total PLEIADES places: {total}")
    print(f"  Passing all filters:   {len(kept)}")
    return kept


# ── Write CSV ─────────────────────���───────────────────────────────────────────

def write_csv(rows: list, fetch_date: str) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        f.write("# pleiades_silk_road_cities.csv — PLEIADES ancient settlements, Silk Road corridor\n")
        f.write("# DO NOT HAND-EDIT. Regenerate: python3 data/scripts/fetch_pleiades_silk_road.py\n")
        f.write(f"# Generated: {fetch_date}\n")
        f.write(f"# Source: {PLEIADES_ZIP_URL}\n")
        f.write(f"# Cite: Bagnall et al. (eds.), Pleiades. https://pleiades.stoa.org  DOI:10.5281/zenodo.596503\n")
        f.write(f"# Licence: CC BY 3.0\n")
        f.write(f"# Filter: lon {LON_MIN}–{LON_MAX}E, lat {LAT_MIN}–{LAT_MAX}N; "
                f"types=settlement/urban variants\n")
        f.write(f"# N = {len(rows)}\n")
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        w.writeheader()
        w.writerows(rows)


# ── Validate mode ─────────────────────────────────────────────────────────────

def load_stored_csv() -> list:
    rows = []
    with open(OUTPUT_CSV, newline="", encoding="utf-8") as f:
        non_comment = (line for line in f if not line.startswith("#"))
        for r in csv.DictReader(non_comment):
            rows.append(r)
    return rows


def run_validate() -> int:
    if not OUTPUT_CSV.exists():
        print(f"  ERROR: {OUTPUT_CSV} does not exist. Run without --validate to fetch.")
        return 1
    rows = load_stored_csv()
    sha  = hashlib.sha256(OUTPUT_CSV.read_bytes()).hexdigest()
    lons = [float(r["lon"]) for r in rows]
    lats = [float(r["lat"]) for r in rows]
    print(f"  File:    {OUTPUT_CSV.relative_to(ROOT)}")
    print(f"  Rows:    {len(rows)}")
    print(f"  SHA-256: {sha}")
    print(f"  Lon:     {min(lons):.2f} – {max(lons):.2f}°E")
    print(f"  Lat:     {min(lats):.2f} – {max(lats):.2f}°N")
    rough = sum(1 for r in rows if r.get("location_precision") == "rough")
    print(f"  Rough precision: {rough} / {len(rows)}")
    print(f"  Sample (first 8):")
    for r in rows[:8]:
        print(f"    {r['pleiades_id']:>7}  {r['title'][:38]:38s}  {float(r['lon']):7.2f}E  {float(r['lat']):6.2f}N")
    return 0


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Fetch PLEIADES gazetteer and filter for Silk Road corridor settlements."
    )
    parser.add_argument("--validate", action="store_true",
                        help="Check the stored CSV without re-downloading.")
    args = parser.parse_args()

    if args.validate:
        sys.exit(run_validate())

    fetch_date = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    print("=" * 72)
    print("  PLEIADES SILK ROAD SETTLEMENT FILTER")
    print(f"  Run at: {fetch_date}")
    print(f"  Source: {PLEIADES_ZIP_URL}")
    print("=" * 72)

    raw = fetch_zip()
    zf  = zipfile.ZipFile(io.BytesIO(raw))

    print("\n  Building place-type index ...")
    type_idx = build_type_index(zf)
    print(f"  Place-type entries: {len(type_idx)}")

    print("\n  Filtering places ...")
    rows = filter_places(zf, type_idx)

    write_csv(rows, fetch_date)

    sha = hashlib.sha256(OUTPUT_CSV.read_bytes()).hexdigest()
    print(f"\n  Written: {OUTPUT_CSV.relative_to(ROOT)}")
    print(f"  Rows:    {len(rows)}")
    print(f"  SHA-256: {sha}")
    print("\n  Done.")


if __name__ == "__main__":
    main()