#!/usr/bin/env python3
"""
fetch_p1435_global.py
=====================
Fetch ~100,000 heritage-designated buildings/structures from Wikidata
(P1435 = heritage designation, P625 = coordinates, P31 = built structure)
in regionally-balanced batches that match the UNESCO Cultural/Mixed
corpus geographic distribution.

WHY REGIONAL BATCHES
--------------------
Wikidata SPARQL times out on full-corpus scans.  Bounding-box filters
on regional subsets complete in < 30 s per page.  We fetch each region
in pages of PAGE_SIZE, sleep politely between requests, and write to a
CSV incrementally so the run can be interrupted and resumed.

STRUCTURE-TYPE FILTER
---------------------
We restrict to items whose P31 (instance of) is in a curated list of
built-structure types — building, temple, castle, mosque, church,
monument, stupa, fortress, archaeological site, etc.  This keeps the
sample comparable to UNESCO Cultural sites and avoids natural
features, political entities, and administrative regions.

TARGET SAMPLE (proportional to UNESCO Cultural corpus, N≈1,011)
---------------------------------------------------------------
  Europe          46%  →  46,000
  Asia-Pacific    18%  →  18,000
  Americas        14%  →  14,000
  Africa+MENA     16%  →  16,000
  Other            6%  →   6,000
  ─────────────────────────────
  Total                 100,000

USAGE
-----
  python3 data/fetch_p1435_global.py            # fetch all regions
  python3 data/fetch_p1435_global.py --region europe
  python3 data/fetch_p1435_global.py --force    # clear cache and re-fetch
  python3 data/fetch_p1435_global.py --dry-run  # print queries only

OUTPUT
------
  data/store/wikidata/p1435_global_control.csv — columns: wikidata_id, lat, lon, region
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import requests

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR   = Path(__file__).parent
OUTPUT_CSV = DATA_DIR / "p1435_global_control.csv"

# ── SPARQL config ─────────────────────────────────────────────────────────────
ENDPOINT    = "https://query.wikidata.org/sparql"
HEADERS     = {
    "User-Agent": (
        "gerizim-beru-control-study/2.0 "
        "(academic; beru-grid hypothesis; Python-requests)"
    ),
    "Accept": "application/sparql-results+json",
}
PAGE_SIZE   = 5_000   # rows per request — safe under Wikidata's 60s timeout
SLEEP_SEC   = 3.0     # polite inter-request delay
MAX_RETRIES = 3

# ── Built-structure P31 types ─────────────────────────────────────────────────
# Curated list of Wikidata Q-IDs for built structures comparable to
# UNESCO Cultural sites.  Avoids natural features, admin regions, people.
STRUCTURE_QIDS = [
    "Q41176",    # building (generic)
    "Q2977",     # cathedral
    "Q16560",    # palace
    "Q483453",   # temple
    "Q44613",    # monastery
    "Q839954",   # archaeological site
    "Q23413",    # castle
    "Q12280",    # bridge
    "Q179049",   # basilica
    "Q34627",    # synagogue
    "Q32815",    # mosque
    "Q108325",   # church building
    "Q4989906",  # monument
    "Q1261524",  # stupa
    "Q131596",   # fortress
    "Q57821",    # fortification
    "Q10529562", # ruins
    "Q1763828",  # amphitheatre
    "Q222139",   # aqueduct
    "Q1210064",  # ancient city
    "Q18758630", # ancient monument
    "Q39614",    # cemetery
    "Q24398318", # religious building
    "Q44494",    # mill
    "Q1081138",  # archaeological park
]

_VALUES = "VALUES ?stype { " + " ".join(f"wd:{q}" for q in STRUCTURE_QIDS) + " }"

# ── Regional bounding boxes ───────────────────────────────────────────────────
# (name, lat_min, lat_max, lon_min, lon_max, target_n)
REGIONS = [
    ("Europe",        35.0,  72.0, -12.0,  45.0, 46_000),
    ("Asia-Pacific", -10.0,  60.0,  60.0, 180.0, 18_000),
    ("Americas-N",   14.0,   75.0,-170.0, -50.0,  7_000),
    ("Americas-S",  -60.0,   14.0, -85.0, -30.0,  7_000),
    ("Africa",      -40.0,   38.0, -20.0,  55.0, 16_000),
    ("Other",        25.0,   45.0,  45.0,  65.0,  6_000),  # Central Asia / MENA
]

# ── SPARQL query ──────────────────────────────────────────────────────────────
_QUERY = """\
SELECT DISTINCT ?item (geof:latitude(?coord) AS ?lat) (geof:longitude(?coord) AS ?lon) WHERE {{
  {values}
  ?item wdt:P31   ?stype ;
        wdt:P1435 [] ;
        wdt:P625  ?coord .
  FILTER(geof:latitude(?coord)  >= {s} && geof:latitude(?coord)  <= {n})
  FILTER(geof:longitude(?coord) >= {w} && geof:longitude(?coord) <= {e})
}}
LIMIT {limit}
OFFSET {offset}
"""


def sparql_page(s, n, w, e, limit, offset) -> list[dict]:
    """Fetch one page; return list of {{wikidata_id, lat, lon}}."""
    query = _QUERY.format(values=_VALUES, s=s, n=n, w=w, e=e,
                          limit=limit, offset=offset)
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.get(ENDPOINT,
                             params={"query": query, "format": "json"},
                             headers=HEADERS, timeout=90)
            r.raise_for_status()
            rows = []
            for b in r.json()["results"]["bindings"]:
                try:
                    rows.append({
                        "wikidata_id": b["item"]["value"].rsplit("/", 1)[-1],
                        "lat": float(b["lat"]["value"]),
                        "lon": float(b["lon"]["value"]),
                    })
                except (KeyError, ValueError):
                    continue
            return rows
        except requests.exceptions.Timeout:
            print(f"      timeout (attempt {attempt}/{MAX_RETRIES})", file=sys.stderr)
            time.sleep(SLEEP_SEC * 4)
        except requests.exceptions.HTTPError as exc:
            code = exc.response.status_code
            print(f"      HTTP {code} (attempt {attempt}/{MAX_RETRIES})", file=sys.stderr)
            time.sleep(SLEEP_SEC * 5 if code == 429 else SLEEP_SEC * 2)
        except Exception as exc:
            print(f"      error: {exc} (attempt {attempt}/{MAX_RETRIES})", file=sys.stderr)
            time.sleep(SLEEP_SEC * 2)
    return []


def load_existing() -> tuple[list[dict], set[str]]:
    """Load the output CSV if it exists; return (rows, seen_qids)."""
    if not OUTPUT_CSV.exists():
        return [], set()
    rows, qids = [], set()
    with open(OUTPUT_CSV, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)
            qids.add(row["wikidata_id"])
    return rows, qids


def fetch_region(name, s, n, w, e, target, seen_qids, writer, fh,
                 dry_run=False) -> int:
    """
    Fetch up to `target` new items for one region, skipping QIDs already
    in `seen_qids`.  Writes directly to `writer`/`fh` for incremental saves.
    Returns count of new rows added.
    """
    print(f"\n  [{name}]  target={target:,}  "
          f"lat[{s},{n}]  lon[{w},{e}]")

    if dry_run:
        print(_QUERY.format(values="VALUES ?stype { ... }",
                            s=s, n=n, w=w, e=e,
                            limit=PAGE_SIZE, offset=0))
        return 0

    fetched = 0
    offset  = 0

    while fetched < target:
        page = sparql_page(s, n, w, e, PAGE_SIZE, offset)
        new  = 0
        for row in page:
            if row["wikidata_id"] not in seen_qids:
                seen_qids.add(row["wikidata_id"])
                writer.writerow({"wikidata_id": row["wikidata_id"],
                                 "lat": row["lat"], "lon": row["lon"],
                                 "region": name})
                new += 1
        fh.flush()
        fetched += new
        print(f"    offset={offset:>6}  page={len(page):>5}  "
              f"new={new:>5}  region_total={fetched:>6}")

        if len(page) < PAGE_SIZE:
            print(f"    region exhausted after {fetched:,} new sites")
            break
        if fetched >= target:
            break

        offset += PAGE_SIZE
        time.sleep(SLEEP_SEC)

    return fetched


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--region",  default=None,
                    choices=[r[0] for r in REGIONS],
                    help="Fetch only this region")
    ap.add_argument("--force",   action="store_true",
                    help="Delete cache and re-fetch from scratch")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print SPARQL queries without fetching")
    args = ap.parse_args()

    if args.force and OUTPUT_CSV.exists():
        OUTPUT_CSV.unlink()
        print("Cache cleared.")

    existing_rows, seen_qids = load_existing()
    print(f"Existing cache: {len(existing_rows):,} rows.")

    regions = [r for r in REGIONS if args.region is None or r[0] == args.region]

    mode = "a" if OUTPUT_CSV.exists() else "w"
    with open(OUTPUT_CSV, mode, newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["wikidata_id", "lat", "lon", "region"])
        if mode == "w":
            writer.writeheader()

        totals: dict[str, int] = {}
        for (name, s, n, w, e, target) in regions:
            added = fetch_region(name, s, n, w, e, target,
                                 seen_qids, writer, fh,
                                 dry_run=args.dry_run)
            totals[name] = added
            print(f"  → {name}: {added:,} new rows written")

    print("\n" + "=" * 50)
    print("  FETCH COMPLETE")
    print("=" * 50)
    for rname, cnt in totals.items():
        print(f"  {rname:<15}  {cnt:>7,}")
    print(f"  {'CACHE TOTAL':<15}  {len(seen_qids):>7,}")
    print(f"\n  Output → {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
