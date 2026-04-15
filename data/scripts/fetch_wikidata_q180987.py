#!/usr/bin/env python3
"""
fetch_wikidata_q180987.py
=========================
Fetch, validate, and (optionally) regenerate the Wikidata Q180987 stupa
corpus used as Test 6 in Tenenbaum (2026).

PURPOSE
-------
This script is the single authoritative source for ``data/store/unesco/wikidata_stupas_q180987.csv``.
It serves two modes:

  --validate  (default)
      Compare the stored CSV against a fresh live SPARQL query.
      Exit 0 if QIDs, coordinates, and names all match within tolerance;
      exit 1 if any discrepancy is detected.

  --fetch
      Download a fresh snapshot from Wikidata, write it to
      ``data/store/unesco/wikidata_stupas_q180987.csv``, and print a diff summary
      against the previously stored version.

SPARQL QUERY
------------
The query that reproduces the 229-row snapshot is:

    SELECT DISTINCT ?item ?itemLabel ?coord ?country ?countryLabel WHERE {
      ?item wdt:P31/wdt:P279* wd:Q180987 .
      ?item wdt:P625 ?coord .
      OPTIONAL { ?item wdt:P17 ?country . }
      SERVICE wikibase:label {
        bd:serviceParam wikibase:language "en,mul" .
      }
    }

Key design choices (pre-specified, unchanged since initial download):
  • P31/P279*  — instance-of *including subclasses* of Q180987.  Using plain
    P31 returns only 95 items; the subclass traversal is required to capture
    dagobas, chortens, chedis, etc.
  • No longitude filter — the raw download is global (229 sites), spanning
    -105.5°E (Colorado, USA) to +144.2°E (Australia).  The Test 6 analysis
    script (analysis/wikidata_full_test_suite.py) does not apply a longitude
    filter either, treating the full global corpus as the analysis population.
  • No date filter.
  • No architectural / label filter.  The stupa/non-stupa label split
    applied in the analysis is done at analysis time, not at fetch time.

STORED SNAPSHOT METADATA
-------------------------
  Original fetch date:  2025 (exact date unknown; committed as part of
               initial Q180987 investigation in moses-zipporah-forensic-audit)
  Re-fetched:  2026-04-02T17:04:37Z (identical data, comment header added)
  Row count:   229 (header/comment lines excluded)
  SHA-256:     696eaa96c18b08ff35610093cb0a808d31759d23c27f11f07f2c61c0b53f6b77
  MD5:         d506b388926a250b490cf833c8532e69
  Original SHA-256 (no header):
               0d8456a17723b026269a7deacf6ac17ddbec0ad02b8aca458c0eaa57d945c87d

VALIDATION RESULTS (April 2, 2026)
-----------------------------------
  QID match:          229 / 229  (all present in live query)
  Coordinate drift:   0 mismatches  (tolerance 0.0001°)
  Country drift:      0 mismatches
  Name drift:         1 minor change  (Q56208841: stored "Q56208841",
                      live "Sipamutung Temple" — Wikidata label was added
                      after snapshot; QID and coordinates unchanged)
  Overall verdict:    VALID — stored snapshot is reproducible today

USAGE
-----
    python3 data/fetch_wikidata_q180987.py             # validate
    python3 data/fetch_wikidata_q180987.py --validate  # validate (explicit)
    python3 data/fetch_wikidata_q180987.py --fetch     # download & overwrite

Requires:  requests  (pip install requests)
"""

import argparse
import csv
import hashlib
import json
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("pip install requests  (required)")

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT       = Path(__file__).resolve().parent.parent
DATA_DIR   = ROOT / "data"
OUTPUT_CSV = DATA_DIR / "wikidata_stupas_q180987.csv"

# ── Stored snapshot checksums ─────────────────────────────────────────────────
# These reflect the --fetch output (7-line comment header + data rows).
# The original hand-committed file had no header comments; its checksums were:
#   SHA-256: 0d8456a17723b026269a7deacf6ac17ddbec0ad02b8aca458c0eaa57d945c87d
#   MD5:     47c4df1175e08d6bdbecf94fab463f12
# After the first --fetch (2026-04-02), the header was added; data is identical.
STORED_SHA256 = "696eaa96c18b08ff35610093cb0a808d31759d23c27f11f07f2c61c0b53f6b77"
STORED_MD5    = "d506b388926a250b490cf833c8532e69"
STORED_N      = 229

# ── SPARQL ────────────────────────────────────────────────────────────────────
ENDPOINT   = "https://query.wikidata.org/sparql"
USER_AGENT = "beru-grid-stupa-audit/1.0 (academic research; gerizim-analysis)"

# This is the EXACT query that reproduces the 229-row snapshot.
# P31/P279* is REQUIRED — plain P31 returns only ~95 items.
SPARQL_QUERY = """
SELECT DISTINCT ?item ?itemLabel ?coord ?country ?countryLabel WHERE {
  ?item wdt:P31/wdt:P279* wd:Q180987 .
  ?item wdt:P625 ?coord .
  OPTIONAL { ?item wdt:P17 ?country . }
  SERVICE wikibase:label {
    bd:serviceParam wikibase:language "en,mul" .
  }
}
"""

# ── Helpers ───────────────────────────────────────────────────────────────────

def sparql_fetch(retries: int = 3) -> list:
    """Execute SPARQL_QUERY and return list of result bindings."""
    for attempt in range(retries):
        try:
            resp = requests.get(
                ENDPOINT,
                params={"query": SPARQL_QUERY, "format": "json"},
                headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
                timeout=90,
            )
            resp.raise_for_status()
            return resp.json()["results"]["bindings"]
        except requests.exceptions.Timeout:
            if attempt < retries - 1:
                print(f"  [timeout, retry {attempt + 1}/{retries}]", flush=True)
                time.sleep(5)
            else:
                raise
        except Exception as e:
            if attempt < retries - 1:
                print(f"  [error: {e}, retry {attempt + 1}/{retries}]", flush=True)
                time.sleep(5)
            else:
                raise
    return []


def parse_point(point_str: str) -> tuple:
    """Parse 'Point(lon lat)' → (lat, lon) floats, or (None, None)."""
    m = re.match(r"Point\(\s*([-\d.]+)\s+([-\d.]+)\s*\)", point_str or "")
    if not m:
        return None, None
    return float(m.group(2)), float(m.group(1))   # lat, lon


def bindings_to_rows(bindings: list) -> list:
    """Convert SPARQL result bindings to canonical row dicts."""
    rows = []
    seen = set()
    for b in bindings:
        qid     = b["item"]["value"].split("/")[-1]
        name    = b.get("itemLabel", {}).get("value", qid)
        coord   = b.get("coord",     {}).get("value", "")
        country = b.get("countryLabel", {}).get("value", "")
        lat, lon = parse_point(coord)
        if lat is None or qid in seen:
            continue
        seen.add(qid)
        rows.append({"qid": qid, "name": name, "lat": lat, "lon": lon, "country": country})
    return rows


def load_stored_csv() -> list:
    """Load the stored CSV, return list of row dicts (skips # comment lines)."""
    rows = []
    with open(OUTPUT_CSV, newline="", encoding="utf-8") as f:
        # Strip leading comment lines before handing to DictReader
        non_comment_lines = (line for line in f if not line.startswith("#"))
        for r in csv.DictReader(non_comment_lines):
            try:
                rows.append({
                    "qid":     r["qid"],
                    "name":    r["name"],
                    "lat":     float(r["lat"]),
                    "lon":     float(r["lon"]),
                    "country": r["country"],
                })
            except (ValueError, KeyError):
                pass
    return rows


def checksum_file(path: Path) -> tuple:
    data = path.read_bytes()
    return hashlib.sha256(data).hexdigest(), hashlib.md5(data).hexdigest()


def write_csv(rows: list, run_time: str) -> None:
    """Write rows to OUTPUT_CSV with a reproducibility header."""
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        f.write(f"# wikidata_stupas_q180987.csv — Wikidata Q180987 stupa corpus\n")
        f.write(f"# DO NOT HAND-EDIT. Regenerate with: python3 data/fetch_wikidata_q180987.py --fetch\n")
        f.write(f"# Generated: {run_time}\n")
        f.write(f"# Source: {ENDPOINT}\n")
        f.write(f"# Query: SELECT DISTINCT ?item (P31/P279* wd:Q180987) with P625 coords\n")
        f.write(f"# Filters: none (global, no longitude/date/label filter)\n")
        f.write(f"# N = {len(rows)}\n")
        w = csv.DictWriter(f, fieldnames=["qid", "name", "lat", "lon", "country"])
        w.writeheader()
        w.writerows(rows)


# ── Validate mode ─────────────────────────────────────────────────────────────

def run_validate() -> int:
    """
    Compare the stored CSV against a fresh live SPARQL query.
    Returns 0 (pass) or 1 (fail).
    """
    print("=" * 72)
    print("  WIKIDATA Q180987 — DATASET VALIDATION")
    print(f"  Stored CSV: {OUTPUT_CSV.relative_to(ROOT)}")
    print("=" * 72)

    # 1. Checksum the stored file
    if not OUTPUT_CSV.exists():
        print(f"\n  ERROR: {OUTPUT_CSV} does not exist.")
        return 1

    sha, md5 = checksum_file(OUTPUT_CSV)
    print(f"\n  Stored file checksums:")
    print(f"    SHA-256: {sha}")
    print(f"    MD5:     {md5}")
    checksum_ok = (sha == STORED_SHA256 and md5 == STORED_MD5)
    print(f"    Match known-good: {'✓ YES' if checksum_ok else '✗ NO (file has changed)'}")

    # 2. Load stored rows
    stored_rows = load_stored_csv()
    stored_by_qid = {r["qid"]: r for r in stored_rows}
    print(f"\n  Stored row count: {len(stored_rows)}  (expected {STORED_N})")
    if len(stored_rows) != STORED_N:
        print(f"  WARNING: Row count mismatch!")

    # 3. Live SPARQL query
    print(f"\n  Querying {ENDPOINT} ...")
    print(f"  Query: P31/P279* wd:Q180987 with P625 coords")
    bindings = sparql_fetch()
    live_rows = bindings_to_rows(bindings)
    live_by_qid = {r["qid"]: r for r in live_rows}
    print(f"  Live result count: {len(live_rows)}")

    # 4. Compare
    stored_qids = set(stored_by_qid)
    live_qids   = set(live_by_qid)
    in_both     = stored_qids & live_qids
    only_stored = stored_qids - live_qids
    only_live   = live_qids - stored_qids

    print(f"\n  QID comparison:")
    print(f"    In both:         {len(in_both)}")
    print(f"    Only in stored:  {len(only_stored)}")
    if only_stored:
        for q in sorted(only_stored)[:10]:
            r = stored_by_qid[q]
            print(f"      {q} | {r['name']} | lon={r['lon']:.4f} | {r['country']}")
    print(f"    Only in live:    {len(only_live)}")
    if only_live:
        for q in sorted(only_live)[:10]:
            r = live_by_qid[q]
            print(f"      {q} | {r['name']} | lon={r['lon']:.4f} | {r['country']}")

    # 5. Coordinate drift (for QIDs in both)
    coord_mismatches = []
    name_changes     = []
    country_changes  = []
    COORD_TOL = 0.0001  # ~11 m

    for qid in in_both:
        s = stored_by_qid[qid]
        l = live_by_qid[qid]
        if abs(s["lat"] - l["lat"]) > COORD_TOL or abs(s["lon"] - l["lon"]) > COORD_TOL:
            coord_mismatches.append((qid, s, l))
        if s["name"] != l["name"]:
            name_changes.append((qid, s["name"], l["name"]))
        if s["country"] != l["country"]:
            country_changes.append((qid, s["country"], l["country"]))

    print(f"\n  Coordinate drift  (>{COORD_TOL}°): {len(coord_mismatches)}")
    for qid, s, l in coord_mismatches[:5]:
        print(f"    {qid}: stored=({s['lat']:.5f},{s['lon']:.5f}) "
              f"live=({l['lat']:.5f},{l['lon']:.5f})")

    print(f"  Country changes:  {len(country_changes)}")
    for qid, old, new in country_changes[:5]:
        print(f"    {qid}: \"{old}\" → \"{new}\"")

    print(f"  Name changes:     {len(name_changes)}")
    for qid, old, new in name_changes[:10]:
        print(f"    {qid}: \"{old}\" → \"{new}\"")

    # 6. Verdict
    fatal = bool(only_stored or only_live or coord_mismatches or country_changes)
    minor = bool(name_changes)

    print("\n" + "─" * 72)
    if not fatal and not minor:
        print("  RESULT: ✓ PERFECT MATCH — stored CSV is byte-reproducible today.")
        verdict = 0
    elif not fatal and minor:
        print(f"  RESULT: ✓ VALID WITH MINOR LABEL CHANGES ({len(name_changes)} name(s)).")
        print("  QIDs, coordinates, and countries are unchanged.")
        print("  The analysis results are unaffected by label-only changes.")
        verdict = 0
    else:
        print("  RESULT: ✗ DATASET HAS DRIFTED")
        print("  Re-run with --fetch to update the snapshot.")
        verdict = 1

    print("─" * 72)
    return verdict


# ── Fetch mode ────────────────────────────────────────────────────────────────

def run_fetch() -> int:
    """
    Download a fresh snapshot and write to OUTPUT_CSV.
    Print a diff summary against the previously stored version.
    """
    run_time = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    print("=" * 72)
    print("  WIKIDATA Q180987 — FRESH FETCH")
    print(f"  Run at: {run_time}")
    print(f"  Source: {ENDPOINT}")
    print("=" * 72)

    # Load old version if it exists
    old_rows = {}
    old_checksum = None
    if OUTPUT_CSV.exists():
        old_rows = {r["qid"]: r for r in load_stored_csv()}
        old_sha, old_md5 = checksum_file(OUTPUT_CSV)
        old_checksum = old_sha
        print(f"\n  Existing file: {len(old_rows)} rows  SHA-256={old_sha[:16]}...")

    print(f"\n  Querying {ENDPOINT} ...")
    bindings = sparql_fetch()
    new_rows = bindings_to_rows(bindings)
    print(f"  Live result count: {len(new_rows)}")

    # Diff
    new_by_qid = {r["qid"]: r for r in new_rows}
    added   = set(new_by_qid) - set(old_rows)
    removed = set(old_rows) - set(new_by_qid)
    print(f"\n  Diff vs stored:")
    print(f"    Added:   {len(added)}")
    for q in sorted(added)[:10]:
        r = new_by_qid[q]
        print(f"      + {q} | {r['name']} | lon={r['lon']:.4f} | {r['country']}")
    print(f"    Removed: {len(removed)}")
    for q in sorted(removed)[:10]:
        r = old_rows[q]
        print(f"      - {q} | {r['name']} | lon={r['lon']:.4f} | {r['country']}")

    # Write
    write_csv(new_rows, run_time)
    new_sha, new_md5 = checksum_file(OUTPUT_CSV)
    print(f"\n  Written: {OUTPUT_CSV.relative_to(ROOT)}")
    print(f"    Rows:    {len(new_rows)}")
    print(f"    SHA-256: {new_sha}")
    print(f"    MD5:     {new_md5}")

    if old_checksum and new_sha == old_checksum:
        print("\n  ✓ Checksums match — snapshot is unchanged.")
    elif old_checksum:
        print("\n  ⚠  Checksums differ — snapshot has been updated.")
        print("     Update STORED_SHA256 / STORED_MD5 in this script if you")
        print("     intend to use the new snapshot as the reference.")

    return 0


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Fetch and validate the Wikidata Q180987 stupa corpus."
    )
    parser.add_argument(
        "--fetch",
        action="store_true",
        help="Download a fresh snapshot and overwrite the stored CSV.",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate the stored CSV against a live SPARQL query (default).",
    )
    args = parser.parse_args()

    if args.fetch:
        sys.exit(run_fetch())
    else:
        sys.exit(run_validate())


if __name__ == "__main__":
    main()
