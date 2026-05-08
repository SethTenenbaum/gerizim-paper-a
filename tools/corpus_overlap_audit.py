#!/usr/bin/env python3
"""
corpus_overlap_audit.py — Estimate site-level overlap between the three
monument corpora: UNESCO dome/stupa subset, Wikidata Q180987, and OSM stupa.

Independence of these corpora rests primarily on curation process:
  - UNESCO: context-validated keyword match against WHC XML descriptions
  - Wikidata: SPARQL instance-of:Q180987 query, no researcher involvement
  - OSM:     crowd-sourced building/historic/ruins=stupa tags, Overpass API

This script provides a best-effort empirical check using two proxies:

  1. Name match  — case-insensitive partial string overlap between site
                   name fields (where available). OSM and Wikidata often
                   lack English names; this check is therefore incomplete
                   and will underestimate overlap when names are absent
                   or in different scripts.

  2. GPS proximity — pairs of sites (one from each corpus) with
                   geodesic distance < PROXIMITY_THRESHOLD_KM.

Results should be interpreted cautiously: GPS proximity below 1 km
suggests the same physical site, but corpora can independently represent
the same place without sharing records, and name absence in OSM/Wikidata
means name-match recall is low.

Usage:
    cd /path/to/gerizim-paper-a
    python3 tools/corpus_overlap_audit.py
"""

from __future__ import annotations
import sys
import csv
import math
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.dome_filter import is_dome_site

# ── Configuration ─────────────────────────────────────────────────────────────
PROXIMITY_THRESHOLD_KM = 1.0   # < 1 km  → likely the same physical site
EARTH_RADIUS_KM        = 6371.0

WIKIDATA_CSV = PROJECT_ROOT / "data" / "store" / "unesco" / "wikidata_stupas_q180987.csv"
OSM_CSV      = PROJECT_ROOT / "data" / "store" / "unesco" / "osm_stupas.csv"


# ── Helpers ───────────────────────────────────────────────────────────────────
def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Geodesic distance in km via Haversine formula."""
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi  = math.radians(lat2 - lat1)
    dlam  = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return 2 * EARTH_RADIUS_KM * math.asin(math.sqrt(a))


def name_overlap(a: str, b: str) -> bool:
    """Return True if names share a significant non-trivial token (≥ 5 chars)."""
    a_clean = a.lower().strip()
    b_clean = b.lower().strip()
    if not a_clean or not b_clean:
        return False
    # Exact or substring match
    if a_clean in b_clean or b_clean in a_clean:
        return True
    # Token-level match for tokens ≥ 5 characters
    a_tokens = {t for t in a_clean.split() if len(t) >= 5}
    b_tokens = {t for t in b_clean.split() if len(t) >= 5}
    return bool(a_tokens & b_tokens)


# ── Load corpora ──────────────────────────────────────────────────────────────
def load_unesco_dome() -> list[dict]:
    corpus  = load_corpus()
    sites   = cultural_sites_with_coords(corpus)
    return [
        {"name": s.site, "lat": getattr(s, "latitude", None), "lon": s.longitude}
        for s in sites
        if is_dome_site(s)
           and getattr(s, "latitude", None) is not None
    ]


def load_csv_corpus(path: Path) -> list[dict]:
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            try:
                rows.append({
                    "name": row.get("name", ""),
                    "lat":  float(row["lat"]),
                    "lon":  float(row["lon"]),
                })
            except (ValueError, KeyError):
                pass
    return rows


# ── Pairwise comparison ───────────────────────────────────────────────────────
def find_overlaps(
    corpus_a: list[dict],
    corpus_b: list[dict],
    label_a: str,
    label_b: str,
) -> dict:
    """Return overlap statistics between two corpora."""
    prox_pairs: list[tuple[str, str, float]] = []
    name_pairs: list[tuple[str, str, float]] = []

    for a in corpus_a:
        for b in corpus_b:
            dist = haversine_km(a["lat"], a["lon"], b["lat"], b["lon"])
            if dist < PROXIMITY_THRESHOLD_KM:
                prox_pairs.append((a["name"], b["name"], dist))
            if name_overlap(a["name"], b["name"]):
                name_pairs.append((a["name"], b["name"],
                                   haversine_km(a["lat"], a["lon"], b["lat"], b["lon"])))

    return {
        "label_a": label_a, "label_b": label_b,
        "n_a": len(corpus_a), "n_b": len(corpus_b),
        "proximity_pairs": prox_pairs,
        "name_pairs": name_pairs,
    }


def report(result: dict) -> None:
    la, lb = result["label_a"], result["label_b"]
    na, nb = result["n_a"], result["n_b"]
    pp = result["proximity_pairs"]
    np_ = result["name_pairs"]

    print(f"\n{'─'*60}")
    print(f"  {la}  (N={na})  ×  {lb}  (N={nb})")
    print(f"{'─'*60}")
    print(f"  GPS proximity < {PROXIMITY_THRESHOLD_KM} km:  {len(pp)} pair(s)")
    if pp:
        for a_name, b_name, dist in sorted(pp, key=lambda x: x[2])[:10]:
            print(f"    {dist:.3f} km  |  '{a_name[:40]}'  ×  '{b_name[:40]}'")
        if len(pp) > 10:
            print(f"    ... and {len(pp)-10} more")
    print(f"  Name overlap:           {len(np_)} pair(s)")
    if np_:
        for a_name, b_name, dist in sorted(np_, key=lambda x: x[2])[:10]:
            print(f"    {dist:.1f} km  |  '{a_name[:40]}'  ×  '{b_name[:40]}'")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Corpus overlap audit")
    print(f"  Proximity threshold: {PROXIMITY_THRESHOLD_KM} km")
    print("=" * 60)
    print()
    print("Loading corpora…")

    unesco  = load_unesco_dome()
    wiki    = load_csv_corpus(WIKIDATA_CSV)
    osm     = load_csv_corpus(OSM_CSV)

    print(f"  UNESCO dome/stupa:  N = {len(unesco)}")
    print(f"  Wikidata Q180987:   N = {len(wiki)}")
    print(f"  OSM stupa:          N = {len(osm)}")

    r1 = find_overlaps(unesco, wiki,  "UNESCO dome", "Wikidata Q180987")
    r2 = find_overlaps(unesco, osm,   "UNESCO dome", "OSM stupa")
    r3 = find_overlaps(wiki,   osm,   "Wikidata Q180987", "OSM stupa")

    report(r1)
    report(r2)
    report(r3)

    print()
    print("=" * 60)
    print("NOTE: Name-match recall is low when OSM/Wikidata records lack")
    print("English names. Independence of these corpora rests primarily on")
    print("their distinct curation processes, not on zero site overlap.")
    print("=" * 60)
