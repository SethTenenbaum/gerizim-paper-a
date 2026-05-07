"""
fetch_osm_stupas.py
===================
Fetches stupa features from OpenStreetMap via the Overpass API and writes a
deduplicated CSV to data/store/unesco/osm_stupas.csv.

Tags queried (nodes, ways, and relations):
  - man_made=stupa        (primary tag)
  - historic=stupa
  - building=stupa
  - ruins=stupa

For ways and relations, the centroid is computed from the bounding box
midpoint returned by Overpass (out center).

Deduplication: features within 0.01° (~1 km) of an already-accepted
feature are dropped (keeps the first encountered).

USAGE
-----
    cd /path/to/gerizim-paper-a
    python3 data/scripts/fetch_osm_stupas.py

OUTPUT
------
    data/store/unesco/osm_stupas.csv

LICENSE NOTE
------------
OpenStreetMap data © OpenStreetMap contributors, ODbL 1.0.
https://www.openstreetmap.org/copyright
"""

import csv
import json
import math
import time
import urllib.request
import urllib.parse
from pathlib import Path

OUT_PATH = Path(__file__).resolve().parents[2] / "data" / "store" / "unesco" / "osm_stupas.csv"
OVERPASS_URL = "https://overpass-api.de/api/interpreter"

# Tags to query — each generates its own Overpass union clause.
# man_made=stupa is intentionally excluded: it is applied in OSM to everything
# from garden ornaments to modern Western Buddhist centres and yields ~2,200
# largely unnamed, non-heritage features that are not comparable to the
# Wikidata Q180987 corpus.  building=stupa, historic=stupa, and ruins=stupa
# all require a mapped structural footprint or heritage designation, giving a
# corpus of ~240 features that is geographically and typologically comparable
# to the Wikidata dataset.
STUPA_TAGS = [
    ("building",  "stupa"),
    ("historic",  "stupa"),
    ("ruins",     "stupa"),
]

OVERPASS_QUERY = """
[out:json][timeout:120];
(
  node["historic"="stupa"];
  way["historic"="stupa"];
  relation["historic"="stupa"];
  node["building"="stupa"];
  way["building"="stupa"];
  relation["building"="stupa"];
  node["ruins"="stupa"];
  way["ruins"="stupa"];
  relation["ruins"="stupa"];
);
out center tags;
"""

DEDUP_THRESHOLD_DEG = 0.01   # ~1 km


def fetch_overpass(query: str, retries: int = 3) -> dict:
    data = urllib.parse.urlencode({"data": query}).encode()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                OVERPASS_URL,
                data=data,
                headers={"User-Agent": "gerizim-paper-a/1.0 (research; seth@fourthtemple.com)"},
            )
            with urllib.request.urlopen(req, timeout=180) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception as e:
            if attempt < retries - 1:
                print(f"  Attempt {attempt+1} failed ({e}), retrying in 10s...")
                time.sleep(10)
            else:
                raise


def centroid(element: dict):
    """Return (lat, lon) for a node, or bounding-box center for way/relation."""
    etype = element["type"]
    if etype == "node":
        return element.get("lat"), element.get("lon")
    # ways and relations: Overpass 'out center' puts midpoint in element["center"]
    center = element.get("center")
    if center:
        return center.get("lat"), center.get("lon")
    return None, None


def deduplicate(rows: list[dict], threshold: float) -> list[dict]:
    """Spatial deduplication: drop rows within `threshold` degrees of an accepted row."""
    kept = []
    for row in rows:
        lat, lon = float(row["lat"]), float(row["lon"])
        too_close = False
        for k in kept:
            dlat = lat - float(k["lat"])
            dlon = lon - float(k["lon"])
            if math.sqrt(dlat**2 + dlon**2) < threshold:
                too_close = True
                break
        if not too_close:
            kept.append(row)
    return kept


def main():
    print("Fetching OSM stupa features via Overpass API...")
    print(f"  Query tags: {[f'{k}={v}' for k,v in STUPA_TAGS]}")
    result = fetch_overpass(OVERPASS_QUERY)
    elements = result.get("elements", [])
    print(f"  Raw elements returned: {len(elements)}")

    rows = []
    skipped_no_coord = 0
    tag_counts = {}

    for el in elements:
        lat, lon = centroid(el)
        if lat is None or lon is None:
            skipped_no_coord += 1
            continue
        tags = el.get("tags", {})
        osm_id = f"{el['type']}/{el['id']}"
        name = tags.get("name") or tags.get("name:en") or ""
        country = tags.get("addr:country") or tags.get("is_in:country") or ""
        # Record which tag matched
        matched_tag = ""
        for k, v in STUPA_TAGS:
            if tags.get(k) == v:
                matched_tag = f"{k}={v}"
                tag_counts[matched_tag] = tag_counts.get(matched_tag, 0) + 1
                break
        rows.append({
            "osm_id":      osm_id,
            "name":        name,
            "lat":         f"{lat:.7f}",
            "lon":         f"{lon:.7f}",
            "country":     country,
            "matched_tag": matched_tag,
        })

    print(f"  With coordinates: {len(rows)}  (skipped no-coord: {skipped_no_coord})")
    print("  Tag breakdown:")
    for tag, count in sorted(tag_counts.items(), key=lambda x: -x[1]):
        print(f"    {tag}: {count}")

    rows_dedup = deduplicate(rows, DEDUP_THRESHOLD_DEG)
    print(f"  After deduplication (threshold {DEDUP_THRESHOLD_DEG}°): {len(rows_dedup)}")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", newline="", encoding="utf-8") as fh:
        fh.write("# osm_stupas.csv — OpenStreetMap stupa corpus\n")
        fh.write("# DO NOT HAND-EDIT. Regenerate with: python3 data/scripts/fetch_osm_stupas.py\n")
        fh.write(f"# Source: OpenStreetMap via Overpass API (overpass-api.de)\n")
        fh.write(f"# Tags: building=stupa, historic=stupa, ruins=stupa (man_made=stupa excluded — see script)\n")
        fh.write(f"# License: OpenStreetMap data (c) OpenStreetMap contributors, ODbL 1.0\n")
        fh.write(f"# N (deduplicated) = {len(rows_dedup)}\n")
        writer = csv.DictWriter(fh, fieldnames=["osm_id", "name", "lat", "lon", "country", "matched_tag"])
        writer.writeheader()
        writer.writerows(rows_dedup)

    print(f"\n  Written: {OUT_PATH}")
    print(f"  Final N = {len(rows_dedup)}")


if __name__ == "__main__":
    main()
