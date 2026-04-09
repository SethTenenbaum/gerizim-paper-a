#!/usr/bin/env python3
"""
fetch_extended.py — Download extended descriptions from UNESCO web pages.
=========================================================================

For each site in `data/store/unesco/unesco.xml`, fetches the full "Outstanding Universal
Value" text from its UNESCO web page (http_url) and caches it in
`data/store/unesco/extended_cache.json`.

The extended text contains:
  - Brief Synthesis
  - Criterion-by-criterion justification
  - Integrity statement
  - Authenticity statement
  - Protection and management requirements

This text is MUCH richer than the XML `short_description`. For example,
the word "stupa" appears in the Kathmandu Valley extended text but NOT in
its XML short_description.

USAGE
─────
  python3 data/fetch_extended.py              # fetch all (skip cached)
  python3 data/fetch_extended.py --force      # re-fetch everything
  python3 data/fetch_extended.py --limit 50   # fetch only 50 sites
  python3 data/fetch_extended.py --id 121     # fetch one specific site

RATE LIMITING
─────────────
  UNESCO's server is respected with a 1-second delay between requests.
  The full corpus (~1,248 sites) takes ~25 minutes.

OUTPUT
──────
  data/store/unesco/extended_cache.json — dict keyed by id_number:
    {
      "121": {
        "id_number": "121",
        "site": "Kathmandu Valley",
        "http_url": "https://whc.unesco.org/en/list/121",
        "full_text": "Located in the foothills of the Himalayas...",
        "brief_synthesis": "Located in the foothills...",
        "criteria_detail": "Criterion (iii): The seven monument...",
        "fetched_at": "2025-01-15T12:34:56"
      },
      ...
    }
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

import cloudscraper
from bs4 import BeautifulSoup

# ── Paths ─────────────────────────────────────────────────────────────────────
_DATA_DIR      = Path(__file__).parent
XML_PATH       = _DATA_DIR / "unesco.xml"
EXTENDED_CACHE = _DATA_DIR / "extended_cache.json"

# ── HTTP config ───────────────────────────────────────────────────────────────
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/120.0.0.0 Safari/537.36"
    ),
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.9",
    "Accept-Encoding": "gzip, deflate, br",
    "Connection": "keep-alive",
    "Upgrade-Insecure-Requests": "1",
}
REQUEST_TIMEOUT = 30
DELAY_SECONDS   = 1.5   # polite delay between requests
MAX_RETRIES     = 2      # retry failed requests once


# ── Parsing helpers ───────────────────────────────────────────────────────────
def extract_ouv_text(html: str) -> Dict[str, str]:
    """
    Extract the "Outstanding Universal Value" section from a UNESCO page.

    Returns a dict with:
      - full_text:        The entire OUV section as plain text
      - brief_synthesis:  Just the Brief Synthesis paragraph(s)
      - criteria_detail:  The Criterion (i), (ii), ... paragraphs
    """
    soup = BeautifulSoup(html, "html.parser")

    result = {
        "full_text": "",
        "brief_synthesis": "",
        "criteria_detail": "",
    }

    # The OUV is typically under an <h2> "Outstanding Universal Value"
    # or within the main description area
    full_parts = []
    brief_parts = []
    criteria_parts = []

    # Strategy 1: Find the OUV heading and grab all following text
    ouv_heading = None
    for h2 in soup.find_all("h2"):
        if "outstanding universal value" in h2.get_text().lower():
            ouv_heading = h2
            break

    if ouv_heading:
        # Walk through siblings after the OUV heading
        in_brief = False
        in_criteria = False

        for sibling in ouv_heading.find_next_siblings():
            # Stop at next h2
            if sibling.name == "h2":
                break

            text = sibling.get_text(separator=" ", strip=True)
            if not text:
                continue

            full_parts.append(text)

            # Track sections
            text_lower = text.lower()
            if "brief synthesis" in text_lower or "brief description" in text_lower:
                in_brief = True
                in_criteria = False
            elif re.search(r"criterion\s+\(", text_lower):
                in_brief = False
                in_criteria = True
                criteria_parts.append(text)
            elif any(kw in text_lower for kw in [
                "integrity", "authenticity", "protection and management"
            ]):
                in_brief = False
                in_criteria = False
            elif in_brief:
                brief_parts.append(text)
            elif in_criteria:
                criteria_parts.append(text)

    # Strategy 2: Fallback — grab all paragraph text from the page body
    # (catches pages with different HTML structure)
    if not full_parts:
        # Look for the main content area
        main = soup.find("article") or soup.find("main") or soup.find("div", class_="main")
        if main is None:
            main = soup

        for p in main.find_all("p"):
            text = p.get_text(separator=" ", strip=True)
            if len(text) > 50:  # skip trivial paragraphs
                full_parts.append(text)

    result["full_text"] = "\n\n".join(full_parts)
    result["brief_synthesis"] = "\n\n".join(brief_parts)
    result["criteria_detail"] = "\n\n".join(criteria_parts)

    return result


def fetch_site_extended(url: str) -> Optional[Dict[str, str]]:
    """Fetch a single UNESCO page and extract extended text."""
    scraper = cloudscraper.create_scraper()
    for attempt in range(MAX_RETRIES + 1):
        try:
            resp = scraper.get(url, timeout=REQUEST_TIMEOUT)
            resp.raise_for_status()
            return extract_ouv_text(resp.text)
        except Exception as e:
            if attempt < MAX_RETRIES:
                time.sleep(2 * (attempt + 1))
                continue
            print(f"  ✗ Failed: {e}")
            return None


# ── Main ──────────────────────────────────────────────────────────────────────
def load_sites_from_xml() -> list:
    """Parse XML to get basic site info (id, name, url)."""
    tree = ET.parse(str(XML_PATH))
    root = tree.getroot()
    sites = []
    for row in root.findall("row"):
        def get(tag):
            el = row.find(tag)
            return (el.text or "").strip() if el is not None and el.text else ""

        sites.append({
            "id_number": get("id_number"),
            "site": get("site"),
            "http_url": get("http_url"),
        })
    return sites


def load_cache() -> Dict[str, dict]:
    """Load existing cache."""
    if EXTENDED_CACHE.exists():
        with open(EXTENDED_CACHE, "r", encoding="utf-8") as f:
            return json.load(f)
    return {}


def save_cache(cache: Dict[str, dict]) -> None:
    """Save cache to disk."""
    with open(EXTENDED_CACHE, "w", encoding="utf-8") as f:
        json.dump(cache, f, ensure_ascii=False, indent=2)


def main():
    parser = argparse.ArgumentParser(description="Fetch extended UNESCO descriptions")
    parser.add_argument("--force", action="store_true",
                        help="Re-fetch all sites, even if cached")
    parser.add_argument("--limit", type=int, default=0,
                        help="Maximum number of sites to fetch (0 = all)")
    parser.add_argument("--id", type=str, default="",
                        help="Fetch only a specific site by id_number")
    args = parser.parse_args()

    sites = load_sites_from_xml()
    cache = load_cache()

    # Deduplicate by id_number (transnational entries may share IDs)
    seen = set()
    unique_sites = []
    for s in sites:
        if s["id_number"] not in seen:
            seen.add(s["id_number"])
            unique_sites.append(s)

    # Filter to specific site if requested
    if args.id:
        unique_sites = [s for s in unique_sites if s["id_number"] == args.id]
        if not unique_sites:
            print(f"No site found with id_number={args.id}")
            sys.exit(1)

    # Determine which sites need fetching
    to_fetch = []
    for s in unique_sites:
        if args.force or s["id_number"] not in cache:
            to_fetch.append(s)

    if not to_fetch:
        print(f"All {len(unique_sites)} sites already cached. "
              f"Use --force to re-fetch.")
        print(f"Cache: {EXTENDED_CACHE}")
        return

    if args.limit > 0:
        to_fetch = to_fetch[:args.limit]

    print(f"UNESCO Extended Description Fetcher")
    print(f"───────────────────────────────────")
    print(f"Total unique sites:   {len(unique_sites)}")
    print(f"Already cached:       {len(cache)}")
    print(f"To fetch:             {len(to_fetch)}")
    print(f"Delay between reqs:   {DELAY_SECONDS}s")
    print(f"Estimated time:       ~{len(to_fetch) * (DELAY_SECONDS + 1):.0f}s")
    print()

    fetched = 0
    failed = 0

    for i, s in enumerate(to_fetch, 1):
        url = s["http_url"]
        name = s["site"]
        id_num = s["id_number"]

        print(f"[{i}/{len(to_fetch)}] {id_num}: {name[:60]}...", end=" ", flush=True)

        result = fetch_site_extended(url)
        if result:
            cache[id_num] = {
                "id_number": id_num,
                "site": name,
                "http_url": url,
                **result,
                "fetched_at": datetime.utcnow().isoformat(),
            }
            text_len = len(result.get("full_text", ""))
            print(f"✓ ({text_len:,} chars)")
            fetched += 1
        else:
            failed += 1

        # Save periodically (every 25 sites)
        if i % 25 == 0:
            save_cache(cache)
            print(f"  [saved {len(cache)} entries to cache]")

        # Rate limit
        if i < len(to_fetch):
            time.sleep(DELAY_SECONDS)

    # Final save
    save_cache(cache)

    print()
    print(f"Done! Fetched {fetched}, failed {failed}")
    print(f"Total cached: {len(cache)} sites")
    print(f"Cache saved: {EXTENDED_CACHE}")

    # Quick stats
    with_text = sum(1 for v in cache.values() if len(v.get("full_text", "")) > 100)
    print(f"Sites with substantial text: {with_text}")


if __name__ == "__main__":
    main()
