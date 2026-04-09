"""
unesco_corpus.py — Single source of truth for the UNESCO World Heritage corpus.
================================================================================

Every script in this project should import from this module instead of
re-parsing the XML independently.

CAPABILITIES
────────────
  1. Parse `data/store/unesco/unesco.xml` → list of UNESCOSite dataclass objects
  2. Load cached **extended descriptions** scraped from UNESCO web pages
     (Brief Synthesis, Criteria, Integrity, Authenticity, full OUV text).
  3. Keyword search across BOTH the XML short_description AND the extended
     web-page text.
  4. Common beru-deviation math used by all analysis scripts.

EXTENDED DATA
─────────────
  The XML `short_description` is often abbreviated. The UNESCO web page for
  each site (e.g. https://whc.unesco.org/en/list/121) contains a much richer
  "Outstanding Universal Value" section, which includes keywords like "stupa"
  that do not appear in the XML. Run `data/fetch_extended.py` to download
  and cache this data into `data/store/unesco/extended_cache.json`.

USAGE
─────
  from data.unesco_corpus import load_corpus, search_sites, beru_deviation

  sites = load_corpus()                     # all 1,248 sites
  hits  = search_sites(sites, ["stupa"])    # keyword search (XML + extended)
"""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

# ── Paths ─────────────────────────────────────────────────────────────────────
_DATA_DIR       = Path(__file__).parent.parent / "store" / "unesco"
XML_PATH        = _DATA_DIR / "unesco.xml"
EXTENDED_CACHE  = _DATA_DIR / "extended_cache.json"

# ── Constants ─────────────────────────────────────────────────────────────────
# These are now the canonical source from lib/beru.py and config.json.
# Backward-compatible aliases are provided so existing analysis scripts
# that import from data.unesco_corpus continue to work during migration.
# NEW CODE should import from lib.beru instead.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).parent.parent))

from lib.beru import (
    GERIZIM, BERU, TIER_APP, TIER_APLUS, TIER_A_MAX, TIER_B_MAX,
    P_NULL_AP, P_NULL_A,
    deviation as _deviation_impl,
    tier_label, is_aplus, is_a_or_better, dev_to_km,
)


# Legacy aliases
def beru_deviation(lon: float, anchor: float = GERIZIM, beru: float = BERU) -> float:
    """Compute the deviation of a longitude from the nearest 0.1-beru harmonic."""
    return _deviation_impl(lon, anchor=anchor, beru=beru)


def tier_km(dev: float) -> float:
    """Convert beru deviation to approximate km (at equator)."""
    return dev_to_km(dev)


def sig_label(p: float) -> str:
    """Return a significance label for a p-value."""
    from lib.stats import significance_label
    return significance_label(p)


# ── Dataclass ─────────────────────────────────────────────────────────────────
@dataclass
class UNESCOSite:
    """A single UNESCO World Heritage site parsed from the XML + extended cache."""

    # Core XML fields
    id_number:       str = ""
    unique_number:   str = ""
    site:            str = ""
    http_url:        str = ""
    image_url:       str = ""
    iso_code:        str = ""
    category:        str = ""           # Cultural / Natural / Mixed
    criteria_txt:    str = ""
    date_inscribed:  str = ""
    states:          str = ""
    regions:         str = ""
    location:        str = ""
    transnational:   str = "0"
    danger:          str = ""
    extension:       str = "0"
    secondary_dates: str = ""
    dossier:         str = ""

    # Text fields (HTML-stripped)
    short_description: str = ""         # from XML <short_description>
    justification:     str = ""         # from XML <justification>

    # Geocoordinates (first POI)
    latitude:  Optional[float] = None
    longitude: Optional[float] = None

    # All POIs
    all_pois: List[Dict[str, float]] = field(default_factory=list)

    # Extended web-page text (loaded from cache)
    extended_description: str = ""      # Full OUV text from UNESCO page
    brief_synthesis:      str = ""      # Just the Brief Synthesis section
    criteria_detail:      str = ""      # Criterion-by-criterion detail

    # ── Computed helpers ──────────────────────────────────────────────────
    @property
    def has_coords(self) -> bool:
        return self.latitude is not None and self.longitude is not None

    @property
    def is_cultural_or_mixed(self) -> bool:
        return self.category in ("Cultural", "Mixed")

    @property
    def year(self) -> Optional[int]:
        try:
            return int(self.date_inscribed.strip())
        except (ValueError, AttributeError):
            return None

    @property
    def full_text(self) -> str:
        """All available text, lowercased, for keyword searching."""
        parts = [
            self.site,
            self.short_description,
            self.justification,
            self.extended_description,
        ]
        return "\n".join(p for p in parts if p).lower()

    @property
    def xml_text(self) -> str:
        """Only XML-sourced text (short_description + justification + name), lowercased."""
        parts = [self.site, self.short_description, self.justification]
        return "\n".join(p for p in parts if p).lower()


# ── HTML stripping ────────────────────────────────────────────────────────────
_HTML_RE = re.compile(r"<[^>]+>")

def strip_html(text: str) -> str:
    """Remove HTML tags from a string."""
    if not text:
        return ""
    return _HTML_RE.sub(" ", text).strip()


# ── XML parsing ───────────────────────────────────────────────────────────────
def _parse_xml(xml_path: Path = XML_PATH) -> List[UNESCOSite]:
    """Parse the UNESCO XML into a list of UNESCOSite objects."""
    tree = ET.parse(str(xml_path))
    root = tree.getroot()
    sites: List[UNESCOSite] = []

    for row in root.findall("row"):
        def get(tag: str) -> str:
            el = row.find(tag)
            return (el.text or "").strip() if el is not None and el.text else ""

        site = UNESCOSite(
            id_number       = get("id_number"),
            unique_number   = get("unique_number"),
            site            = get("site"),
            http_url        = get("http_url"),
            image_url       = get("image_url"),
            iso_code        = get("iso_code"),
            category        = get("category"),
            criteria_txt    = get("criteria_txt"),
            date_inscribed  = get("date_inscribed"),
            states          = get("states"),
            regions         = get("regions"),
            location        = get("location"),
            transnational   = get("transnational"),
            danger          = get("danger"),
            extension       = get("extension"),
            secondary_dates = get("secondary_dates"),
            dossier         = get("dossier"),
            short_description = strip_html(get("short_description")),
            justification     = strip_html(get("justification")),
        )

        # Parse geolocations — first POI becomes lat/lon
        geo = row.find("geolocations")
        if geo is not None:
            for poi_el in geo.findall("poi"):
                try:
                    lat = float(poi_el.find("latitude").text)
                    lon = float(poi_el.find("longitude").text)
                    iso2 = (poi_el.find("iso2").text or "").strip() if poi_el.find("iso2") is not None else ""
                    site.all_pois.append({"lat": lat, "lon": lon, "iso2": iso2})
                except (TypeError, ValueError, AttributeError):
                    continue

            if site.all_pois:
                site.latitude = site.all_pois[0]["lat"]
                site.longitude = site.all_pois[0]["lon"]

        sites.append(site)

    return sites


# ── Extended cache loading ────────────────────────────────────────────────────
def _load_extended_cache(cache_path: Path = EXTENDED_CACHE) -> Dict[str, dict]:
    """Load the extended description cache (id_number → extended data)."""
    if not cache_path.exists():
        return {}
    with open(cache_path, "r", encoding="utf-8") as f:
        return json.load(f)


def _merge_extended(sites: List[UNESCOSite], cache: Dict[str, dict]) -> None:
    """Merge extended web-page data into site objects (in-place)."""
    for site in sites:
        ext = cache.get(site.id_number)
        if ext is None:
            continue
        site.extended_description = ext.get("full_text", "")
        site.brief_synthesis      = ext.get("brief_synthesis", "")
        site.criteria_detail      = ext.get("criteria_detail", "")


# ── Main entry point ──────────────────────────────────────────────────────────
_CORPUS_CACHE: Optional[List[UNESCOSite]] = None


def load_corpus(
    xml_path:   Path = XML_PATH,
    cache_path: Path = EXTENDED_CACHE,
    force_reload: bool = False,
) -> List[UNESCOSite]:
    """
    Load the full UNESCO corpus from XML, enriched with extended web data.

    Results are cached in-process; call with force_reload=True to re-parse.

    Returns:
        List of UNESCOSite objects (typically ~1,248 entries).
    """
    global _CORPUS_CACHE
    if _CORPUS_CACHE is not None and not force_reload:
        return _CORPUS_CACHE

    sites = _parse_xml(xml_path)
    cache = _load_extended_cache(cache_path)
    if cache:
        _merge_extended(sites, cache)
        _n_ext = sum(1 for s in sites if s.extended_description)
        print(f"[unesco_corpus] Loaded {len(sites)} sites, "
              f"{_n_ext} with extended descriptions")
    else:
        print(f"[unesco_corpus] Loaded {len(sites)} sites "
              f"(no extended cache — run data/fetch_extended.py to download)")

    _CORPUS_CACHE = sites
    return sites


# ── Search / filter helpers ───────────────────────────────────────────────────
def search_sites(
    sites: List[UNESCOSite],
    keywords: Sequence[str],
    *,
    use_extended: bool = True,
    word_boundary: bool = True,
    cultural_only: bool = False,
    with_coords: bool = False,
) -> List[UNESCOSite]:
    """
    Return sites whose text matches ANY of the given keywords.

    Args:
        sites:         List of UNESCOSite objects.
        keywords:      List of keywords to search for (case-insensitive).
        use_extended:  If True, search full_text (XML + extended web data).
                       If False, search only xml_text.
        word_boundary: If True, use \\b word boundaries around each keyword.
        cultural_only: If True, exclude Natural-category sites.
        with_coords:   If True, exclude sites without geocoordinates.

    Returns:
        Filtered list of matching UNESCOSite objects.
    """
    if word_boundary:
        patterns = [re.compile(r"\b" + re.escape(kw.lower()) + r"\b") for kw in keywords]
    else:
        patterns = [re.compile(re.escape(kw.lower())) for kw in keywords]

    results = []
    for site in sites:
        if cultural_only and not site.is_cultural_or_mixed:
            continue
        if with_coords and not site.has_coords:
            continue

        text = site.full_text if use_extended else site.xml_text
        if any(p.search(text) for p in patterns):
            results.append(site)

    return results


def cultural_sites_with_coords(sites: List[UNESCOSite]) -> List[UNESCOSite]:
    """Return only Cultural/Mixed sites that have geocoordinates."""
    return [s for s in sites if s.is_cultural_or_mixed and s.has_coords]


# ── CLI self-test ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import sys

    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)

    print(f"\nTotal sites:           {len(corpus)}")
    print(f"Cultural/Mixed:        {sum(1 for s in corpus if s.is_cultural_or_mixed)}")
    print(f"With coordinates:      {sum(1 for s in corpus if s.has_coords)}")
    print(f"Cultural + coords:     {len(cultural)}")
    print(f"With extended text:    {sum(1 for s in corpus if s.extended_description)}")

    # Example keyword search
    test_kw = sys.argv[1:] if len(sys.argv) > 1 else ["stupa"]
    print(f"\nKeyword search: {test_kw}")

    # XML-only search
    xml_hits = search_sites(corpus, test_kw, use_extended=False,
                            cultural_only=True, with_coords=True)
    print(f"  XML-only matches:    {len(xml_hits)}")
    for s in xml_hits[:5]:
        print(f"    - {s.site} ({s.states})")

    # Full search (XML + extended)
    full_hits = search_sites(corpus, test_kw, use_extended=True,
                             cultural_only=True, with_coords=True)
    print(f"  Full-text matches:   {len(full_hits)}")
    for s in full_hits[:10]:
        print(f"    - {s.site} ({s.states})")

    # Show sites found only via extended
    xml_ids = {s.id_number for s in xml_hits}
    ext_only = [s for s in full_hits if s.id_number not in xml_ids]
    if ext_only:
        print(f"\n  Found ONLY via extended text ({len(ext_only)}):")
        for s in ext_only:
            print(f"    + {s.site} ({s.states})")
