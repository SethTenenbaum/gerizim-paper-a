# data package — common UNESCO corpus library
# Re-export everything from the canonical module in data/scripts/
"""
Shim so that analysis scripts can do:
    from data.unesco_corpus import load_corpus, cultural_sites_with_coords
"""
from data.scripts.unesco_corpus import *          # noqa: F401,F403
from data.scripts.unesco_corpus import (          # explicit for IDE support
    load_corpus,
    cultural_sites_with_coords,
    search_sites,
    strip_html,
    beru_deviation,
    UNESCOSite,
    XML_PATH,
    EXTENDED_CACHE,
)
