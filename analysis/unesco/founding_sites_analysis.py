"""
founding_sites_analysis.py
==========================
Thematic breakdown of the Tier-A+ sites detected in the corpus.

The 57 sites within ±6.7 km of a beru harmonic are classified by
UNESCO's own inscription language into five categories.

IMPORTANT — NO HAND-PICKING
─────────────────────────────
Every site in this analysis is included because it passes the same
mechanical beru-deviation threshold (≤ 0.002 beru) applied uniformly
to all 1011 cultural/mixed UNESCO sites. No site is added or removed
by hand. The old APLUS_CATALOGUE list (which assigned categories by
hand to pre-selected sites) has been removed entirely.

Classification uses keyword lists from config.json, applied to the
UNESCO short_description text. The same keyword lists are used in
founding_enrichment_test.py so results are directly comparable.

Category codes (from config.json):
  F = founding_capital     — first/original seat of a kingdom/empire/state
  S = sacred_origin        — birthplace/founding sanctuary of a religion
  M = founding_monument    — prototype monument that defined a tradition
  X = founding_axis        — primary imperial road, wall, or trade corridor
  L = ancient_landscape    — oldest continuous cultural landscape
  ? = unclassified by any keyword (shown last, for transparency)

STATISTICAL STRUCTURE
─────────────────────
  H₀: A+ sites are a random draw from all UNESCO Cultural sites
  H₁: A+ sites are disproportionately "founding" sites

  Operationalization: UNESCO's OWN inscription year as a proxy for
  importance. Sites inscribed in 1978–1984 (the committee's strictest
  phase) represent the international consensus "canon" of heritage.

  Prediction: If the grid marks founding importance, the earliest
  inscribed sites should have the highest A+ enrichment, and this
  enrichment should dilute as the list expands to include lesser sites.

ANCHOR: Mount Gerizim 35.274°E  |  BERU = 30.0°  |  1 beru = 80 nodes
"""

import sys
from pathlib import Path
from scipy.stats import binomtest, fisher_exact

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_A_MAX, TIER_B_MAX, P_NULL_AP,
    deviation as beru_dev, tier_label as tier, is_aplus, is_a_or_better,
)
from lib.founding_filter import (
    classify_site as _classify_site_obj,
    CATEGORY_LABELS as CAT_LABELS,
    PRIORITY,
    primary_category,
)
from lib.stats import significance_label as sig
from lib.reporting import (
    print_header, print_thematic_breakdown, print_unclassified_sites,
    print_top_sites, print_enrichment_table, print_fisher_inline,
)

# ── No keyword lists loaded directly — context-aware filter handles all ───────


def classify_site_text(site_obj) -> set:
    """
    Return a set of category codes using the context-aware founding filter.

    Uses lib.founding_filter.classify_site() which performs:
      1. Unambiguous keyword matching (always accepted)
      2. Ambiguous keyword matching with sentence-level context validation
    """
    cats = _classify_site_obj(site_obj)
    return set(cats.keys()) if cats else set()


def primary_cat(s):
    for c in PRIORITY:
        if c in s["cats"]:
            return c
    return "?"


# ── Load corpus ───────────────────────────────────────────────────────────────
corpus = load_corpus()
_cultural = cultural_sites_with_coords(corpus)

sites = []
for s in _cultural:
    yr = s.year
    if yr is None:
        continue
    d = beru_dev(s.longitude)
    cats = classify_site_text(s)
    sites.append({
        "name": s.site, "lon": s.longitude, "yr": yr,
        "dev": d, "tier": tier(d), "km": d * BERU * 111.0,
        "cats": cats, "text": s.full_text,
    })

# ── Derive the A+ set mechanically from the corpus ───────────────────────────
aplus_sites = sorted(
    [s for s in sites if is_aplus(s["tier"])],
    key=lambda s: s["dev"],
)

# ── 1. Thematic breakdown table ───────────────────────────────────────────────
print_header(
    f"BERU-GRID A+ SITES: THEMATIC ANALYSIS\n"
    f"  Tier-A+ = within ±6.7 km of a beru harmonic  "
    f"(1/80 of Earth's circumference)\n"
    f"  Total A+ sites detected from corpus: {len(aplus_sites)}\n"
    f"  Classification: keyword-driven from config.json, applied uniformly",
    width=95,
)

by_cat = {}
for s in aplus_sites:
    by_cat.setdefault(primary_cat(s), []).append(s)

# KW_LISTS for backward compat with print_thematic_breakdown (shows matched kw)
from lib.beru import load_keywords
KW_LISTS = {
    cat: load_keywords(key) if key else []
    for cat, key in [("F", "founding_capital"), ("S", "sacred_origin"),
                     ("M", "founding_monument"), ("X", "founding_axis"),
                     ("L", "ancient_landscape"), ("?", None)]
}
print_thematic_breakdown(by_cat, CAT_LABELS, KW_LISTS, PRIORITY)
print_unclassified_sites(by_cat.get("?", []))

# ── 2. Enrichment by inscription cohort ──────────────────────────────────────
def cohort_stats(subset):
    n   = len(subset)
    nap = sum(1 for s in subset if is_aplus(s["tier"]))
    p   = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue if n else 1.0
    return n, nap, p

early = [s for s in sites if s["yr"] <= 1984]
pre2k = [s for s in sites if s["yr"] <  2000]
late  = [s for s in sites if s["yr"] >= 2000]

ne,  nap_e,  pe  = cohort_stats(early)
np2, nap_p2, pp2 = cohort_stats(pre2k)
nl,  nap_l,  pl  = cohort_stats(late)
na2, nap_a,  pa  = cohort_stats(sites)

or_fl, p_fl = fisher_exact(
    [[nap_e, ne - nap_e], [nap_l, nl - nap_l]],
    alternative="greater",
)

print_header(
    "QUANTITATIVE SUMMARY: ENRICHMENT BY INSCRIPTION COHORT\n"
    "  (UNESCO inscription year = proxy for importance/strictness)",
    width=95,
)
print_enrichment_table([
    {"label": "1978–1984  (first 7 yr)",  "n": ne,  "n_aplus": nap_e,  "p": pe},
    {"label": "1978–1999  (pre-2000)",    "n": np2, "n_aplus": nap_p2, "p": pp2},
    {"label": "2000–2025  (modern era)",  "n": nl,  "n_aplus": nap_l,  "p": pl},
    {"label": "ALL (1978–2025)",          "n": na2, "n_aplus": nap_a,  "p": pa},
], null_rate=P_NULL_AP, label_width=24)

print()
print_fisher_inline("Fisher exact (1978–1984 vs 2000–2025", or_fl, p_fl)

# ── 3. Thematic count summary bar chart ──────────────────────────────────────
print_header(
    f"THEMATIC BREAKDOWN OF ALL {len(aplus_sites)} A+ SITES (keyword-classified)",
    width=95,
)
counts = {c: len(by_cat.get(c, [])) for c in PRIORITY}
total_ap = len(aplus_sites)
for cat in PRIORITY:
    n = counts[cat]
    if n == 0:
        continue
    pct = 100 * n / total_ap
    bar = "█" * n
    print(f"  {CAT_LABELS[cat]:<45}  {n:>3}  ({pct:4.1f}%)  {bar}")

founding_n = sum(counts.get(c, 0) for c in "FSMX")
print(f"""
  COMBINED 'FOUNDING' (F+S+M+X):  {founding_n} / {total_ap}  ({100*founding_n/total_ap:.0f}%)
  UNCLASSIFIED (?):                {counts.get('?', 0)} / {total_ap}  ({100*counts.get('?',0)/total_ap:.0f}%)

  The classification is mechanical: a site receives a category if and only
  if its UNESCO description contains one of the config.json keyword strings
  for that category. No site is assigned by hand.

  Sites in the '?' category match no keyword — these are not excluded from
  the enrichment tests (they remain in the full corpus), but they do not
  contribute to the founding-site count. The keyword lists in config.json
  can be extended to classify them if appropriate.
""")

# ── 4. Top 10 most precise hits ───────────────────────────────────────────────
print_top_sites(aplus_sites, n=10)


