"""
meta_keyword_test.py
=====================
Meta-test: of all possible single-keyword searches on the UNESCO corpus,
what fraction achieves enrichment at least as strong as our curated
circular-monument keyword list?

This directly addresses the peer-review concern:
  "Did you tune your keyword list to find the result?"

METHOD
======
1. Extract every unique English word (4+ letters) that appears in
   UNESCO short_descriptions (Cultural/Mixed sites with coordinates).
2. For each word, build the implied dataset (all UNESCO sites whose
   description contains that word).
3. Compute the same Tier-A binomial p-value and Tier-A hit rate.
4. Report:
   - Distribution of hit rates across all possible keywords
   - How many keywords achieve p ≤ our observed p-value
   - Exact percentile rank of our result
   - Expected number of "significant" keyword searches by chance

This gives a rigorous answer to: "Is the circular-monument result
special, or would any architectural keyword find this?"

ARCHITECTURE KEYWORD BASELINE
==============================
Our curated set (N=17) achieves:
  Tier-A:  5/17 = 29.4%, p = 0.2418  (underpowered)
  Tier-A+: 2/17 = 11.8%, p = 0.1465  (underpowered)
  χ²:      p = 0.027 (significant)

The meta-test asks: what fraction of all keywords achieves χ² ≤ 0.027?
"""

import re
import sys
import numpy as np
from pathlib import Path
from scipy.stats import binomtest, chisquare
from collections import Counter, defaultdict
import json

# ── Use common corpus library ──────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_A_MAX, TIER_B_MAX, P_NULL_AP, P_NULL_A,
    deviation as beru_deviation, tier_label, is_aplus,
    load_keywords,
)

# ── Config ────────────────────────────────────────────────────────────────────

# Our observed values from spherical_monument_test.py (N=17)
OBS_N     = 17
OBS_N_A   = 5
OBS_N_AP  = 2
OBS_P_A   = 0.2418
OBS_P_AP  = 0.1465
OBS_CHI_P = 0.0272

# ── Helpers ───────────────────────────────────────────────────────────────────
def compute_stats(lons):
    """Given a list of longitudes, compute tier counts and p-values."""
    n = len(lons)
    if n == 0:
        return None
    devs = [beru_deviation(lon) for lon in lons]
    n_ap = sum(1 for d in devs if d <= TIER_APLUS)
    n_a  = sum(1 for d in devs if d <= TIER_A_MAX)

    bt_ap = binomtest(n_ap, n, P_NULL_AP, alternative="greater")
    bt_a  = binomtest(n_a,  n, P_NULL_A,  alternative="greater")

    # Chi-square uniform (5 bins)
    obs_bins = [0] * 5
    for d in devs:
        obs_bins[min(int(d / 0.010), 4)] += 1
    try:
        chi_stat, chi_p = chisquare(obs_bins, f_exp=[n / 5.0] * 5)
    except Exception:
        chi_p = 1.0

    return {
        "n": n,
        "n_ap": n_ap,
        "n_a": n_a,
        "rate_ap": n_ap / n,
        "rate_a":  n_a  / n,
        "p_ap": bt_ap.pvalue,
        "p_a":  bt_a.pvalue,
        "chi_p": chi_p,
    }


# ── Load all Cultural/Mixed sites with coordinates ────────────────────────────
print("Loading UNESCO corpus...")
corpus = load_corpus()
_cultural = cultural_sites_with_coords(corpus)

all_sites = []
for s in _cultural:
    all_sites.append({
        "name":    s.site,
        "lon":     s.longitude,
        "desc_lc": s.short_description.lower(),
        "name_lc": s.site.lower(),
    })

print(f"Loaded {len(all_sites)} Cultural/Mixed sites with coordinates")

# ── Build inverted index: word → list of lons ────────────────────────────────
print("Building inverted index of words...")
word_to_lons = defaultdict(list)

for site in all_sites:
    # Use both description and site name
    text = site["desc_lc"] + " " + site["name_lc"]
    words = set(re.findall(r"[a-z]{4,}", text))
    for w in words:
        word_to_lons[w].append(site["lon"])

print(f"Unique words (4+ letters): {len(word_to_lons)}")

# ── Compute stats for our curated keyword set ─────────────────────────────────
print("\nComputing stats for CURATED keyword set...")
# ── Load curated keyword set from config ─────────────────────────────────────
FORM_KEYWORDS = load_keywords("dome_forms")
# Word-boundary regexes to avoid false matches (e.g., "dome" inside "domestic")
FORM_KW_RES = {kw: re.compile(r"\b" + re.escape(kw) + r"\b") for kw in FORM_KEYWORDS}

# Replicate the exact curated set match (word-boundary)
curated_lons = []
for site in all_sites:
    text = site["desc_lc"] + " " + site["name_lc"]
    if any(FORM_KW_RES[k].search(text) for k in FORM_KEYWORDS):
        curated_lons.append(site["lon"])

curated_stats = compute_stats(curated_lons)
print(f"  N={curated_stats['n']}, Tier-A={curated_stats['n_a']}/{curated_stats['n']} "
      f"({100*curated_stats['rate_a']:.1f}%), p_A={curated_stats['p_a']:.4f}, "
      f"χ²p={curated_stats['chi_p']:.4f}")

# ── Meta-test: scan all words with N in range [5, 200] ───────────────────────
print("\nRunning meta-test over all keyword candidates...")
print("(Testing every word that matches 5–200 UNESCO sites)")
print()

MIN_N = 5    # minimum sites to have power; below this too underpowered
MAX_N = 200  # above this, word is too common to be architectural

results = []
for word, lons in word_to_lons.items():
    n = len(lons)
    if n < MIN_N or n > MAX_N:
        continue
    st = compute_stats(lons)
    if st:
        st["word"] = word
        results.append(st)

print(f"Tested {len(results)} words with {MIN_N} ≤ N ≤ {MAX_N}")

# ── Analysis ──────────────────────────────────────────────────────────────────
results.sort(key=lambda x: x["chi_p"])

# Distribution of chi-square p-values
chi_ps  = np.array([r["chi_p"] for r in results])
p_a_vals = np.array([r["p_a"]  for r in results])
p_ap_vals = np.array([r["p_ap"] for r in results])
rate_a_vals = np.array([r["rate_a"] for r in results])

n_beat_chi   = np.sum(chi_ps  <= OBS_CHI_P)
n_beat_p_a   = np.sum(p_a_vals <= OBS_P_A)
n_beat_rate_a = np.sum(rate_a_vals >= OBS_N_A / OBS_N)

pct_chi   = float(n_beat_chi)  / len(results) * 100
pct_p_a   = float(n_beat_p_a)  / len(results) * 100
pct_rate_a = float(n_beat_rate_a) / len(results) * 100

print()
print("=" * 80)
print("  META-TEST RESULTS")
print("=" * 80)
print(f"""
  Our curated circular-monument keyword set result:
    N = {curated_stats['n']}, Tier-A = {curated_stats['n_a']}/{curated_stats['n']} = {100*curated_stats['rate_a']:.1f}%
    Binomial p (Tier-A) = {curated_stats['p_a']:.4f}
    χ²-uniform p        = {curated_stats['chi_p']:.4f}

  Meta-test population:
    {len(results)} keyword candidates (matched {MIN_N}–{MAX_N} sites each)

  How many random keywords achieve each threshold?
  ─────────────────────────────────────────────────────────────────────
  Metric                          Threshold      N keywords    Fraction
  ─────────────────────────────────────────────────────────────────────
  χ²-uniform p ≤ {OBS_CHI_P:.4f}           {OBS_CHI_P:.4f}        {n_beat_chi:>8d}    {pct_chi:>6.1f}%
  Tier-A binomial p ≤ {OBS_P_A:.4f}       {OBS_P_A:.4f}        {n_beat_p_a:>8d}    {pct_p_a:>6.1f}%
  Tier-A rate ≥ {100*OBS_N_A/OBS_N:.1f}%               {OBS_N_A/OBS_N:.4f}        {n_beat_rate_a:>8d}    {pct_rate_a:>6.1f}%
""")

# Expected by chance
exp_chi_sig  = len(results) * 0.05
exp_p_a_sig  = len(results) * 0.05
print(f"  Under pure chance (uniform null):")
print(f"    Expected keywords with χ²p ≤ 0.05:  {exp_chi_sig:.0f}  (we see {np.sum(chi_ps<=0.05)})")
print(f"    Expected keywords with p_A ≤ 0.05:  {exp_p_a_sig:.0f}  (we see {np.sum(p_a_vals<=0.05)})")
print()

# Top 30 words by χ² p-value (these are the "most enriched" architectural words)
print("  Top 30 keywords by χ²-uniform p (most strongly enriched):")
print(f"  {'Word':<25s}  {'N':>5}  {'Tier-A':>6}  {'Rate-A':>7}  {'p_A':>8}  {'χ²p':>8}")
print(f"  {'─'*70}")
for r in results[:30]:
    print(f"  {r['word']:<25s}  {r['n']:>5d}  {r['n_a']:>6d}  "
          f"{r['rate_a']:>7.3f}  {r['p_a']:>8.4f}  {r['chi_p']:>8.4f}")

# Where does "stupa" alone rank?
print()
stupa_results = [r for r in results if r["word"] == "stupa"]
if stupa_results:
    stupa_r = stupa_results[0]
    stupa_rank = results.index(stupa_r) + 1
    print(f"  'stupa' alone: N={stupa_r['n']}, Tier-A={stupa_r['n_a']}/{stupa_r['n']} "
          f"({100*stupa_r['rate_a']:.1f}%), p_A={stupa_r['p_a']:.4f}, χ²p={stupa_r['chi_p']:.4f}")
    print(f"  Rank by χ²p: #{stupa_rank} of {len(results)}")

# Where does "circular" alone rank?
for kw in ["circular", "tholos", "domed", "dome", "sphere", "spherical", "round", "stupa"]:
    kw_results = [r for r in results if r["word"] == kw]
    if kw_results:
        r = kw_results[0]
        rank = results.index(r) + 1
        print(f"  '{kw}': N={r['n']}, Tier-A={r['n_a']}/{r['n']} ({100*r['rate_a']:.1f}%), "
              f"p_A={r['p_a']:.4f}, χ²p={r['chi_p']:.4f}  [rank #{rank}]")

# ── Minimal keyword test: just "circular", "sphere", "dome" ──────────────────
print()
print("=" * 80)
print("  MINIMAL KEYWORD TEST: {'stupa/stupas', 'tholos', 'dome/domed/domes', 'spherical'}")
print("  (Minimal morphological set: round/hemispherical monumental architecture)")
print("  (Word-boundary matching to avoid false positives from 'domestic' etc.)")
print("=" * 80)

MINIMAL_KWS = FORM_KEYWORDS  # loaded from config; same set as the curated test above
MINIMAL_KW_RES = {kw: re.compile(r"\b" + re.escape(kw) + r"\b") for kw in MINIMAL_KWS}
minimal_lons = []
minimal_sites = []
for site in all_sites:
    text = site["desc_lc"] + " " + site["name_lc"]
    matched = [k for k in MINIMAL_KWS if MINIMAL_KW_RES[k].search(text)]
    if matched:
        minimal_lons.append(site["lon"])
        minimal_sites.append((site["name"], site["lon"], matched))

minimal_stats = compute_stats(minimal_lons)
if minimal_stats:
    print(f"\n  N={minimal_stats['n']}, "
          f"Tier-A+={minimal_stats['n_ap']}/{minimal_stats['n']} ({100*minimal_stats['rate_ap']:.1f}%), "
          f"Tier-A={minimal_stats['n_a']}/{minimal_stats['n']} ({100*minimal_stats['rate_a']:.1f}%)")
    print(f"  Binomial p (Tier-A+): {minimal_stats['p_ap']:.4f}")
    print(f"  Binomial p (Tier-A):  {minimal_stats['p_a']:.4f}")
    print(f"  χ²-uniform p:         {minimal_stats['chi_p']:.4f}")
    print()
    
    # Show per-site table
    devs = [(name, lon, beru_deviation(lon)) for name, lon, _ in minimal_sites]
    devs.sort(key=lambda x: x[2])
    print(f"  {'Site':<55s}  {'Lon':>8}  {'Dev':>8}  {'km':>6}  T")
    print(f"  {'─'*90}")
    for name, lon, dev in devs:
        km = dev * BERU * 111.0
        tier = tier_label(dev)
        mark = " ◀◀ A+" if tier in ("A++", "A+") else (" ◀ A" if tier == "A" else "")
        print(f"  {name:<55s}  {lon:>8.4f}  {dev:>8.5f}  {km:>6.1f}  {tier}{mark}")

# ── Save results for further analysis ────────────────────────────────────────
out = {
    "meta_test_n_keywords_tested": len(results),
    "meta_test_min_n": MIN_N,
    "meta_test_max_n": MAX_N,
    "curated_chi_p": curated_stats["chi_p"],
    "curated_p_a": curated_stats["p_a"],
    "curated_n": curated_stats["n"],
    "curated_n_a": curated_stats["n_a"],
    "n_keywords_beating_chi": int(n_beat_chi),
    "n_keywords_beating_p_a": int(n_beat_p_a),
    "fraction_beating_chi": pct_chi,
    "fraction_beating_p_a": pct_p_a,
    "minimal_n": minimal_stats["n"] if minimal_stats else None,
    "minimal_n_a": minimal_stats["n_a"] if minimal_stats else None,
    "minimal_p_a": minimal_stats["p_a"] if minimal_stats else None,
    "minimal_chi_p": minimal_stats["chi_p"] if minimal_stats else None,
    "top30": results[:30],
}
out_path = Path(__file__).parent / "meta_keyword_results.json"
with open(out_path, "w") as f:
    json.dump(out, f, indent=2)
print(f"\n  Results saved to {out_path}")

print()
print("=" * 80)
print("  INTERPRETATION")
print("=" * 80)
print(f"""
  If {pct_chi:.1f}% of all possible keywords achieve χ²p ≤ {OBS_CHI_P:.4f},
  then the significance of our curated result (χ²p={curated_stats['chi_p']:.4f}) 
  corrected for this search over the keyword space has an effective p-value of
  approximately {pct_chi/100:.4f} (i.e., {pct_chi:.1f}% of random keywords do as well).

  This is the critical number for the "multiple comparisons over keyword space"
  objection. If this fraction is < 5%, the curated result is in the top 5% of
  all possible keyword searches — which is meaningful even after this correction.
  If it's > 50%, the result is unremarkable.
""")
