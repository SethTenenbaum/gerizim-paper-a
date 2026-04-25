"""
meta_keyword_audit.py
=====================
Reviewer-facing full listing of all candidate words from the meta-keyword
test (§5.4 of the paper).

WHAT THIS PRODUCES
------------------
  1. All N candidate words (4+ letters, matched 5–200 UNESCO sites) with stats:
       word | N | A+ | A+% | A | A% | χ²p
  2. The 57 words (3.3%) that beat the curated dome keyword χ² threshold —
     marked with *** in the listing.
  3. Where each dome keyword (stupa, tholos, dome, …) ranks among all words.
  4. Sorted output: alphabetical AND by χ² p-value ascending.
  5. Plain-text audit file written to analysis/unesco/meta_keyword_audit.txt.

USAGE
-----
    cd /path/to/gerizim-paper-a
    python3 analysis/unesco/meta_keyword_audit.py
"""

import re
import sys
import json
from pathlib import Path
from collections import defaultdict
from scipy.stats import binomtest, chisquare

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_A_MAX, P_NULL_AP, P_NULL_A,
    deviation as beru_deviation,
    load_keywords,
)

# ── Config ─────────────────────────────────────────────────────────────────────
MIN_N = 5
MAX_N = 200
SEP  = "═" * 80
SEP2 = "─" * 80


# ── Load corpus ────────────────────────────────────────────────────────────────
print("Loading UNESCO corpus…")
corpus   = load_corpus()
_cultural = cultural_sites_with_coords(corpus)

all_sites = []
for s in _cultural:
    all_sites.append({
        "name":    s.site,
        "lon":     s.longitude,
        "desc_lc": s.short_description.lower(),
        "name_lc": s.site.lower(),
    })
print(f"  {len(all_sites)} Cultural/Mixed sites with coordinates")


# ── Build inverted index ───────────────────────────────────────────────────────
print("Building word index…")
word_to_lons = defaultdict(list)
for site in all_sites:
    text  = site["desc_lc"] + " " + site["name_lc"]
    words = set(re.findall(r"[a-z]{4,}", text))
    for w in words:
        word_to_lons[w].append(site["lon"])

print(f"  {len(word_to_lons)} unique 4+ letter words found")


# ── Compute stats per word ─────────────────────────────────────────────────────
def compute_stats(lons):
    n = len(lons)
    if n == 0:
        return None
    devs = [beru_deviation(lon) for lon in lons]
    n_ap = sum(1 for d in devs if d <= TIER_APLUS)
    n_a  = sum(1 for d in devs if d <= TIER_A_MAX)
    p_ap = binomtest(n_ap, n, P_NULL_AP, alternative="greater").pvalue
    p_a  = binomtest(n_a,  n, P_NULL_A,  alternative="greater").pvalue
    obs_bins = [0] * 5
    for d in devs:
        obs_bins[min(int(d / TIER_A_MAX), 4)] += 1
    try:
        _, chi_p = chisquare(obs_bins, f_exp=[n / 5.0] * 5)
    except Exception:
        chi_p = 1.0
    return {"n": n, "n_ap": n_ap, "n_a": n_a, "p_ap": p_ap, "p_a": p_a, "chi_p": chi_p}


print("Computing stats for all candidate words…")
results = []
for word, lons in word_to_lons.items():
    if MIN_N <= len(lons) <= MAX_N:
        st = compute_stats(lons)
        if st:
            st["word"] = word
            results.append(st)

print(f"  {len(results)} words with {MIN_N} ≤ N ≤ {MAX_N}")


# ── Curated dome keyword set stats ────────────────────────────────────────────
FORM_KEYWORDS = load_keywords("dome_forms")
FORM_KW_RES   = {kw: re.compile(r"\b" + re.escape(kw) + r"\b") for kw in FORM_KEYWORDS}

curated_lons = []
for site in all_sites:
    text = site["desc_lc"] + " " + site["name_lc"]
    if any(FORM_KW_RES[k].search(text) for k in FORM_KEYWORDS):
        curated_lons.append(site["lon"])

curated_stats = compute_stats(curated_lons)
curated_chi_p = curated_stats["chi_p"]

# Mark words that beat the curated threshold
results_by_chi = sorted(results, key=lambda r: r["chi_p"])
n_beat         = sum(1 for r in results if r["chi_p"] <= curated_chi_p)
pct_beat       = 100.0 * n_beat / len(results)

# Rank lookup (1-based, sorted by chi_p ascending = best first)
chi_rank = {r["word"]: i + 1 for i, r in enumerate(results_by_chi)}


# ── Build report ───────────────────────────────────────────────────────────────
lines = []
def w(s=""): lines.append(s)

w(SEP)
w("  META-KEYWORD CANDIDATE WORDS AUDIT")
w("  meta_keyword_audit.py  —  auto-generated")
w(SEP)
w()
w(f"  UNESCO corpus:              {len(all_sites)} Cultural/Mixed sites with coordinates")
w(f"  Unique 4+ letter words:     {len(word_to_lons)}")
w(f"  Candidate words (5≤N≤200):  {len(results)}")
w()
w(f"  Curated dome keyword set:   N={curated_stats['n']}  χ²p={curated_chi_p:.4f}")
w(f"  Words beating curated χ²p:  {n_beat}  ({pct_beat:.1f}%)")
w()
w(f"  Dome keywords: {', '.join(sorted(FORM_KEYWORDS))}")
w()

# ── Section 1: Dome keyword rankings ──────────────────────────────────────────
w(SEP2)
w("  1. DOME KEYWORD INDIVIDUAL RANKINGS")
w(SEP2)
w(f"  {'keyword':<20}  {'N':>5}  {'A+':>4}  {'A+%':>6}  {'A':>4}  {'A%':>6}  {'χ²p':>8}  {'rank':>6}")
w(f"  {'─'*20}  {'─'*5}  {'─'*4}  {'─'*6}  {'─'*4}  {'─'*6}  {'─'*8}  {'─'*6}")

for kw in sorted(FORM_KEYWORDS):
    row = next((r for r in results if r["word"] == kw), None)
    if row:
        rank = chi_rank[kw]
        w(f"  {kw:<20}  {row['n']:>5}  {row['n_ap']:>4}  {100*row['n_ap']/row['n']:>5.1f}%  "
          f"{row['n_a']:>4}  {100*row['n_a']/row['n']:>5.1f}%  {row['chi_p']:>8.4f}  #{rank}")
    else:
        w(f"  {kw:<20}  — (not in candidate range)")

w()

# ── Section 2: Top 50 words by χ² p-value ─────────────────────────────────────
w(SEP2)
w("  2. TOP 50 WORDS BY χ² p-VALUE  (*** = beats curated dome result)")
w(SEP2)
w(f"  {'rank':>5}  {'word':<22}  {'N':>5}  {'A+':>4}  {'A+%':>6}  {'A':>4}  {'A%':>6}  {'χ²p':>8}")
w(f"  {'─'*5}  {'─'*22}  {'─'*5}  {'─'*4}  {'─'*6}  {'─'*4}  {'─'*6}  {'─'*8}")

for i, r in enumerate(results_by_chi[:50]):
    flag = " ***" if r["chi_p"] <= curated_chi_p else ""
    w(f"  {i+1:>5}  {r['word']:<22}  {r['n']:>5}  {r['n_ap']:>4}  "
      f"{100*r['n_ap']/r['n']:>5.1f}%  {r['n_a']:>4}  {100*r['n_a']/r['n']:>5.1f}%  "
      f"{r['chi_p']:>8.4f}{flag}")

w()

# ── Section 3: All words that beat the curated threshold ──────────────────────
w(SEP2)
w(f"  3. ALL {n_beat} WORDS BEATING CURATED χ²p ≤ {curated_chi_p:.4f}  (χ² ranked)")
w(SEP2)
w(f"  {'rank':>5}  {'word':<22}  {'N':>5}  {'A+':>4}  {'A+%':>6}  {'A':>4}  {'A%':>6}  {'χ²p':>8}  dome?")
w(f"  {'─'*5}  {'─'*22}  {'─'*5}  {'─'*4}  {'─'*6}  {'─'*4}  {'─'*6}  {'─'*8}  {'─'*5}")

for i, r in enumerate(results_by_chi):
    if r["chi_p"] > curated_chi_p:
        break
    is_dome = "YES" if r["word"] in FORM_KEYWORDS else ""
    w(f"  {i+1:>5}  {r['word']:<22}  {r['n']:>5}  {r['n_ap']:>4}  "
      f"{100*r['n_ap']/r['n']:>5.1f}%  {r['n_a']:>4}  {100*r['n_a']/r['n']:>5.1f}%  "
      f"{r['chi_p']:>8.4f}  {is_dome}")

w()

# ── Section 4: Full alphabetical listing ──────────────────────────────────────
w(SEP2)
w(f"  4. FULL ALPHABETICAL LISTING  ({len(results)} candidate words)")
w(f"     *** = beats curated χ²p,  dome = curated keyword")
w(SEP2)
w(f"  {'word':<22}  {'N':>5}  {'A+':>4}  {'A+%':>6}  {'A':>4}  {'A%':>6}  {'χ²p':>8}  {'χ²rank':>7}  flags")
w(f"  {'─'*22}  {'─'*5}  {'─'*4}  {'─'*6}  {'─'*4}  {'─'*6}  {'─'*8}  {'─'*7}  {'─'*10}")

for r in sorted(results, key=lambda x: x["word"]):
    rank  = chi_rank[r["word"]]
    flags = []
    if r["chi_p"] <= curated_chi_p:
        flags.append("***")
    if r["word"] in FORM_KEYWORDS:
        flags.append("dome")
    flag_str = " ".join(flags)
    w(f"  {r['word']:<22}  {r['n']:>5}  {r['n_ap']:>4}  "
      f"{100*r['n_ap']/r['n']:>5.1f}%  {r['n_a']:>4}  {100*r['n_a']/r['n']:>5.1f}%  "
      f"{r['chi_p']:>8.4f}  #{rank:<6}  {flag_str}")

w()
w(SEP)

report = "\n".join(lines)
print(report)

# ── Write audit file ───────────────────────────────────────────────────────────
out_path = Path(__file__).parent / "meta_keyword_audit.txt"
out_path.write_text(report)
print(f"\n  Audit written to: {out_path}")
