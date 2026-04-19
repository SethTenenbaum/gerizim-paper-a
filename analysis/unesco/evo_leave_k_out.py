"""
evo_leave_k_out.py
==================
Leave-2-out and leave-3-out sensitivity analysis for the hemispherical
mound-evolution corpus (Test 2b, raw keyword sweep, N=110, A+=13).

Mirrors dome_leave_k_out.py but applies to the combined mound+stupa+dome
population used in Test 2b.

Usage:
    python3 analysis/unesco/evo_leave_k_out.py
"""

from __future__ import annotations

import itertools
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from data.unesco_corpus import load_corpus
from lib.beru import BERU, TIER_APLUS, P_NULL_AP, deviation, tier_label, is_aplus
from scipy.stats import binomtest

# ── Load keyword sets ────────────────────────────────────────────────────────
_KW_PATH = Path(__file__).parent.parent.parent / "keywords.json"
with open(_KW_PATH) as _f:
    _KW = json.load(_f)

_evo = _KW["mound_evolution"]
MOUND_KEYWORDS = _evo["mound_unambiguous"] + _evo["mound_ambiguous"]
STUPA_KEYWORDS = _evo["stupa"]
DOME_KEYWORDS  = _evo["dome"]
ALL_KEYWORDS   = MOUND_KEYWORDS + STUPA_KEYWORDS + DOME_KEYWORDS

KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    for kw in ALL_KEYWORDS
}


def _load_evo_sites():
    corpus = load_corpus()
    sites = []
    for site in corpus:
        if site.category == "Natural":
            continue
        if not site.has_coords:
            continue
        txt = site.full_text
        if not any(r.search(txt) for r in KEYWORD_RES.values()):
            continue
        dev = deviation(site.longitude)
        dev_km = dev * BERU * 111.0
        ap = is_aplus(tier_label(dev))
        sites.append({
            "name":   site.site,
            "lon":    site.longitude,
            "dev":    dev,
            "dev_km": dev_km,
            "ap":     ap,
        })
    return sites


def _lko_pvalue(all_sites, removed_indices):
    remaining = [s for i, s in enumerate(all_sites) if i not in removed_indices]
    n   = len(remaining)
    nap = sum(1 for s in remaining if s["ap"])
    p   = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue
    rate = 100.0 * nap / n if n else 0.0
    return n, nap, rate, p


def main():
    sites = _load_evo_sites()
    n_total = len(sites)
    ap_sites = [(i, s) for i, s in enumerate(sites) if s["ap"]]
    n_ap = len(ap_sites)
    ap_sites_sorted = sorted(ap_sites, key=lambda x: x[1]["dev"])

    print("=" * 72)
    print(f"  evo_leave_k_out.py")
    print(f"  Evolution corpus (Test 2b raw sweep): N={n_total}, A+={n_ap}")
    print("=" * 72)

    print(f"\n  A+ sites ranked by closeness to harmonic:")
    for rank, (idx, s) in enumerate(ap_sites_sorted, 1):
        print(f"  {rank:>3}.  {s['dev_km']:>6.2f} km  {s['name']}")

    top3 = ap_sites_sorted[:3]
    idx1, s1 = top3[0]
    idx2, s2 = top3[1]
    idx3, s3 = top3[2]

    lko2_n, lko2_ap, lko2_rate, lko2_p = _lko_pvalue(sites, {idx1, idx2})
    lko3_n, lko3_ap, lko3_rate, lko3_p = _lko_pvalue(sites, {idx1, idx2, idx3})

    print(f"\n  Leave-2-out (drop closest 2): N={lko2_n}, A+={lko2_ap}, p={lko2_p:.4f}")
    print(f"  Leave-3-out (drop closest 3): N={lko3_n}, A+={lko3_ap}, p={lko3_p:.4f}")

    all_ap_idx = [i for i, _ in ap_sites]

    worst2_p = max(
        _lko_pvalue(sites, set(combo))[3]
        for combo in itertools.combinations(all_ap_idx, 2)
    )
    worst3_p = max(
        _lko_pvalue(sites, set(combo))[3]
        for combo in itertools.combinations(all_ap_idx, 3)
    )

    n2 = sum(1 for _ in itertools.combinations(all_ap_idx, 2))
    n3 = sum(1 for _ in itertools.combinations(all_ap_idx, 3))

    print(f"\n  Worst-case leave-2-out (all {n2} pairs):    p = {worst2_p:.4f}")
    print(f"  Worst-case leave-3-out (all {n3} triples):  p = {worst3_p:.4f}")
    print(f"\n  Bonferroni threshold (k=3): p < 0.0167")
    print(f"  Leave-3-out {'PASSES' if worst3_p < 0.0167 else 'FAILS'} Bonferroni at k=3")

    print(f"\n  % LaTeX macros (GROUP EVOLKO):")
    def fmt(p): return "$< 0.001$" if p < 0.001 else f"{p:.4f}"
    print(f"  \\newcommand{{\\EvoLKOworstTwoP}}{{{fmt(worst2_p)}}}  % evo LKO worst C({n_ap},2)={n2}")
    print(f"  \\newcommand{{\\EvoLKOworstThreeP}}{{{fmt(worst3_p)}}}  % evo LKO worst C({n_ap},3)={n3}")
    print(f"  \\newcommand{{\\EvoLKOnCombosTwo}}{{{n2}}}  % evo C({n_ap},2)")
    print(f"  \\newcommand{{\\EvoLKOnCombosThree}}{{{n3}}}  % evo C({n_ap},3)")


if __name__ == "__main__":
    main()
