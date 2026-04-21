"""
dome_leave_k_out.py
===================
Leave-2-out and leave-3-out sensitivity analysis for the dome/spherical
monument corpus (Test 2, raw keyword sweep, N=90, A+=11).

Targets the three closest-to-harmonic (lowest beru-deviation) A+ sites
in the dome corpus.  Removes them in combinations of 2 and 3, re-runs
the one-sided binomial test against the geometric null (P_NULL_AP from config), and reports
the worst-case (highest) p-value across all combinations.

MACROS EMITTED (GROUP LOO2)
---------------------------
  Site deviation macros (km):
    \\LKOsiteOneDevKm     — closest A+ site: km deviation
    \\LKOsiteTwoDevKm     — 2nd-closest A+ site: km deviation
    \\LKOsiteThreeDevKm   — 3rd-closest A+ site: km deviation
    \\LKOsiteOneName      — closest A+ site: short name (text)
    \\LKOsiteTwoName      — 2nd-closest A+ site: short name (text)
    \\LKOsiteThreeName    — 3rd-closest A+ site: short name (text)

  Leave-2-out (drop 2 closest):
    \\LKOtwoN             — N after removing 2 closest sites
    \\LKOtwoAp            — A+ after removing 2 closest sites
    \\LKOtwoRate          — A+ rate (%) after removing 2 closest
    \\LKOtwoP             — p-value after removing 2 closest (worst-case)

  Leave-3-out (drop all 3):
    \\LKOthreeN           — N after removing 3 closest sites
    \\LKOthreeAp          — A+ after removing 3 closest sites
    \\LKOthreeRate        — A+ rate (%) after removing 3 closest
    \\LKOthreeP           — p-value after removing 3 closest

  Worst-case across ALL C(11,2) and C(11,3) combinations:
    \\LKOworstTwoP        — worst p-value across all C(11,2)=55 leave-2-out combos
    \\LKOworstThreeP      — worst p-value across all C(11,3)=165 leave-3-out combos

USAGE
-----
    cd /path/to/gerizim-paper-a
    python3 analysis/unesco/dome_leave_k_out.py

REPRODUCIBILITY
---------------
  Deterministic — exact binomial, no random seed.
  Keyword list from lib/dome_filter.FORM_KEYWORD_RES (keywords.json).
  Tier thresholds from lib/beru.py (config.json).
"""

from __future__ import annotations

__version__ = "1.0.0"

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from data.unesco_corpus import load_corpus
from lib.beru import BERU, TIER_APLUS, P_NULL_AP, deviation, tier_label, is_aplus
from lib.dome_filter import FORM_KEYWORD_RES
from lib.results_store import ResultsStore
from scipy.stats import binomtest


# ── Short-name helpers ────────────────────────────────────────────────────────

def _short_name(full: str, max_len: int = 40) -> str:
    """Return a concise site name suitable for a LaTeX macro value."""
    # Strip parenthetical suffixes that bulk up the name
    for sep in (" (", " and the ", " and its"):
        if sep in full:
            full = full[: full.index(sep)]
    return full.strip()[:max_len]


# ── Load corpus and classify dome sites ──────────────────────────────────────

def _load_dome_sites():
    corpus = load_corpus()
    sites = []
    for site in corpus:
        if site.category == "Natural":
            continue
        if not site.has_coords:
            continue
        txt = site.full_text
        if not any(r.search(txt) for r in FORM_KEYWORD_RES.values()):
            continue
        dev = deviation(site.longitude)
        dev_km = dev * BERU * 111.0          # equatorial km per degree
        ap = is_aplus(tier_label(dev))
        sites.append({
            "name":    site.site,
            "lon":     site.longitude,
            "dev":     dev,
            "dev_km":  dev_km,
            "ap":      ap,
        })
    return sites


# ── Leave-k-out core ──────────────────────────────────────────────────────────

def _lko_pvalue(all_sites: list[dict], removed_indices: set[int]) -> tuple[int, int, float, float]:
    remaining = [s for i, s in enumerate(all_sites) if i not in removed_indices]
    n   = len(remaining)
    nap = sum(1 for s in remaining if s["ap"])
    p   = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue
    rate = 100.0 * nap / n if n else 0.0
    return n, nap, rate, p


def main():
    sites = _load_dome_sites()
    n_total = len(sites)
    ap_sites = [(i, s) for i, s in enumerate(sites) if s["ap"]]
    n_ap = len(ap_sites)

    # Sort A+ sites by deviation (ascending = closest to harmonic)
    ap_sites_sorted = sorted(ap_sites, key=lambda x: x[1]["dev"])

    print("=" * 72)
    print(f"  dome_leave_k_out.py  v{__version__}")
    print(f"  Dome corpus: N={n_total}, A+={n_ap}")
    print("=" * 72)

    print(f"\n  A+ sites ranked by closeness to harmonic node (beru deviation):")
    print(f"  {'Rank':>4}  {'Dev (km)':>9}  {'Name'}")
    print(f"  {'----':>4}  {'---------':>9}  {'----'}")
    for rank, (idx, s) in enumerate(ap_sites_sorted, 1):
        print(f"  {rank:>4}  {s['dev_km']:>9.2f}  {s['name']}")

    # ── Top 3 closest A+ sites ────────────────────────────────────────────────
    top3 = ap_sites_sorted[:3]
    site_one_idx,   site_one   = top3[0]
    site_two_idx,   site_two   = top3[1]
    site_three_idx, site_three = top3[2]

    # ── Leave-2-out: remove sites ranked 1 & 2 ───────────────────────────────
    lko2_n, lko2_ap, lko2_rate, lko2_p = _lko_pvalue(
        sites, {site_one_idx, site_two_idx}
    )

    # ── Leave-3-out: remove sites ranked 1, 2, 3 ─────────────────────────────
    lko3_n, lko3_ap, lko3_rate, lko3_p = _lko_pvalue(
        sites, {site_one_idx, site_two_idx, site_three_idx}
    )

    print(f"\n  Leave-2-out (drop rank 1 + rank 2):")
    print(f"    N={lko2_n}, A+={lko2_ap}, rate={lko2_rate:.1f}%, p={lko2_p:.4f}")
    print(f"\n  Leave-3-out (drop rank 1 + rank 2 + rank 3):")
    print(f"    N={lko3_n}, A+={lko3_ap}, rate={lko3_rate:.1f}%, p={lko3_p:.4f}")

    # ── Worst-case across ALL C(n_ap, 2) and C(n_ap, 3) combos ──────────────
    all_ap_indices = [i for i, _ in ap_sites]

    worst2_p = 0.0
    for combo in itertools.combinations(all_ap_indices, 2):
        _, _, _, p = _lko_pvalue(sites, set(combo))
        if p > worst2_p:
            worst2_p = p

    worst3_p = 0.0
    for combo in itertools.combinations(all_ap_indices, 3):
        _, _, _, p = _lko_pvalue(sites, set(combo))
        if p > worst3_p:
            worst3_p = p

    n_combos2 = sum(1 for _ in itertools.combinations(all_ap_indices, 2))
    n_combos3 = sum(1 for _ in itertools.combinations(all_ap_indices, 3))
    print(f"\n  Worst-case leave-2-out  (all {n_combos2} pairs):   p = {worst2_p:.4f}")
    print(f"  Worst-case leave-3-out  (all {n_combos3} triples): p = {worst3_p:.4f}")

    # ── LaTeX macros ──────────────────────────────────────────────────────────
    def fmt_p(p: float) -> str:
        if p < 0.001:
            return "$< 0.001$"
        return f"{p:.4f}"

    print(f"\n  % LaTeX macros (GROUP LOO2):")
    # Site names
    print(f"  \\newcommand{{\\LKOsiteOneName}}{{{_short_name(site_one['name'])}}}  % LKO: closest A+ site name")
    print(f"  \\newcommand{{\\LKOsiteTwoName}}{{{_short_name(site_two['name'])}}}  % LKO: 2nd-closest A+ site name")
    print(f"  \\newcommand{{\\LKOsiteThreeName}}{{{_short_name(site_three['name'])}}}  % LKO: 3rd-closest A+ site name")
    # Site deviations in km
    print(f"  \\newcommand{{\\LKOsiteOneDevKm}}{{{site_one['dev_km']:.2f}}}  % LKO: closest A+ site deviation (km)")
    print(f"  \\newcommand{{\\LKOsiteTwoDevKm}}{{{site_two['dev_km']:.2f}}}  % LKO: 2nd-closest A+ site deviation (km)")
    print(f"  \\newcommand{{\\LKOsiteThreeDevKm}}{{{site_three['dev_km']:.2f}}}  % LKO: 3rd-closest A+ site deviation (km)")
    # Leave-2-out stats
    print(f"  \\newcommand{{\\LKOtwoN}}{{{lko2_n}}}  % LKO2: N after removing 2 closest")
    print(f"  \\newcommand{{\\LKOtwoAp}}{{{lko2_ap}}}  % LKO2: A+ after removing 2 closest")
    print(f"  \\newcommand{{\\LKOtwoRate}}{{{lko2_rate:.1f}}}  % LKO2: A+ rate (%) after removing 2 closest")
    print(f"  \\newcommand{{\\LKOtwoP}}{{{fmt_p(lko2_p)}}}  % LKO2: p-value after removing 2 closest")
    # Leave-3-out stats
    print(f"  \\newcommand{{\\LKOthreeN}}{{{lko3_n}}}  % LKO3: N after removing 3 closest")
    print(f"  \\newcommand{{\\LKOthreeAp}}{{{lko3_ap}}}  % LKO3: A+ after removing 3 closest")
    print(f"  \\newcommand{{\\LKOthreeRate}}{{{lko3_rate:.1f}}}  % LKO3: A+ rate (%) after removing 3 closest")
    print(f"  \\newcommand{{\\LKOthreeP}}{{{fmt_p(lko3_p)}}}  % LKO3: p-value after removing 3 closest")
    # Worst-case across all combos
    print(f"  \\newcommand{{\\LKOworstTwoP}}{{{fmt_p(worst2_p)}}}  % LKO: worst p across all C({n_ap},2)={n_combos2} leave-2-out combos")
    print(f"  \\newcommand{{\\LKOworstThreeP}}{{{fmt_p(worst3_p)}}}  % LKO: worst p across all C({n_ap},3)={n_combos3} leave-3-out combos")
    # Combinatorial counts (so manuscript math can reference them)
    print(f"  \\newcommand{{\\LKOnAp}}{{{n_ap}}}  % LKO: total A+ sites in dome corpus (n for C(n,k))")
    print(f"  \\newcommand{{\\LKOnCombosTwo}}{{{n_combos2}}}  % LKO: C({n_ap},2) leave-2-out combinations")
    print(f"  \\newcommand{{\\LKOnCombosThree}}{{{n_combos3}}}  % LKO: C({n_ap},3) leave-3-out combinations")

    print(f"\n" + "=" * 72)
    print(f"  DONE — all GROUP LOO2 macros printed above.")
    print("=" * 72)

    # ── Write to results store ────────────────────────────────────────────────
    store = ResultsStore()
    store.write_many({
        "LKOtwoN":         lko2_n,
        "LKOtwoAp":        lko2_ap,
        "LKOtwoRate":      round(lko2_rate, 1),
        "LKOtwoP":         round(lko2_p, 4),
        "LKOthreeN":       lko3_n,
        "LKOthreeAp":      lko3_ap,
        "LKOthreeRate":    round(lko3_rate, 1),
        "LKOthreeP":       round(lko3_p, 4),
        "LKOworstTwoP":    round(worst2_p, 4),
        "LKOworstThreeP":  round(worst3_p, 4),
        "LKOnAp":          n_ap,
        "LKOnCombosTwo":   n_combos2,
        "LKOnCombosThree": n_combos3,
    })
    print("\n  ✓ Written to data/store/results.json")


if __name__ == "__main__":
    main()
