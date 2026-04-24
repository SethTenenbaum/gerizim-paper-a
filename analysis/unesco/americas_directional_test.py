"""
americas_directional_test.py
==============================
Formal one-sided directional test of beru-harmonic depletion in the
UNESCO Americas sub-corpus.

MOTIVATION
----------
The beru-alignment hypothesis predicts that A+ rates should be lower in
non-Eurasian sub-corpora — regions with no documented connection to
Babylonian metrology.  The Americas is the natural negative-control
region: the most geographically remote from the Levantine corridor.

This script frames that prediction as a pre-specified one-sided
directional test and computes the formal p-value for the paper's
Limitations section.

DESIGN
------
Population: all UNESCO Cultural/Mixed sites with coordinates whose
  UNESCO regional tag contains "Latin America and the Caribbean" or
  "North America", OR whose longitude falls in the Western Hemisphere
  (lon < -30°), whichever set is larger.  (The regional-tag approach
  is primary; the longitude heuristic is the fallback if region tags
  are unavailable.)

H₀: P(Tier-A+) = 4%  (geometric null)
H₁: P(Tier-A+) < 4%  (one-sided depletion, the directional prediction)

Test: one-sided exact binomial (scipy.stats.binom_test or binomtest).

USAGE
-----
    python3 analysis/unesco/americas_directional_test.py

MANUSCRIPT MACROS PRODUCED
---------------------------
    \\AmericasApCount        — observed A+ in Americas sub-corpus
    \\AmericasApRate         — A+ rate (%) in Americas
    \\AmericasOneSidedP      — one-sided binomial p-value (depletion)
    \\AmericasOneSidedStar   — significance label
    \\AmericasTwoSidedP      — two-sided binomial p-value (deviation from null)
    \\AmericasTwoSidedStar   — significance label
    \\AmericasDirectional    — "in" or "opposite to" (confirms direction)
    \\AmericasN              — Americas sub-corpus N
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import TIER_APLUS, P_NULL_AP
from lib.beru import deviation as _beru_deviation
from lib.results_store import ResultsStore

# Longitude cutoff for Western Hemisphere fallback
AMERICAS_LON_MAX = -30.0

# UNESCO regional strings that identify Americas sites
AMERICAS_REGION_STRINGS = [
    "latin america",
    "caribbean",
    "north america",
]


def beru_dev(lon: float) -> float:
    return _beru_deviation(lon)


def is_americas_site(site) -> bool:
    """Return True if the site is in the Americas by region tag or longitude."""
    region = (getattr(site, "region", "") or "").lower()
    if any(tag in region for tag in AMERICAS_REGION_STRINGS):
        return True
    # Fallback: Western-Hemisphere longitude
    lon = getattr(site, "longitude", None)
    if lon is not None and lon < AMERICAS_LON_MAX:
        return True
    return False


def sig_label(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    if p < 0.10:
        return "~"
    return "ns"


def binomial_p_less(n_ap: int, n_total: int) -> float:
    """One-sided binomial p-value for depletion (H₁: rate < P_NULL_AP)."""
    try:
        from scipy.stats import binom_test
        return float(binom_test(n_ap, n_total, P_NULL_AP, alternative="less"))
    except (TypeError, ImportError):
        from scipy.stats import binomtest as _bt
        return float(_bt(n_ap, n_total, P_NULL_AP, alternative="less").pvalue)


def binomial_p_two_sided(n_ap: int, n_total: int) -> float:
    """Two-sided binomial p-value for deviation from the null."""
    try:
        from scipy.stats import binom_test
        return float(binom_test(n_ap, n_total, P_NULL_AP, alternative="two-sided"))
    except (TypeError, ImportError):
        from scipy.stats import binomtest as _bt
        return float(_bt(n_ap, n_total, P_NULL_AP, alternative="two-sided").pvalue)


def main():
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)

    americas_sites = [s for s in cultural if is_americas_site(s)]
    N_americas = len(americas_sites)
    americas_lons = [s.longitude for s in americas_sites]

    n_ap = sum(1 for lon in americas_lons if beru_dev(lon) <= TIER_APLUS)
    rate = 100.0 * n_ap / N_americas if N_americas > 0 else 0.0
    p_one_sided = binomial_p_less(n_ap, N_americas)
    p_two_sided = binomial_p_two_sided(n_ap, N_americas)
    label_one_sided = sig_label(p_one_sided)
    label_two_sided = sig_label(p_two_sided)

    directional = "in" if rate < 100 * P_NULL_AP else "opposite to"

    print("=" * 80)
    print("  AMERICAS UNESCO SUB-CORPUS — BINOMIAL TESTS")
    print(f"  H₀: P(A+) = {100*P_NULL_AP:.0f}%")
    print("=" * 80)
    print(f"\n  Americas sub-corpus N = {N_americas}")
    print(f"  Observed A+ = {n_ap}  ({rate:.1f}%)")
    print(f"  One-sided binomial p (depletion) = {p_one_sided:.4f}  {label_one_sided}")
    print(f"  Two-sided binomial p = {p_two_sided:.4f}  {label_two_sided}")
    print(f"  Direction: {directional} the pre-specified direction")

    print(f"\n  INTERPRETATION:")
    print(f"    The geometric null predicts {100*P_NULL_AP:.0f}% A+; the Americas shows {rate:.1f}%.")
    if rate < 100 * P_NULL_AP:
        print(f"    The observed rate is below the null, but not significantly so in the conservative")
        print(f"    two-sided test (p = {p_two_sided:.4f}, {label_two_sided}).")
    else:
        print(f"    The observed rate is not below the null, and the two-sided test is not significant")
        print(f"    (p = {p_two_sided:.4f}, {label_two_sided}).")

    print("\n" + "=" * 80)
    print("  LATEX MACROS (Americas binomial tests):")
    print("=" * 80)
    print(f"  \\newcommand{{\\AmericasN}}{{{N_americas}}}")
    print(f"  \\newcommand{{\\AmericasApCount}}{{{n_ap}}}")
    print(f"  \\newcommand{{\\AmericasApRate}}{{{rate:.1f}}}")
    print(f"  \\newcommand{{\\AmericasOneSidedP}}{{{p_one_sided:.4f}}}")
    print(f"  \\newcommand{{\\AmericasOneSidedStar}}{{{label_one_sided}}}")
    print(f"  \\newcommand{{\\AmericasTwoSidedP}}{{{p_two_sided:.4f}}}")
    print(f"  \\newcommand{{\\AmericasTwoSidedStar}}{{{label_two_sided}}}")
    print(f"  \\newcommand{{\\AmericasDirectional}}{{{directional}}}")

    ResultsStore().write_many({
        "AmericasN":            float(N_americas),
        "AmericasApCount":      float(n_ap),
        "AmericasApRate":       round(rate, 1),
        "AmericasOneSidedP":    round(p_one_sided, 4),
        "AmericasTwoSidedP":    round(p_two_sided, 4),
        "AmericasDirectional":  directional,
    })
    print("\nResults written to data/store/results.json")


if __name__ == "__main__":
    main()
