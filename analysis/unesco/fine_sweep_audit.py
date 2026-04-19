"""
Fine unit sweep audit: ±1% of 0.10 beru.

Reports A+ enrichment for dome and full-corpus populations at each
spacing step. The canonical 0.10-beru row is highlighted. A ±0.3%
shift collapses joint significance in at least one population.

Replaces the former Online Supplementary Table S3 (removed when the
coarse unit-sweep table was consolidated into prose).

Usage:
    python analysis/unesco/fine_sweep_audit.py
"""

import sys
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import TIER_APLUS, deviation_at_spacing
from lib.dome_filter import is_dome_site
from lib.stats import significance_label as _sig_label, null_rate_at_spacing

FINE_SPACINGS = [
    0.0970, 0.0980, 0.0990, 0.0993, 0.0995,
    0.1000,
    0.1001, 0.1003, 0.1010, 0.1020, 0.1030,
]


def sig(p):
    label = _sig_label(p)
    return "~" if label == "~" else label


def count_hits(lons, spacing):
    return sum(1 for lon in lons if deviation_at_spacing(lon, spacing) <= TIER_APLUS)


def main():
    corpus = load_corpus()
    sites = cultural_sites_with_coords(corpus)

    all_lons = [s.longitude for s in sites]
    dome_lons = [s.longitude for s in sites if is_dome_site(s)]

    N_all = len(all_lons)
    N_dome = len(dome_lons)

    header = f"{'Spacing':>8}  {'Null%':>6}  {'Dome A+':>8}  {'Dome p':>10}  {'Sig':>4}  {'Full A+':>8}  {'Full p':>10}  {'Sig':>4}"
    sep = "-" * len(header)

    print(f"Fine unit sweep (±1% of canonical 0.10 beru)")
    print(f"Dome N = {N_dome}, Full N = {N_all}")
    print(sep)
    print(header)
    print(sep)

    for sp in FINE_SPACINGS:
        nr = null_rate_at_spacing(TIER_APLUS, sp)

        d_hits = count_hits(dome_lons, sp)
        d_p = binomtest(d_hits, N_dome, nr, alternative="greater").pvalue

        f_hits = count_hits(all_lons, sp)
        f_p = binomtest(f_hits, N_all, nr, alternative="greater").pvalue

        marker = ">>>" if abs(sp - 0.100) < 1e-6 else "   "
        d_p_str = f"{d_p:.4f}" if d_p >= 0.0001 else f"{d_p:.2e}"
        f_p_str = f"{f_p:.4f}" if f_p >= 0.0001 else f"{f_p:.2e}"
        print(
            f"{marker} {sp:.4f}  {nr*100:5.2f}%  {d_hits:>3}/{N_dome}   {d_p_str:>10}  {sig(d_p):>4}  "
            f"{f_hits:>3}/{N_all}   {f_p_str:>10}  {sig(f_p):>4}"
        )

    print(sep)
    print(">>> = canonical 0.10-beru row")
    print("Joint significance requires both dome and full-corpus p < 0.05.")
    print("A ±0.3% shift from canonical collapses joint significance in at least one population.")


if __name__ == "__main__":
    main()
