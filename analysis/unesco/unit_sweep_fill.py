
import re
import sys
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, TIER_APLUS, deviation_at_spacing
from lib.stats import significance_label as _sig_label

# shared with spherical_monument_test.py)
from lib.dome_filter import is_dome_site
DOME_EXCLUDE = set()  # no exclusions (see spherical_monument_test.py for rationale)

COARSE_SPACINGS = [0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.15, 0.20]

# Fine spacings: ±1% of 0.10
FINE_SPACINGS = [0.099, 0.0993, 0.0995, 0.100, 0.1001, 0.1003, 0.101]

def load_sites():
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)

    sites = []
    for s in cultural:
        # consistent with spherical_monument_test.py
        is_dome = is_dome_site(s) if s.site not in DOME_EXCLUDE else False

        sites.append({"name": s.site, "lon": s.longitude, "is_dome": is_dome})

    return sites

def count_hits(lons, spacing):
    hits = 0
    for lon in lons:
        dev = deviation_at_spacing(lon, spacing)
        if dev <= TIER_APLUS:
            hits += 1
    return hits

def null_rate(spacing):
    from lib.stats import null_rate_at_spacing
    return null_rate_at_spacing(TIER_APLUS, spacing)

def sig_label(p):
    label = _sig_label(p)
    return "\\sim" if label == "~" else label

def main():
    sites = load_sites()
    all_lons = [s["lon"] for s in sites]
    dome_lons = [s["lon"] for s in sites if s["is_dome"]]

    N_all = len(all_lons)
    N_dome = len(dome_lons)

    print(f"Loaded {N_all} Cultural/Mixed sites, {N_dome} dome/spherical sites")
    print()

    # ── Coarse sweep ─────────────────────────────────────────────────────────
    print("=" * 100)
    print("  COARSE UNIT SWEEP (Table 6 in manuscript)")
    print("=" * 100)
    print(f"{'Spacing':>8} {'Null%':>6} {'Dome hits':>10} {'Dome N':>7} {'p(Dome)':>12} {'Sig':>5} "
          f"{'Full hits':>10} {'Full N':>7} {'p(Full)':>12} {'Sig':>5}")
    print("-" * 100)

    coarse_rows = []
    for sp in COARSE_SPACINGS:
        nr = null_rate(sp)

        # Dome
        d_hits = count_hits(dome_lons, sp)
        d_res = binomtest(d_hits, N_dome, nr, alternative='greater')
        d_p = d_res.pvalue

        # Full corpus
        f_hits = count_hits(all_lons, sp)
        f_res = binomtest(f_hits, N_all, nr, alternative='greater')
        f_p = f_res.pvalue

        bold = ">>> " if sp == 0.10 else "    "
        print(f"{bold}{sp:.3f}  {nr*100:5.1f}%  {d_hits:>10}/{N_dome}  {d_p:12.4f} {sig_label(d_p):>5}  "
              f"{f_hits:>10}/{N_all}  {f_p:12.4f} {sig_label(f_p):>5}")

        coarse_rows.append((sp, nr, d_hits, d_p, f_hits, f_p))

    # ── Fine sweep ───────────────────────────────────────────────────────────
    print()
    print("=" * 100)
    print("  FINE UNIT SWEEP (Table 7 in manuscript)")
    print("=" * 100)
    print(f"{'Spacing':>8} {'Null%':>6} {'Dome hits':>10} {'Dome N':>7} {'p(Dome)':>12} {'Sig':>5} "
          f"{'Full hits':>10} {'Full N':>7} {'p(Full)':>12} {'Sig':>5}")
    print("-" * 100)

    fine_rows = []
    for sp in FINE_SPACINGS:
        nr = null_rate(sp)

        d_hits = count_hits(dome_lons, sp)
        d_res = binomtest(d_hits, N_dome, nr, alternative='greater')
        d_p = d_res.pvalue

        f_hits = count_hits(all_lons, sp)
        f_res = binomtest(f_hits, N_all, nr, alternative='greater')
        f_p = f_res.pvalue

        bold = ">>> " if sp == 0.100 else "    "
        print(f"{bold}{sp:.4f} {nr*100:5.2f}%  {d_hits:>10}/{N_dome}  {d_p:12.4f} {sig_label(d_p):>5}  "
              f"{f_hits:>10}/{N_all}  {f_p:12.4f} {sig_label(f_p):>5}")

        fine_rows.append((sp, nr, d_hits, d_p, f_hits, f_p))

    # ── LaTeX output ─────────────────────────────────────────────────────────
    print()
    print("=" * 100)
    print("  LATEX TABLE ROWS (COARSE)")
    print("=" * 100)
    for sp, nr, d_hits, d_p, f_hits, f_p in coarse_rows:
        bold_s = "\\textbf{" if sp == 0.10 else ""
        bold_e = "}" if sp == 0.10 else ""
        d_sig = sig_label(d_p)
        f_sig = sig_label(f_p)
        d_str = f"{d_hits}/{N_dome}"
        f_str = f"{f_hits}/{N_all}"
        d_p_str = f"{d_p:.4f}" if d_p >= 0.0001 else f"{d_p:.1e}"
        f_p_str = f"{f_p:.4f}" if f_p >= 0.0001 else f"{f_p:.1e}"
        print(f"{bold_s}{sp:.2f}{bold_e} & {bold_s}{nr*100:.1f}\\%{bold_e} "
              f"& {bold_s}{d_str}{bold_e} & {bold_s}{d_p_str}~{d_sig}{bold_e} "
              f"& {bold_s}{f_str}{bold_e} & {bold_s}{f_p_str}~{f_sig}{bold_e} \\\\")

    print()
    print("  LATEX TABLE ROWS (FINE)")
    for sp, nr, d_hits, d_p, f_hits, f_p in fine_rows:
        bold_s = "\\textbf{" if sp == 0.100 else ""
        bold_e = "}" if sp == 0.100 else ""
        d_sig = sig_label(d_p)
        f_sig = sig_label(f_p)
        d_str = f"{d_hits}/{N_dome}"
        f_str = f"{f_hits}/{N_all}"
        d_p_str = f"{d_p:.4f}" if d_p >= 0.0001 else f"{d_p:.1e}"
        f_p_str = f"{f_p:.4f}" if f_p >= 0.0001 else f"{f_p:.1e}"
        print(f"{bold_s}{sp:.4f}{bold_e} & {bold_s}{nr*100:.2f}\\%{bold_e} "
              f"& {bold_s}{d_str}{bold_e} & {bold_s}{d_p_str}~{d_sig}{bold_e} "
              f"& {bold_s}{f_str}{bold_e} & {bold_s}{f_p_str}~{f_sig}{bold_e} \\\\")

if __name__ == "__main__":
    main()
