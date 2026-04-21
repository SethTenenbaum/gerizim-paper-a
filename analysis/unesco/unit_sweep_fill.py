import re
import sys
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, TIER_APLUS, deviation_at_spacing
from lib.stats import significance_label as _sig_label
from lib.results_store import ResultsStore

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

    # ── Overlap between 0.08 and 0.10 spacings ──────────────────────────────
    def hit_mask(lons, spacing):
        return np.array([deviation_at_spacing(lon, spacing) <= TIER_APLUS for lon in lons], dtype=bool)

    dome_08_mask = hit_mask(dome_lons, 0.08)
    dome_10_mask = hit_mask(dome_lons, 0.10)
    full_08_mask = hit_mask(all_lons, 0.08)
    full_10_mask = hit_mask(all_lons, 0.10)

    dome_08 = int(np.count_nonzero(dome_08_mask))
    dome_10 = int(np.count_nonzero(dome_10_mask))
    dome_both = int(np.count_nonzero(dome_08_mask & dome_10_mask))
    full_08 = int(np.count_nonzero(full_08_mask))
    full_10 = int(np.count_nonzero(full_10_mask))
    full_both = int(np.count_nonzero(full_08_mask & full_10_mask))

    print()
    print("=" * 100)
    print("  OVERLAP BETWEEN 0.08 AND 0.10 SPACINGS")
    print("=" * 100)
    print(f"{'Set':>8} {'Hits@0.08':>10} {'Hits@0.10':>10} {'Hits@Both':>10} {'%Both of 0.08':>14} {'%Both of 0.10':>14}")
    print("-" * 100)
    print(f"  {'DOME':>6} {dome_08:>10} {dome_10:>10} {dome_both:>10} {100*dome_both/dome_08 if dome_08 else 0:14.1f}% {100*dome_both/dome_10 if dome_10 else 0:14.1f}%")
    print(f"  {'FULL':>6} {full_08:>10} {full_10:>10} {full_both:>10} {100*full_both/full_08 if full_08 else 0:14.1f}% {100*full_both/full_10 if full_10 else 0:14.1f}%")

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

    # ── LaTeX macros for adjacent-spacing prose references ─────────────────
    # The 0.08-beru row is cited in the robustness section to show that
    # non-harmonic adjacent spacings do not produce joint significance.
    print()
    print("  % LaTeX macros (unit_sweep_fill.py):")

    # Map spacing -> CamelCase tag used in macro names
    COARSE_TAG = {
        0.05: "ZeroFive",   0.06: "ZeroSix",   0.07: "ZeroSeven",
        0.08: "ZeroEight",  0.09: "ZeroNine",  0.10: "ZeroTen",
        0.11: "ZeroEleven", 0.12: "ZeroTwelve",0.15: "ZeroFifteen",
        0.20: "ZeroTwenty",
    }
    FINE_TAG = {
        0.0990: "NineNineZero",  0.0993: "NineNineThree",
        0.0995: "NineNineFive",  0.1001: "OneZeroZeroOne",
        0.1003: "OneZeroZeroThree", 0.1010: "OneZeroOneZero",
    }

    def _fmt(p):
        return f"{p:.4f}" if p >= 0.0001 else "<0.0001"

    for sp, nr, d_hits, d_p, f_hits, f_p in coarse_rows:
        tag = COARSE_TAG.get(round(sp, 4), "")
        if not tag or sp == 0.10:   # canonical row uses named macros elsewhere
            if abs(sp - 0.08) < 1e-9:   # keep backward-compat names
                print(f"  \\newcommand{{\\sweepEightDomeP}}{{{_fmt(d_p)}}}   % p(Dome), 0.08-beru spacing")
                print(f"  \\newcommand{{\\sweepEightFullP}}{{{_fmt(f_p)}}}   % p(Full), 0.08-beru spacing")
            continue
        print(f"  \\newcommand{{\\sweep{tag}Null}}{{{nr*100:.1f}\\%}}   % null rate, {sp:.2f}-beru")
        print(f"  \\newcommand{{\\sweep{tag}DomeHits}}{{{d_hits}}}   % dome hits, {sp:.2f}-beru")
        print(f"  \\newcommand{{\\sweep{tag}DomeP}}{{{_fmt(d_p)}}}   % p(Dome), {sp:.2f}-beru")
        print(f"  \\newcommand{{\\sweep{tag}FullHits}}{{{f_hits}}}   % full hits, {sp:.2f}-beru")
        print(f"  \\newcommand{{\\sweep{tag}FullP}}{{{_fmt(f_p)}}}   % p(Full), {sp:.2f}-beru")

    print(f"  \\newcommand{{\\sweepEightDomeP}}{{{_fmt(next(d_p for sp,nr,d_hits,d_p,f_hits,f_p in coarse_rows if abs(sp-0.08)<1e-9))}}}   % p(Dome), 0.08-beru")
    print(f"  \\newcommand{{\\sweepEightFullP}}{{{_fmt(next(f_p for sp,nr,d_hits,d_p,f_hits,f_p in coarse_rows if abs(sp-0.08)<1e-9))}}}   % p(Full), 0.08-beru")
    # Persist sweepEight p-values so emit_sig_macros.py can emit their Sig companions
    _eight_d_p = next(d_p for sp,nr,d_hits,d_p,f_hits,f_p in coarse_rows if abs(sp-0.08)<1e-9)
    _eight_f_p = next(f_p for sp,nr,d_hits,d_p,f_hits,f_p in coarse_rows if abs(sp-0.08)<1e-9)
    ResultsStore().write_many({"sweepEightDomeP": _eight_d_p, "sweepEightFullP": _eight_f_p})
    # Overlap macros: 0.08 vs 0.10 spacings
    dome_08_pct_of_ten = round(100 * dome_both / dome_10, 1) if dome_10 else 0.0
    dome_08_pct_of_eight = round(100 * dome_both / dome_08, 1) if dome_08 else 0.0
    full_08_pct_of_ten = round(100 * full_both / full_10, 1) if full_10 else 0.0
    full_08_pct_of_eight = round(100 * full_both / full_08, 1) if full_08 else 0.0
    print(f"  \\newcommand{{\\sweepEightTenDomeBoth}}{{{dome_both}}}   % DOME hits shared by 0.08- and 0.10-beru")
    print(f"  \\newcommand{{\\sweepEightTenDomePctOfTen}}{{{dome_08_pct_of_ten}}}   % pct of 0.10-beru dome hits also in 0.08-beru")
    print(f"  \\newcommand{{\\sweepEightTenDomePctOfEight}}{{{dome_08_pct_of_eight}}}   % pct of 0.08-beru dome hits also in 0.10-beru")
    print(f"  \\newcommand{{\\sweepEightTenFullBoth}}{{{full_both}}}   % full-corpus hits shared by 0.08- and 0.10-beru")
    print(f"  \\newcommand{{\\sweepEightTenFullPctOfTen}}{{{full_08_pct_of_ten}}}   % pct of 0.10-beru full hits also in 0.08-beru")
    print(f"  \\newcommand{{\\sweepEightTenFullPctOfEight}}{{{full_08_pct_of_eight}}}   % pct of 0.08-beru full hits also in 0.10-beru")

    for sp, nr, d_hits, d_p, f_hits, f_p in fine_rows:
        tag = FINE_TAG.get(round(sp, 4), "")
        if not tag:
            continue
        print(f"  \\newcommand{{\\fine{tag}Null}}{{{nr*100:.2f}\%}}   % null rate, {sp:.4f}-beru")
        print(f"  \\newcommand{{\\fine{tag}DomeHits}}{{{d_hits}}}   % dome hits, {sp:.4f}-beru")
        print(f"  \\newcommand{{\\fine{tag}DomeP}}{{{_fmt(d_p)}}}   % p(Dome), {sp:.4f}-beru")
        print(f"  \\newcommand{{\\fine{tag}FullHits}}{{{f_hits}}}   % full hits, {sp:.4f}-beru")
        print(f"  \\newcommand{{\\fine{tag}FullP}}{{{_fmt(f_p)}}}   % p(Full), {sp:.4f}-beru")

if __name__ == "__main__":
    main()
