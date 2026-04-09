#!/usr/bin/env python3
"""
generate_figures.py — Generate all three manuscript figures for Paper A.

Figures produced:
    1. fig_devhist.pdf      — Histogram of beru deviations (δ) for all UNESCO sites
    2. fig_temporal.pdf      — Temporal gradient: A+ rate by inscription-year cohort
    3. fig_unitsweep.pdf     — Unit sensitivity: −log₁₀(p) vs harmonic spacing

All constants are loaded from config.json via lib/beru.py.
Data is loaded via data.unesco_corpus.

Usage:
    cd gerizim-analysis
    python3 manuscript/generate_figures.py
"""

import sys
from pathlib import Path

# ── Project root on sys.path ─────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for PDF
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.stats import binomtest

from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_A_MAX, P_NULL_AP,
    deviation as beru_dev, tier_label as tier, deviation_at_spacing,
)
from lib.stats import (
    significance_label as sig, cochran_armitage, null_rate_at_spacing,
)
from lib.dome_filter import is_dome_site

# ── Matplotlib house style ───────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.08,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

OUTDIR = Path(__file__).resolve().parent / "figures"
OUTDIR.mkdir(exist_ok=True)

# Colour palette — muted, colour-blind friendly
C_PRIMARY   = "#2c6fbb"  # deep blue
C_HIGHLIGHT = "#d62728"  # red accent for shaded tier
C_NULL      = "#888888"  # grey for null baseline
C_DOME      = "#e6550d"  # orange for dome sub-population
C_FULL      = "#2c6fbb"  # blue for full corpus
C_COHORT    = ["#1b7837", "#5aae61", "#a6dba0", "#d9f0d3", "#e7e1ef"]
C_BAR_SIG   = "#2c6fbb"  # significant cohort bar
C_BAR_NS    = "#bdbdbd"  # non-significant cohort bar


# ══════════════════════════════════════════════════════════════════════════════
#  Load data once
# ══════════════════════════════════════════════════════════════════════════════
print("Loading UNESCO corpus…")
corpus = load_corpus()
cultural = cultural_sites_with_coords(corpus)
N_total = len(cultural)
print(f"  {N_total} Cultural/Mixed sites with coordinates")

# Compute deviations for every site
sites = []
for s in cultural:
    yr = s.year
    d = beru_dev(s.longitude)
    sites.append({
        "name": s.site,
        "lon": s.longitude,
        "yr": yr,
        "dev": d,
        "tier": tier(d),
        "is_ap": d <= TIER_APLUS,
        "is_dome": is_dome_site(s),
    })

deviations = np.array([s["dev"] for s in sites])
print(f"  A+ sites: {sum(s['is_ap'] for s in sites)}")
print()


# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 1: Beru deviation histogram  (fig:devhist)
# ══════════════════════════════════════════════════════════════════════════════
def make_devhist():
    """Histogram of beru deviations with A+ window shaded."""
    fig, ax = plt.subplots(figsize=(7, 4.2))

    bin_width = 0.005
    max_dev = 0.05
    bins = np.arange(0, max_dev + bin_width, bin_width)
    n_bins = len(bins) - 1

    counts, edges, patches = ax.hist(
        deviations, bins=bins,
        color=C_PRIMARY, edgecolor="white", linewidth=0.6,
        alpha=0.85, zorder=3,
    )

    # Shade the A+ window (first bin: 0–0.005 covers δ ≤ 0.002 and a bit more)
    # More precisely: highlight only δ ≤ 0.002 region
    ax.axvspan(0, TIER_APLUS, color=C_HIGHLIGHT, alpha=0.12, zorder=1,
               label=f"Tier-A+ window (δ ≤ {TIER_APLUS})")

    # Uniform null expectation line
    expected_per_bin = N_total / n_bins
    ax.axhline(expected_per_bin, color=C_NULL, linewidth=1.5,
               linestyle="--", zorder=4, label=f"Uniform null ({expected_per_bin:.1f}/bin)")

    # Annotate the leftmost bin — point inward from the right to avoid
    # colliding with the title; keep arrow tip at the bar top
    leftmost_count = int(counts[0])
    y_max_data = counts.max()
    ax.annotate(
        f"{leftmost_count} sites",
        xy=(bin_width / 2, leftmost_count),
        xytext=(0.018, leftmost_count - 18),   # below the bar top, inside the plot
        fontsize=9, fontweight="bold", color=C_HIGHLIGHT,
        arrowprops=dict(arrowstyle="->", color=C_HIGHLIGHT, lw=1.2),
        zorder=5,
        annotation_clip=False,
    )

    # A+ threshold line — label placed low so it never reaches the title
    ax.axvline(TIER_APLUS, color=C_HIGHLIGHT, linewidth=1.0,
               linestyle=":", alpha=0.7, zorder=4)
    ax.text(TIER_APLUS + 0.0005, y_max_data * 0.45,
            f"δ = {TIER_APLUS}", fontsize=8, color=C_HIGHLIGHT,
            rotation=90, va="center")

    ax.set_xlabel("Beru deviation (δ) from nearest 0.1-beru harmonic")
    ax.set_ylabel("Number of sites")
    ax.set_title(f"Distribution of beru deviations — UNESCO Cultural/Mixed (N = {N_total})")
    ax.legend(loc="upper right", framealpha=0.9)
    ax.set_xlim(0, max_dev)
    # Add 20 % headroom above tallest bar so annotations never clip into the title
    ax.set_ylim(0, y_max_data * 1.20)

    fig.tight_layout()
    outpath = OUTDIR / "fig_devhist.pdf"
    fig.savefig(outpath)
    fig.savefig(outpath.with_suffix(".png"))
    plt.close(fig)
    print(f"  ✓ {outpath.name} ({outpath.with_suffix('.png').name})")
    return outpath


# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 2: Temporal gradient  (fig:temporal)
# ══════════════════════════════════════════════════════════════════════════════
def make_temporal():
    """Bar chart of A+ rate by inscription-year cohort."""
    fig, ax = plt.subplots(figsize=(7, 4.2))

    # Five cohorts matching the manuscript
    cohorts = [
        ("1978–1984", 1978, 1984),
        ("1985–1991", 1985, 1991),
        ("1992–1999", 1992, 1999),
        ("2000–2009", 2000, 2009),
        ("2010–2025", 2010, 2025),
    ]

    labels = []
    rates = []
    ns = []
    naps = []
    pvals = []

    # Only sites with a year
    dated = [s for s in sites if s["yr"] is not None]

    for label, y0, y1 in cohorts:
        subset = [s for s in dated if y0 <= s["yr"] <= y1]
        n = len(subset)
        nap = sum(1 for s in subset if s["is_ap"])
        rate = 100.0 * nap / n if n else 0
        p = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue if n else 1.0
        labels.append(label)
        rates.append(rate)
        ns.append(n)
        naps.append(nap)
        pvals.append(p)

    x = np.arange(len(labels))
    bar_colors = [C_BAR_SIG if p < 0.05 else C_BAR_NS for p in pvals]

    bars = ax.bar(x, rates, width=0.65, color=bar_colors, edgecolor="white",
                  linewidth=0.8, zorder=3)

    # 4% null baseline
    ax.axhline(100 * P_NULL_AP, color=C_NULL, linewidth=1.5, linestyle="--",
               zorder=2, label=f"Null rate ({100*P_NULL_AP:.0f}%)")

    # Annotate each bar with count and significance
    for i, (rate, n, nap, p) in enumerate(zip(rates, ns, naps, pvals)):
        sig_txt = sig(p)
        annotation = f"{nap}/{n}"
        y_off = 0.3
        ax.text(i, rate + y_off, annotation,
                ha="center", va="bottom", fontsize=8, color="#333333")
        if sig_txt not in ("ns",):
            ax.text(i, rate + y_off + 0.8, sig_txt,
                    ha="center", va="bottom", fontsize=8,
                    fontweight="bold", color=C_HIGHLIGHT)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_xlabel("UNESCO inscription-year cohort")
    ax.set_ylabel("Tier-A+ rate (%)")
    ax.set_title("Temporal gradient — A+ rate by inscription cohort")
    ax.legend(loc="upper right", framealpha=0.9)
    ax.set_ylim(0, max(rates) + 3)

    # Add Cochran-Armitage annotation
    # Three-cohort test (from manuscript macros)
    dated_for_ca = [s for s in dated if s["yr"] is not None]
    ca_cohorts_3 = [
        (1978, 1984), (1985, 1999), (2000, 2025),
    ]
    ns3, aps3 = [], []
    for y0, y1 in ca_cohorts_3:
        sub = [s for s in dated_for_ca if y0 <= s["yr"] <= y1]
        ns3.append(len(sub))
        aps3.append(sum(1 for s in sub if s["is_ap"]))
    ca3 = cochran_armitage(ns3, aps3, scores=[1, 2, 3])

    ax.text(0.98, 0.65,
            f"Cochran-Armitage (3-cohort)\n"
            f"Z = {ca3.z_statistic:.2f}, p = {ca3.p_value:.4f}",
            transform=ax.transAxes, fontsize=8,
            ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f0f0",
                      edgecolor="#cccccc", alpha=0.9))

    fig.tight_layout()
    outpath = OUTDIR / "fig_temporal.pdf"
    fig.savefig(outpath)
    fig.savefig(outpath.with_suffix(".png"))
    plt.close(fig)
    print(f"  ✓ {outpath.name} ({outpath.with_suffix('.png').name})")
    return outpath


# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 3: Unit sensitivity sweep  (fig:unitsweep)
# ══════════════════════════════════════════════════════════════════════════════
def make_unitsweep():
    """−log₁₀(p) vs harmonic spacing for dome and full corpus populations."""
    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Spacings to test
    coarse = [0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.15, 0.20]

    all_lons = np.array([s["lon"] for s in sites])
    dome_lons = np.array([s["lon"] for s in sites if s["is_dome"]])

    N_all = len(all_lons)
    N_dome = len(dome_lons)
    print(f"  Unit sweep: {N_all} full corpus, {N_dome} dome sites")

    def count_hits(lons, spacing):
        return sum(1 for lon in lons if deviation_at_spacing(lon, spacing) <= TIER_APLUS)

    # Compute p-values for each spacing
    sp_dome_logp = []
    sp_full_logp = []
    sp_dome_pval = []
    sp_full_pval = []

    for sp in coarse:
        nr = null_rate_at_spacing(TIER_APLUS, sp)

        # Dome population
        d_hits = count_hits(dome_lons, sp)
        d_p = binomtest(d_hits, N_dome, nr, alternative="greater").pvalue
        sp_dome_pval.append(d_p)
        sp_dome_logp.append(-np.log10(max(d_p, 1e-15)))

        # Full corpus
        f_hits = count_hits(all_lons, sp)
        f_p = binomtest(f_hits, N_all, nr, alternative="greater").pvalue
        sp_full_pval.append(f_p)
        sp_full_logp.append(-np.log10(max(f_p, 1e-15)))

    # Plot both series
    ax.plot(coarse, sp_dome_logp, "o-",
            color=C_DOME, linewidth=2, markersize=7,
            label=f"Dome/spherical (N = {N_dome})", zorder=4)
    ax.plot(coarse, sp_full_logp, "s-",
            color=C_FULL, linewidth=2, markersize=6,
            label=f"Full corpus (N = {N_all})", zorder=4)

    # Significance thresholds — labels anchored to right axis, inside plot
    y_max_data = max(max(sp_dome_logp), max(sp_full_logp))
    ax.axhline(-np.log10(0.05), color=C_NULL, linewidth=1, linestyle="--",
               alpha=0.7, zorder=2)
    ax.text(0.205, -np.log10(0.05), "p = 0.05",
            fontsize=7.5, color=C_NULL, ha="left", va="center",
            transform=ax.get_yaxis_transform())

    ax.axhline(-np.log10(0.01), color=C_NULL, linewidth=0.8, linestyle=":",
               alpha=0.5, zorder=2)
    ax.text(0.205, -np.log10(0.01), "p = 0.01",
            fontsize=7.5, color=C_NULL, ha="left", va="center",
            transform=ax.get_yaxis_transform())

    ax.axhline(-np.log10(0.001), color=C_NULL, linewidth=0.8, linestyle=":",
               alpha=0.4, zorder=2)
    ax.text(0.205, -np.log10(0.001), "p = 0.001",
            fontsize=7.5, color=C_NULL, ha="left", va="center",
            transform=ax.get_yaxis_transform())

    # Highlight the 0.10 peak — label placed *below* the dome series peak
    # so it never collides with the title
    ax.axvline(0.10, color="#cccccc", linewidth=12, alpha=0.3, zorder=1)
    idx_010 = coarse.index(0.10)
    dome_peak_p = sp_dome_pval[idx_010]
    full_peak_p = sp_full_pval[idx_010]
    dome_lp = sp_dome_logp[idx_010]
    full_lp = sp_full_logp[idx_010]

    # "0.10 beru" badge sits just below the dome peak, never above it
    ax.text(0.10, dome_lp * 0.60,
            "0.10 beru\n(canonical)",
            fontsize=8.5, fontweight="bold", ha="center", va="top",
            color="#444444",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="#cccccc", alpha=0.85),
            zorder=6)

    # p-value annotations with arrows coming from below the data points
    ax.annotate(f"p = {dome_peak_p:.4f}",
                xy=(0.10, dome_lp),
                xytext=(0.128, dome_lp - 0.35),
                fontsize=8, color=C_DOME,
                arrowprops=dict(arrowstyle="->", color=C_DOME, lw=0.9),
                zorder=5)
    ax.annotate(f"p = {full_peak_p:.4f}",
                xy=(0.10, full_lp),
                xytext=(0.128, full_lp + 0.30),
                fontsize=8, color=C_FULL,
                arrowprops=dict(arrowstyle="->", color=C_FULL, lw=0.9),
                zorder=5)

    ax.set_xlabel("Harmonic spacing (beru)")
    ax.set_ylabel("$-\\log_{10}(p)$")
    ax.set_title("Unit sensitivity — signal is confined to the canonical 0.1-beru spacing")
    ax.legend(loc="upper left", framealpha=0.9)
    ax.set_xlim(0.04, 0.21)
    # Explicit y ceiling: 15 % above dome peak so title stays clear
    ax.set_ylim(0, y_max_data * 1.15)
    ax.set_xticks(coarse)
    ax.set_xticklabels([f"{s:.2f}" for s in coarse], fontsize=8, rotation=45)

    fig.tight_layout()
    outpath = OUTDIR / "fig_unitsweep.pdf"
    fig.savefig(outpath)
    fig.savefig(outpath.with_suffix(".png"))
    plt.close(fig)
    print(f"  ✓ {outpath.name} ({outpath.with_suffix('.png').name})")
    return outpath


# ══════════════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 70)
    print("  Generating manuscript figures")
    print("=" * 70)
    print()

    print("Figure 1 — Beru deviation histogram (fig:devhist)")
    p1 = make_devhist()

    print("\nFigure 2 — Temporal gradient (fig:temporal)")
    p2 = make_temporal()

    print("\nFigure 3 — Unit sensitivity sweep (fig:unitsweep)")
    p3 = make_unitsweep()

    print()
    print("=" * 70)
    print(f"  All figures saved to: {OUTDIR}/")
    print(f"  PDF + PNG for each figure.")
    print("=" * 70)
    print()
    print("To include in LaTeX manuscript:")
    print(f"  \\includegraphics[width=0.9\\textwidth]{{figures/fig_devhist}}")
    print(f"  \\includegraphics[width=0.9\\textwidth]{{figures/fig_temporal}}")
    print(f"  \\includegraphics[width=0.9\\textwidth]{{figures/fig_unitsweep}}")
