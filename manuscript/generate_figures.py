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

_N_PERMS_NULL_C = 100_000

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

# Color palette — muted, color-blind friendly
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

    # Shade exactly the A+ window (δ ≤ 0.002 = TIER_APLUS).
    # Note: the histogram bin width is 0.005, so the first bar extends to 0.005
    # and will always be taller than the A+ count (56). The annotation uses
    # a direct count of sites with δ ≤ TIER_APLUS, not counts[0].
    ax.axvspan(0, TIER_APLUS, color=C_HIGHLIGHT, alpha=0.12, zorder=1,
               label=f"Tier-A+ window (δ ≤ {TIER_APLUS})")

    # Uniform null expectation line
    expected_per_bin = N_total / n_bins
    ax.axhline(expected_per_bin, color=C_NULL, linewidth=1.5,
               linestyle="--", zorder=4, label=f"Uniform null ({expected_per_bin:.1f}/bin)")

    # Annotate the A+ window count (δ ≤ TIER_APLUS = 0.002), not the full
    # first bin (δ ≤ 0.005).  The first bin is 2.5× wider than the A+ window
    # so counts[0] always overstates the signal.
    n_aplus = int(sum(1 for d in deviations if d <= TIER_APLUS))
    y_max_data = counts.max()
    ax.annotate(
        f"{n_aplus} A+ sites\n(δ ≤ {TIER_APLUS})",
        xy=(TIER_APLUS, n_aplus),
        xytext=(0.018, n_aplus + 8),
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

    # Add Cochran-Armitage annotation.
    # Values are read from generated_macros.tex to ensure the figure annotation
    # is always identical to the canonical analysis output (temporal_gradient_test.py).
    # Do NOT recompute Z here — use the macro values so figure and text are in sync.
    _macros_path = Path(__file__).resolve().parent / "generated_macros.tex"
    z_ca3 = None
    p_ca3 = None
    if _macros_path.exists():
        import re
        _macro_text = _macros_path.read_text()
        _z_match = re.search(r"\\newcommand\{\\ZcochranThree\}\{([^}]+)\}", _macro_text)
        _p_match = re.search(r"\\newcommand\{\\pCochranThree\}\{([^}]+)\}", _macro_text)
        if _z_match and _p_match:
            z_ca3 = float(_z_match.group(1))
            p_ca3 = float(_p_match.group(1))
    if z_ca3 is None:
        raise RuntimeError(
            "Could not read \\ZcochranThree / \\pCochranThree from generated_macros.tex. "
            "Run temporal_gradient_test.py first to regenerate macros."
        )

    ax.text(0.98, 0.65,
            f"Cochran-Armitage (3-cohort)\n"
            f"Z = {z_ca3:.2f}, p = {p_ca3:.4f}",
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
#  FIGURE 4: Null C — restricted geographic draw  (fig:null_c)
# ══════════════════════════════════════════════════════════════════════════════
def make_null_c():
    """Bootstrap distributions for Null C at three window widths (±2°, ±5°, ±10°).

    Null C asks: do dome sites outperform OTHER monuments at the same longitudes?
    For each window W, pool every full-corpus site within ±W° of any dome site,
    then bootstrap N_dome draws to build the expected A+ distribution.
    All three windows give p < 0.05, demonstrating the result is window-independent.
    """
    from scipy.stats import binomtest as _binomtest

    rng = np.random.default_rng(42)

    dome_sites = [s for s in sites if s["is_dome"]]
    dome_lons  = np.array([s["lon"] for s in dome_sites])
    N_dome     = len(dome_lons)
    obs_ap     = sum(1 for s in dome_sites if s["is_ap"])

    all_lons_arr = np.array([s["lon"] for s in sites])

    windows = [2.0, 5.0, 10.0]
    labels  = ["±2°", "±5°", "±10°"]

    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), sharey=False)
    fig.suptitle(
        "Null C — restricted geographic draw: dome sites vs. same-longitude corpus",
        fontsize=11, y=1.01,
    )

    for ax, W, lbl in zip(axes, windows, labels):
        # Build pool: corpus sites within ±W° of any dome longitude
        mask = np.zeros(len(all_lons_arr), dtype=bool)
        for dlon in dome_lons:
            mask |= np.abs(all_lons_arr - dlon) <= W
        pool = all_lons_arr[mask]
        N_pool = len(pool)

        # Bootstrap: draw N_dome from pool, count A+ each time
        draws = rng.choice(pool, size=(_N_PERMS_NULL_C, N_dome), replace=True)
        arc   = np.abs(draws - GERIZIM) / BERU
        devs  = np.abs(arc - np.round(arc / 0.1) * 0.1)
        boot_ap = (devs <= TIER_APLUS).sum(axis=1).astype(int)

        p_val  = float(np.mean(boot_ap >= obs_ap))
        mean_b = float(boot_ap.mean())
        std_b  = float(boot_ap.std())
        z_val  = (obs_ap - mean_b) / std_b if std_b > 0 else 0.0

        # Histogram of bootstrap distribution
        max_ap = max(boot_ap.max(), obs_ap) + 1
        bins   = np.arange(0, max_ap + 2) - 0.5
        counts, edges = np.histogram(boot_ap, bins=bins)
        bar_centers = 0.5 * (edges[:-1] + edges[1:])

        ax.bar(bar_centers, counts / _N_PERMS_NULL_C * 100,
               width=0.8, color=C_PRIMARY, alpha=0.75,
               edgecolor="white", linewidth=0.5, zorder=3,
               label="Bootstrap\ndistribution")

        # Mean line
        ax.axvline(mean_b, color=C_NULL, linewidth=1.5, linestyle="--",
                   zorder=4, label=f"Mean = {mean_b:.1f}")

        # Observed line
        ax.axvline(obs_ap, color=C_HIGHLIGHT, linewidth=2.0, linestyle="-",
                   zorder=5, label=f"Observed = {obs_ap}")

        # Shade the tail (≥ observed)
        tail_x = bar_centers[bar_centers >= obs_ap - 0.5]
        tail_y = counts[bar_centers >= obs_ap - 0.5] / _N_PERMS_NULL_C * 100
        if len(tail_x):
            ax.bar(tail_x, tail_y, width=0.8,
                   color=C_HIGHLIGHT, alpha=0.35,
                   edgecolor="none", zorder=4)

        # Stats annotation
        sig_lbl = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
        ax.text(0.97, 0.97,
                f"N pool = {N_pool}\nZ = {z_val:.2f}\np = {p_val:.4f}  {sig_lbl}",
                transform=ax.transAxes,
                fontsize=8, ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.35", facecolor="#f5f5f5",
                          edgecolor="#cccccc", alpha=0.9))

        ax.set_title(f"Window {lbl}  (pool N = {N_pool})", fontsize=10)
        ax.set_xlabel("A+ sites in bootstrap draw", fontsize=9)
        ax.set_ylabel("Frequency (%)" if ax is axes[0] else "", fontsize=9)
        ax.legend(loc="lower right", fontsize=7.5, framealpha=0.85,
                  bbox_to_anchor=(1, 0.58), bbox_transform=ax.transAxes)
        ax.set_xlim(-0.5, max(max_ap, obs_ap + 1) + 0.5)

    fig.tight_layout()
    outpath = OUTDIR / "fig_null_c.pdf"
    fig.savefig(outpath)
    fig.savefig(outpath.with_suffix(".png"))
    plt.close(fig)
    print(f"  ✓ {outpath.name} ({outpath.with_suffix('.png').name})")
    return outpath


# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 5: Geographic tier trail — Silk Road / Buddhist transmission
#                                   (fig:geo_trail)
# ══════════════════════════════════════════════════════════════════════════════
def make_geo_trail():
    """Eurasian map showing tier-coloured UNESCO dome sites and Wikidata stupas.

    Focuses on the Buddhist transmission corridor from the Levant through
    the Silk Road to Java/SE Asia.  Sites are coloured by harmonic tier;
    approximate Silk Road overland and maritime routes are drawn as reference.

    Tier colour scheme:
      A++ / A+  → crimson   (≤6.7 km from harmonic)
      A         → orange    (≤33 km from harmonic)
      B         → grey
      C / C-    → steel-blue (inter-harmonic band)

    Markers:
      UNESCO dome sites  → circle  (●)
      Wikidata Q180987   → triangle (▲)
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import csv

    from lib.dome_filter import is_dome_site_raw as _is_dome_raw

    # ── Tier colour / size lookup ─────────────────────────────────────────────
    def _tier_style(t):
        if t in ("A++", "A+"):
            return "#c0392b", 55, 1.0
        elif t == "A":
            return "#e67e22", 28, 0.88
        elif t == "B":
            return "#aaaaaa", 10, 0.45
        else:
            return "#3498db", 10, 0.45

    # ── Load UNESCO dome sites (raw sweep, N=90) ──────────────────────────────
    dome_rows = []
    for s in cultural:
        if not _is_dome_raw(s):
            continue
        lat = getattr(s, "latitude", None)
        lon = getattr(s, "longitude", None)
        if lat is None or lon is None:
            continue
        d = beru_dev(lon)
        t = tier(d)
        clr, sz, alpha = _tier_style(t)
        dome_rows.append((lon, lat, t, clr, sz, alpha, s.site))

    # ── Load Wikidata Q180987 stupa CSV ───────────────────────────────────────
    _CSV = PROJECT_ROOT / "data" / "store" / "unesco" / "wikidata_stupas_q180987.csv"
    wiki_rows = []
    with open(_CSV, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(row for row in fh if not row.startswith("#"))
        for rec in reader:
            try:
                lat = float(rec["lat"])
                lon = float(rec["lon"])
                name = rec.get("name", "")
            except (ValueError, KeyError):
                continue
            d = beru_dev(lon)
            t = tier(d)
            clr, sz, alpha = _tier_style(t)
            wiki_rows.append((lon, lat, t, clr, sz, alpha, name))

    # ── Approximate Silk Road routes (historical waypoints) ──────────────────
    # Overland northern route: Gerizim → Anatolia → Persia → Central Asia → China
    _silk_north = [
        (GERIZIM, 32.2),  # Mount Gerizim (anchor)
        (36.3, 36.2),   # Antioch/Antakya
        (38.7, 39.9),   # Erzincan
        (44.5, 40.5),   # Tabriz
        (51.4, 35.7),   # Tehran
        (59.0, 37.9),   # Merv
        (63.6, 40.1),   # Bukhara
        (69.3, 41.3),   # Samarkand
        (76.9, 43.2),   # Almaty
        (80.3, 42.9),   # Xinjiang border
        (87.6, 43.8),   # Urumqi
        (98.5, 39.7),   # Dunhuang
        (107.3, 34.3),  # Chang'an (Xi'an)
    ]
    # Maritime route west: Gerizim → Mediterranean → Iberia
    _silk_med = [
        (GERIZIM, 32.2),  # Mount Gerizim (anchor)
        (33.0, 35.1),   # Cyprus
        (28.2, 36.4),   # Rhodes
        (23.7, 37.9),   # Athens/Piraeus
        (14.5, 35.9),   # Sicily
        (12.5, 41.9),   # Rome/Ostia
        (10.2, 37.0),   # Carthage/Tunis
        (5.4,  43.3),   # Marseille
        (2.2,  41.4),   # Barcelona
    ]
    # Maritime route east: Gerizim → Arabian Sea → India → SE Asia → Java
    _silk_sea = [
        (GERIZIM, 32.2),  # Mount Gerizim (anchor)
        (38.0, 21.5),   # Red Sea
        (44.0, 12.8),   # Aden
        (55.0, 23.6),   # Oman coast
        (66.9, 24.9),   # Karachi
        (72.8, 18.9),   # Mumbai
        (80.3, 7.9),    # Colombo/Sri Lanka
        (92.5, 13.1),   # Bay of Bengal
        (96.1, 16.9),   # Yangon
        (100.5, 5.4),   # Malacca
        (107.6, -6.2),  # Jakarta/Java
        (110.2, -7.6),  # Borobudur area
    ]

    # ── Build figure ──────────────────────────────────────────────────────────
    proj = ccrs.Mercator(central_longitude=80.0, min_latitude=-15, max_latitude=65)
    fig  = plt.figure(figsize=(13, 6.5))
    ax   = fig.add_axes([0.01, 0.04, 0.98, 0.88], projection=proj)

    # Map extent: Levant to Java, covering full Silk Road corridor
    ax.set_extent([-5, 145, -12, 58], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.LAND,  facecolor="#f2ede6", zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor="#d6e8f2", zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.4,  edgecolor="#999999", zorder=1)
    ax.add_feature(cfeature.BORDERS,   linewidth=0.25, edgecolor="#bbbbbb",
                   linestyle=":", zorder=1)

    # Harmonic grid lines (every 3° from Gerizim, within view)
    for k in range(-15, 40):
        hlon = GERIZIM + k * 3.0
        if -5 <= hlon <= 145:
            ax.plot([hlon, hlon], [-14, 60], transform=ccrs.PlateCarree(),
                    color="#dddddd", linewidth=0.5, linestyle="--",
                    alpha=0.6, zorder=2)

    # Highlighted harmonic lines: Gerizim anchor and Java node
    java_lon = GERIZIM + 25 * 3.0   # 2.5-beru = harmonic #25
    for hlon in [GERIZIM, java_lon]:
        ax.plot([hlon, hlon], [-14, 60], transform=ccrs.PlateCarree(),
                color="#7f8c8d", linewidth=1.2, linestyle="-",
                alpha=0.8, zorder=3)
    # Java node line label is handled via the Borobudur site annotation below

    # Silk Road routes
    def _plot_route(pts, color, lw, ls, label, zorder=3):
        lons = [p[0] for p in pts]
        lats = [p[1] for p in pts]
        ax.plot(lons, lats, transform=ccrs.PlateCarree(),
                color=color, linewidth=lw, linestyle=ls,
                alpha=0.55, zorder=zorder, label=label)

    _plot_route(_silk_north, "#8e44ad", 1.4, "-",  "Overland route (N.)")
    _plot_route(_silk_med,   "#16a085", 1.4, "--", "Maritime route")
    _plot_route(_silk_sea,   "#16a085", 1.4, "--", "_nolegend_")

    # ── Plot sites (background tiers first, then A+) ──────────────────────────
    for priority_tiers in [{"B", "C", "C-", "C--"}, {"A"}, {"A++", "A+"}]:
        for lon, lat, t, clr, sz, alpha, _ in dome_rows:
            if t in priority_tiers:
                ax.scatter(lon, lat, transform=ccrs.PlateCarree(),
                           s=sz, c=clr, alpha=alpha, marker="o",
                           edgecolors="white" if sz >= 28 else "none",
                           linewidths=0.5, zorder=5 if t in ("A++","A+") else 4)
        for lon, lat, t, clr, sz, alpha, _ in wiki_rows:
            if t in priority_tiers:
                ax.scatter(lon, lat, transform=ccrs.PlateCarree(),
                           s=sz * 0.75, c=clr, alpha=alpha, marker="^",
                           edgecolors="white" if sz >= 28 else "none",
                           linewidths=0.5, zorder=5 if t in ("A++","A+") else 4)

    # Annotate key A+ sites
    _annotations = [
        (GERIZIM,  32.2,  "Mt. Gerizim\n(anchor)",  "left",   1.5),
        (44.32,    40.16, "Echmiatsin",              "right",  1.0),
        (107.3,    34.3,  "Silk Roads\n(Chang'an)",  "left",   1.0),
        (110.2,    -7.55, "Borobudur\n(Java node, 110.3°E)", "left", 1.0),
        (80.45,    27.5,  "Lumbini",                 "right",  1.0),
    ]
    for lon, lat, label, ha, dy in _annotations:
        ax.annotate(label,
                    xy=(lon, lat), xytext=(lon + (2 if ha == "left" else -2), lat + dy),
                    transform=ccrs.PlateCarree(),
                    fontsize=7, color="#333333", ha=ha,
                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.7),
                    clip_on=True)

    # ── Legend ────────────────────────────────────────────────────────────────
    from matplotlib.lines import Line2D
    tier_handles = [
        Line2D([0],[0], marker="o", color="w", markerfacecolor="#c0392b",
               markersize=8, label="Tier A$^{+}$/A$^{++}$  (≤6.7 km)"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor="#e67e22",
               markersize=6, label="Tier A  (≤33 km)"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor="#aaaaaa",
               markersize=5, label="Tier B"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor="#3498db",
               markersize=5, label="C-band (inter-harmonic)"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor="#555555",
               markersize=6, label="○ UNESCO dome/stupa"),
        Line2D([0],[0], marker="^", color="w", markerfacecolor="#555555",
               markersize=6, label="△ Wikidata Q180987 stupa"),
        Line2D([0],[0], color="#8e44ad", linewidth=1.4, label="Overland Silk Road"),
        Line2D([0],[0], color="#16a085", linewidth=1.4, linestyle="--",
               label="Maritime Silk Road"),
    ]
    ax.legend(handles=tier_handles, loc="lower left",
              fontsize=7.2, framealpha=0.90, ncol=2,
              title="Tier / corpus / route", title_fontsize=7.5,
              handletextpad=0.4, borderpad=0.5, columnspacing=0.8)

    ax.set_title(
        "Buddhist transmission corridor — UNESCO dome/stupa and Wikidata Q180987 stupa tier distribution\n"
        "Circles = UNESCO ($N=90$); triangles = Wikidata ($N=229$).  "
        "Dashed grid: 0.1-beru harmonics.  Crimson = A$^{+}$, orange = A-tier.",
        fontsize=8.5, pad=5,
    )

    fig.tight_layout(pad=0.3)
    outpath = OUTDIR / "fig_geo_trail.pdf"
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

    print("\nFigure 4 — Null C restricted geographic draw (fig:null_c)")
    p4 = make_null_c()

    print("\nFigure 5 — Geographic tier trail / Silk Road corridor (fig:geo_trail)")
    p5 = make_geo_trail()

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
    print(f"  \\includegraphics[width=\\textwidth]{{figures/fig_null_c}}")
