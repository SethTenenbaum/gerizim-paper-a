"""
plot_period_sweep.py
====================
Generate Figure: Rayleigh vs. Binomial period sweep for the UNESCO dome corpus.

Sweeps candidate periods 1.5°–6.0° in 0.1° steps and computes:
  - Rayleigh permutation p (mean circular concentration)
  - Binomial proximity-enrichment p (Tier-A threshold: ±0.10 × T)

Outputs:
    manuscript/figures/fig_period_sweep.pdf
    manuscript/figures/fig_period_sweep.png

Run from repo root:
    python3 analysis/unesco/plot_period_sweep.py

The figure illustrates the Rayleigh–binomial dissociation: periods where
mean phase concentration (Rayleigh) is significant may show no proximity
enrichment (binomial), and vice versa.  T = 3.0° is the dominant binomial
proximity peak.

N_PERM defaults to 10,000 for speed (sufficient for stable figure curves).
Pass --full to use the config value (100,000) for publication-quality output.
"""

import sys
import json
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import binomtest

_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_ROOT))

from data.unesco_corpus import load_corpus
from lib.dome_filter import is_dome_site_raw as is_dome_site
from lib.beru import GERIZIM

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Generate period sweep figure")
parser.add_argument("--full", action="store_true",
                    help="Use config n_permutations (100k) instead of default 10k")
args = parser.parse_args()

_CFG   = json.loads((_ROOT / "config.json").read_text())
N_PERM = _CFG["simulation"]["n_permutations"] if args.full else 10_000
SEED   = _CFG["simulation"]["random_seed"]

ANCHOR      = GERIZIM
TIER_A_FRAC = 0.10
PERIODS     = np.round(np.arange(1.5, 6.05, 0.1), 2)
ALPHA       = 0.05
TABLE1_PERIODS = {2.4, 3.0, 3.5, 3.6}

# ── Load dome corpus ───────────────────────────────────────────────────────────
corpus   = load_corpus()
cultural = [s for s in corpus if s.category != "Natural" and s.has_coords]
dome     = [s for s in cultural if is_dome_site(s)]
lons     = np.array([s.longitude for s in dome])
N        = len(lons)

print(f"UNESCO dome corpus: N={N}  |  anchor={ANCHOR}°E  |  N_PERM={N_PERM:,}")

# ── Sweep ──────────────────────────────────────────────────────────────────────
rng = np.random.default_rng(SEED + 7)

def nearest_harmonic_dev(lon, T, anchor):
    arc = (lon - anchor) % T
    return min(arc, T - arc)

ray_p_vals  = []
bin_p_vals  = []

for T in PERIODS:
    d      = TIER_A_FRAC * T
    p_null = 2 * TIER_A_FRAC

    # Binomial proximity enrichment
    hits  = sum(1 for lon in lons if nearest_harmonic_dev(lon, T, ANCHOR) <= d)
    bin_p = binomtest(hits, N, p_null, alternative="greater").pvalue

    # Rayleigh (permutation-based)
    phases = 2 * np.pi * ((lons - ANCHOR) % T) / T
    C, S   = np.cos(phases).mean(), np.sin(phases).mean()
    R_obs  = np.sqrt(C**2 + S**2)
    null_R = np.array([
        np.sqrt(np.cos(ph).mean()**2 + np.sin(ph).mean()**2)
        for ph in (rng.uniform(0, 2 * np.pi, (N_PERM, N)))
    ])
    ray_p = (np.sum(null_R >= R_obs) + 1) / (N_PERM + 1)

    ray_p_vals.append(ray_p)
    bin_p_vals.append(bin_p)
    print(f"  T={T:.1f}°  ray_p={ray_p:.4f}  bin_p={bin_p:.4f}")

ray_p_vals = np.array(ray_p_vals)
bin_p_vals = np.array(bin_p_vals)

# ── Figure ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":      "serif",
    "font.size":        10,
    "axes.labelsize":   11,
    "axes.titlesize":   11,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "legend.fontsize":  9,
    "figure.dpi":       300,
    "axes.spines.top":  False,
    "axes.spines.right": False,
})

C_BINOM   = "#2c6fbb"   # deep blue — binomial
C_RAYLEIGH = "#e6550d"  # orange    — Rayleigh
C_THRESH  = "#888888"   # gray      — significance threshold
C_PEAK    = "#d62728"   # red       — 3° peak marker

fig, ax = plt.subplots(figsize=(7.5, 4.5))

neg_log_ray = -np.log10(np.clip(ray_p_vals, 1e-10, 1))
neg_log_bin = -np.log10(np.clip(bin_p_vals, 1e-10, 1))
thresh_line = -np.log10(ALPHA)   # ≈ 1.301

# Shade region of dissociation: binomial significant but Rayleigh not
dissoc_mask = (bin_p_vals < ALPHA) & (ray_p_vals >= ALPHA)
for i, T in enumerate(PERIODS):
    if dissoc_mask[i]:
        ax.axvspan(T - 0.05, T + 0.05, color=C_BINOM, alpha=0.08, linewidth=0)

ax.plot(PERIODS, neg_log_ray, color=C_RAYLEIGH, linewidth=1.6,
        label="Rayleigh (circular concentration)", zorder=3)
ax.plot(PERIODS, neg_log_bin, color=C_BINOM, linewidth=1.6,
        label=r"Binomial (proximity enrichment, $\pm0.10\,T$)", zorder=4)

# Significance threshold
ax.axhline(thresh_line, color=C_THRESH, linewidth=0.9, linestyle="--",
           label=fr"$\alpha$ = 0.05")

# 3° vertical guide
ax.axvline(3.0, color=C_PEAK, linewidth=0.9, linestyle=":", zorder=2)
y_top = max(neg_log_bin.max(), neg_log_ray.max()) * 1.02
ax.text(3.0 + 0.07, y_top, "$T = 3°$", color=C_PEAK, fontsize=8, va="top")

ax.set_xlabel("Candidate period $T$ (degrees)")
ax.set_ylabel(r"$-\log_{10}(p)$")
ax.set_xlim(PERIODS[0], PERIODS[-1])
ax.set_xticks(np.arange(1.5, 6.5, 0.5))
ax.set_xticklabels([f"{v:.1f}" for v in np.arange(1.5, 6.5, 0.5)])
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.10), ncol=3,
          frameon=False)

fig.tight_layout()

OUTDIR = _ROOT / "manuscript" / "figures"
OUTDIR.mkdir(exist_ok=True)

_SAVE_KW = dict(dpi=300, bbox_inches="tight", pad_inches=0.08)
out_pdf = OUTDIR / "fig_period_sweep.pdf"
out_png = OUTDIR / "fig_period_sweep.png"
fig.savefig(out_pdf, **_SAVE_KW)
fig.savefig(out_png, **_SAVE_KW)
plt.close(fig)

print(f"\nSaved: {out_pdf}")
print(f"Saved: {out_png}")

# Confirm dominant peak
peak_idx = np.argmin(bin_p_vals)
print(f"\nDominant binomial peak: T={PERIODS[peak_idx]:.1f}°  "
      f"bin_p={bin_p_vals[peak_idx]:.4f}  "
      f"ray_p={ray_p_vals[peak_idx]:.4f}")
