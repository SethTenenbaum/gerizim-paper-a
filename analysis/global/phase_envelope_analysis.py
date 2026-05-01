"""
phase_envelope_analysis.py
==========================
Phase-marginalized robustness analysis for the 3° periodicity.

For every possible phase origin φ ∈ [0, 3°) at 0.01° resolution, we recompute
the A+ tier count for (a) the full UNESCO Cultural/Mixed corpus and (b) the
dome/spherical sub-corpus. We then report:

  * min, median, max A+ count across all 300 possible phase origins
  * fraction of phase choices for which the count exceeds the geometric null
  * the phase that maximises each count (φ*)
  * Gerizim's phase mod 3° (for context)

This turns the anchor from a chosen parameter into a marginalized nuisance
parameter — the headline statistic is now the *minimum* enrichment across
all phase choices, which is anchor-free by construction.

USAGE
-----
  python3 analysis/global/phase_envelope_analysis.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import numpy as np
from scipy.stats import binomtest

from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, TIER_APLUS_DEG, P_NULL_AP
from lib.dome_filter import is_dome_site_raw
from lib.results_store import ResultsStore

HARMONIC_STEP_DEG = 3.0
PHASE_STEP        = 0.01           # resolution of phase scan in degrees
PHASE_GRID        = np.arange(0.0, HARMONIC_STEP_DEG, PHASE_STEP)
N_PHASES          = len(PHASE_GRID)   # 300


def ap_count(lons: np.ndarray, phase: float) -> int:
    """Count A+ sites when the harmonic grid is phase-shifted by `phase` degrees.

    A site is A+ iff its longitude minus the phase, modulo 3°, is within
    TIER_APLUS_DEG of either 0 or 3° (i.e. close to a grid node).
    """
    folded = (lons - phase) % HARMONIC_STEP_DEG
    dev = np.minimum(folded, HARMONIC_STEP_DEG - folded)
    return int(np.sum(dev <= TIER_APLUS_DEG))


def envelope(lons: np.ndarray) -> dict:
    counts = np.array([ap_count(lons, p) for p in PHASE_GRID])
    n = len(lons)
    null_count = P_NULL_AP * n
    # Per-phase one-sided binomial p (against geometric null)
    pvals = np.array([
        binomtest(int(c), n, P_NULL_AP, alternative="greater").pvalue
        for c in counts
    ])
    return {
        "n_sites":         n,
        "counts":          counts,
        "min":             int(counts.min()),
        "median":          float(np.median(counts)),
        "max":             int(counts.max()),
        "null_expected":   float(null_count),
        "phi_at_max":      float(PHASE_GRID[counts.argmax()]),
        "phi_at_min":      float(PHASE_GRID[counts.argmin()]),
        "frac_above_null": float(np.mean(counts > null_count)),
        "frac_sig_05":     float(np.mean(pvals < 0.05)),
        "frac_sig_001":    float(np.mean(pvals < 0.001)),
        "min_pvalue":      float(pvals.min()),
        "worst_pvalue":    float(pvals.max()),
        "median_pvalue":   float(np.median(pvals)),
    }


# ── Load corpora ─────────────────────────────────────────────────────────────
corpus   = load_corpus()
cultural = cultural_sites_with_coords(corpus)
all_lons = np.array([s.longitude for s in cultural])
dome_lons = np.array([s.longitude for s in cultural if is_dome_site_raw(s)])

# ── Compute envelopes ────────────────────────────────────────────────────────
full = envelope(all_lons)
dome = envelope(dome_lons)

# Gerizim's phase mod 3°
ger_phase = GERIZIM % HARMONIC_STEP_DEG
# A+ count at Gerizim phase is computed by treating Gerizim itself as the
# anchor: a site at longitude λ has deviation |λ - GERIZIM| mod 3°, so the
# equivalent grid origin is GERIZIM mod 3°.
ger_full = ap_count(all_lons,  ger_phase)
ger_dome = ap_count(dome_lons, ger_phase)

# Gerizim percentile across the phase grid (higher = better)
ger_full_pctile = float(100.0 * np.mean(full["counts"] <= ger_full))
ger_dome_pctile = float(100.0 * np.mean(dome["counts"] <= ger_dome))

print("=" * 74)
print("  PHASE-ENVELOPE ANALYSIS — A+ count across all φ ∈ [0, 3°)")
print("=" * 74)
print(f"  Phase grid: {N_PHASES} phases × {PHASE_STEP}° step")
print(f"  Tier-A+ threshold: ±{TIER_APLUS_DEG}°  (null rate {P_NULL_AP:.4f})")
print()
for name, env, ger_c in [("Full corpus", full, ger_full), ("Dome subset", dome, ger_dome)]:
    print(f"  {name}  (N = {env['n_sites']})")
    print(f"    Null expectation       : {env['null_expected']:6.2f}")
    print(f"    Min A+ count           : {env['min']:6d}   at φ = {env['phi_at_min']:.2f}°")
    print(f"    Median A+ count        : {env['median']:6.1f}")
    print(f"    Max  A+ count          : {env['max']:6d}   at φ = {env['phi_at_max']:.2f}°")
    print(f"    Gerizim phase A+ count : {ger_c:6d}   at φ = {ger_phase:.4f}°")
    print(f"    Frac. phases > null    : {env['frac_above_null']*100:5.1f}%")
    print(f"    Frac. phases p<0.05    : {env['frac_sig_05']*100:5.1f}%")
    print(f"    Frac. phases p<0.001   : {env['frac_sig_001']*100:5.1f}%")
    print(f"    Worst-case (max) p     : {env['worst_pvalue']:.4g}")
    print(f"    Median p across phases : {env['median_pvalue']:.4g}")
    print()

# ── LaTeX macros ─────────────────────────────────────────────────────────────
print("  % LaTeX macros (phase_envelope_analysis):")
macros = {
    "phaseFullMin":         full["min"],
    "phaseFullMedian":      f"{full['median']:.1f}",
    "phaseFullMax":         full["max"],
    "phaseFullPhiMax":      f"{full['phi_at_max']:.2f}",
    "phaseFullFracSig":     f"{full['frac_sig_05']*100:.1f}",
    "phaseFullFracSigStrict": f"{full['frac_sig_001']*100:.1f}",
    "phaseFullWorstP":      f"{full['worst_pvalue']:.4g}",
    "phaseFullMedianP":     f"{full['median_pvalue']:.4g}",
    "phaseFullGerCount":    ger_full,
    "phaseFullGerPctile":   f"{ger_full_pctile:.1f}",
    "phaseDomeMin":         dome["min"],
    "phaseDomeMedian":      f"{dome['median']:.1f}",
    "phaseDomeMax":         dome["max"],
    "phaseDomePhiMax":      f"{dome['phi_at_max']:.2f}",
    "phaseDomeFracSig":     f"{dome['frac_sig_05']*100:.1f}",
    "phaseDomeFracSigStrict": f"{dome['frac_sig_001']*100:.1f}",
    "phaseDomeWorstP":      f"{dome['worst_pvalue']:.4g}",
    "phaseDomeMedianP":     f"{dome['median_pvalue']:.4g}",
    "phaseDomeGerCount":    ger_dome,
    "phaseDomeGerPctile":   f"{ger_dome_pctile:.1f}",
    "phaseGerizimPhase":    f"{ger_phase:.4f}",
    "phaseNphases":         N_PHASES,
    "phaseStep":            f"{PHASE_STEP}",
}
for k, v in macros.items():
    print(f"  \\newcommand{{\\{k}}}{{{v}}}")

# ── Persist to results store ─────────────────────────────────────────────────
store = ResultsStore()
store.write_many(macros)

print()
print("  ✓ Phase envelope macros emitted and stored.")
