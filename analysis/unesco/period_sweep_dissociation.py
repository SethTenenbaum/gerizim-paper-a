"""
period_sweep_dissociation.py
============================
Full period sweep of the UNESCO keyword-matched dome/stupa corpus:
tests both Rayleigh concentration and binomial proximity enrichment
at each candidate spacing from 1.5° to 6.0° in 0.1° steps (45 periods).

This script was written to make the sweep underlying Table 1 of the paper
fully reproducible. The four periods shown in Table 1 (2.4°, 3.0°, 3.5°, 3.6°)
were selected from this sweep to represent the full range of outcomes:
  - Rayleigh-significant but no proximity enrichment (2.4°, 3.5°, 3.6°)
  - Proximity-enriched but Rayleigh-non-significant (3.0°)

The Tier-A threshold (±0.10 × T, geometric null p0 = 0.20) is used
throughout, matching Table 1.

Run from repo root:
    python3 analysis/unesco/period_sweep_dissociation.py

Outputs:
    - Console table of all 45 periods
    - Highlights the four Table 1 periods
    - Confirms 3.0° is the dominant proximity peak
"""

import sys
import json
import numpy as np
from pathlib import Path
from scipy.stats import binomtest

_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(_ROOT))

from data.unesco_corpus import load_corpus
from lib.dome_filter import is_dome_site_raw as is_dome_site
from lib.beru import GERIZIM
from lib.results_store import ResultsStore

_CFG   = json.loads((_ROOT / "config.json").read_text())
N_PERM = _CFG["simulation"]["n_permutations"]
SEED   = _CFG["simulation"]["random_seed"]

ANCHOR      = GERIZIM          # Mount Gerizim longitude (degrees)
TIER_A_FRAC = 0.10             # Tier-A threshold as fraction of T
PERIODS     = np.round(np.arange(1.5, 6.05, 0.1), 2)   # 1.5, 1.6, ..., 6.0
ALPHA       = 0.05
TABLE1_PERIODS = {2.4, 3.0, 3.5, 3.6}

# ── Load dome corpus ──────────────────────────────────────────────────────────
corpus   = load_corpus()
cultural = [s for s in corpus if s.category != "Natural" and s.has_coords]
dome     = [s for s in cultural if is_dome_site(s)]
lons     = np.array([s.longitude for s in dome])
N        = len(lons)

print()
print("=" * 80)
print("  PERIOD SWEEP DISSOCIATION: Rayleigh vs. Binomial Proximity")
print(f"  UNESCO dome/stupa corpus  N = {N}  |  anchor = {ANCHOR}°E")
print(f"  Periods: {PERIODS[0]}° to {PERIODS[-1]}° in 0.1° steps  ({len(PERIODS)} total)")
print(f"  Tier-A threshold: ±{TIER_A_FRAC} × T   (geometric null p0 = {2*TIER_A_FRAC:.2f})")
print(f"  Rayleigh permutation draws: {N_PERM:,}")
print("=" * 80)

rng = np.random.default_rng(SEED + 7)

# ── Helper functions ──────────────────────────────────────────────────────────

def nearest_harmonic_dev(lon, T, anchor):
    """Deviation from nearest harmonic of period T anchored at anchor (degrees)."""
    arc = (lon - anchor) % T
    return min(arc, T - arc)

def rayleigh_perm_p(lons, T, anchor, n_perm, rng_):
    """Permutation p-value for Rayleigh R at period T."""
    phases = 2 * np.pi * ((lons - anchor) % T) / T
    C, S = np.cos(phases).mean(), np.sin(phases).mean()
    R_obs = np.sqrt(C**2 + S**2)
    count = 0
    for _ in range(n_perm):
        ph = rng_.uniform(0, 2 * np.pi, len(lons))
        if np.sqrt(np.cos(ph).mean()**2 + np.sin(ph).mean()**2) >= R_obs:
            count += 1
    return R_obs, (count + 1) / (n_perm + 1)

def sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "** "
    if p < 0.05:  return "*  "
    if p < 0.10:  return "~  "
    return "ns "

# ── Sweep ─────────────────────────────────────────────────────────────────────

print()
print(f"  {'T':>5}  {'R':>6}  {'Ray-p':>7}  {'Ray-sig':>7}  "
      f"{'Hits':>5}  {'Exp':>5}  {'Bin-p':>7}  {'Bin-sig':>7}")
print("  " + "-" * 66)

results = []
for T in PERIODS:
    d = TIER_A_FRAC * T
    p_null = 2 * TIER_A_FRAC

    # Binomial
    hits = sum(1 for lon in lons if nearest_harmonic_dev(lon, T, ANCHOR) <= d)
    exp  = N * p_null
    bin_p = binomtest(hits, N, p_null, alternative="greater").pvalue

    # Rayleigh (permutation)
    R_obs, ray_p = rayleigh_perm_p(lons, T, ANCHOR, N_PERM, rng)

    marker = " <-- Table 1" if T in TABLE1_PERIODS else ""
    print(f"  {T:>5.1f}  {R_obs:>6.4f}  {ray_p:>7.4f}  {sig_label(ray_p):>7}  "
          f"{hits:>5}  {exp:>5.1f}  {bin_p:>7.4f}  {sig_label(bin_p):>7}{marker}")

    results.append(dict(T=T, R=R_obs, ray_p=ray_p, hits=hits, exp=exp, bin_p=bin_p))

# ── Summary: find dominant proximity peak and adjacent non-sig spacings ───────
sig_binom = [(r["T"], r["bin_p"]) for r in results if r["bin_p"] < ALPHA]
min_binom = min(results, key=lambda r: r["bin_p"])

print()
print("=" * 80)
print(f"  Dominant proximity peak:  T = {min_binom['T']}°   "
      f"binom-p = {min_binom['bin_p']:.4f}  {sig_label(min_binom['bin_p'])}")
print(f"  All significant proximity periods ({ALPHA}): "
      + ", ".join(f"{t}° (p={p:.4f})" for t, p in sig_binom))

# Check adjacent non-harmonic spacings around 3.0°
check_adj = [1.8, 2.1, 2.7, 3.3]
print()
print("  Adjacent non-harmonic spacings:")
for T_adj in check_adj:
    T_adj = round(T_adj, 1)
    r = next((x for x in results if abs(x["T"] - T_adj) < 0.05), None)
    if r:
        print(f"    T = {r['T']}°:  binom-p = {r['bin_p']:.4f}  {sig_label(r['bin_p'])}")

print("=" * 80)
print()
