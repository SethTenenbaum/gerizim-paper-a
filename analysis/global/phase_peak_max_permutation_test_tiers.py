"""
x18_max_permutation_test_tiers.py
GROUP: 11d

Discovery-corrected max-permutation test run for all three A tiers:
A++ (TIER_APP), A+ (TIER_APLUS), A (TIER_A_MAX).

Same null model as x18_max_permutation_test.py (GROUP 11c):
  - Observed: maximum tier-X count across all 36,000 trial anchors.
  - Null A:   N phases drawn uniformly from [0, 3°); record max tier count.
  - Null B:   bootstrap resample from observed phases (conservative).
"""

import sys
import time
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import BERU, TIER_APP, TIER_APLUS, TIER_A_MAX
from lib.results_store import ResultsStore

N_PERMS      = 10_000
SEED         = 42
ANCHOR_STEP  = 0.01
HARMONIC_STEP = 3.0

GERIZIM_ID = "gerizim_synthetic"
corpus   = load_corpus()
from data.unesco_corpus import cultural_sites_with_coords_extended
cultural = cultural_sites_with_coords_extended(corpus)
lons     = np.array([s.longitude for s in cultural if s.id_number != GERIZIM_ID])
N        = len(lons)

phases_obs   = lons % HARMONIC_STEP
anchor_phases = np.arange(0.0, HARMONIC_STEP, ANCHOR_STEP)

def max_count_fast(site_phases: np.ndarray, thresh_deg: float) -> int:
    diff   = np.abs(site_phases[:, None] - anchor_phases[None, :])
    diff   = np.minimum(diff, HARMONIC_STEP - diff)
    counts = np.sum(diff <= thresh_deg, axis=0)
    return int(counts.max())

def _sig(p):
    return ("***" if p < 0.001 else "**" if p < 0.01 else
            "*"   if p < 0.05  else "†"  if p < 0.10 else "ns")

tiers = [
    ("A++", TIER_APP   * BERU),
    ("A+",  TIER_APLUS * BERU),
    ("A",   TIER_A_MAX * BERU),
]

rng = np.random.default_rng(SEED)

print("=" * 70)
print("  x.18° BAND — MAX-PERMUTATION TEST — ALL A TIERS")
print(f"  N={N}  {N_PERMS:,} permutations (seed={SEED})")
print("=" * 70)

rs = ResultsStore()
store_out = {}

for tier_label, thresh_deg in tiers:
    obs_max = max_count_fast(phases_obs, thresh_deg)

    perm_max_A = np.empty(N_PERMS, dtype=int)
    perm_max_B = np.empty(N_PERMS, dtype=int)

    t0 = time.time()
    rng2 = np.random.default_rng(SEED)   # reset seed per tier for reproducibility
    for i in range(N_PERMS):
        rand_phases    = rng2.uniform(0.0, HARMONIC_STEP, size=N)
        perm_max_A[i]  = max_count_fast(rand_phases, thresh_deg)
        boot_phases    = rng2.choice(phases_obs, size=N, replace=True)
        perm_max_B[i]  = max_count_fast(boot_phases, thresh_deg)

    p_A  = float((np.sum(perm_max_A >= obs_max) + 1) / (N_PERMS + 1))
    mu_A = float(perm_max_A.mean())
    sd_A = float(perm_max_A.std(ddof=1))
    z_A  = float((obs_max - mu_A) / sd_A) if sd_A > 0 else float("nan")

    p_B  = float((np.sum(perm_max_B >= obs_max) + 1) / (N_PERMS + 1))
    mu_B = float(perm_max_B.mean())
    sd_B = float(perm_max_B.std(ddof=1))

    elapsed = time.time() - t0
    print(f"\n  Tier {tier_label}  (thresh = {thresh_deg:.4f}°)")
    print(f"    Observed max count          : {obs_max}")
    print(f"    Null A mean / SD            : {mu_A:.2f} / {sd_A:.2f}")
    print(f"    Z (Null A)                  : {z_A:.2f}")
    print(f"    p (Null A, discovery-corr.) : {p_A:.4f}  {_sig(p_A)}")
    print(f"    p (Null B, bootstrap)       : {p_B:.4f}  {_sig(p_B)}")
    print(f"    ({elapsed:.1f}s)")

    label_clean = tier_label.replace("+", "p").replace("-", "m")
    store_out[f"anchorMaxPerm{label_clean}ObsMax"] = obs_max
    store_out[f"anchorMaxPerm{label_clean}NullMu"] = mu_A
    store_out[f"anchorMaxPerm{label_clean}NullSD"] = sd_A
    store_out[f"anchorMaxPerm{label_clean}Z"]      = z_A
    store_out[f"anchorMaxPerm{label_clean}P"]      = p_A
    store_out[f"anchorMaxPerm{label_clean}BootP"]  = p_B

print()
print("=" * 70)
rs.write_many(store_out)
print("Results written to data/store/results.json")
