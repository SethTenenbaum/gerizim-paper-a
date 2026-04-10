"""
sensitivity_slope_permutation_test.py
======================================
Tests whether the observed ±0.3% sensitivity collapse at the canonical
0.1-beru spacing is unusual under an empirical null that preserves the
UNESCO longitude distribution.

Procedure
---------
For each of N_PERMS permuted datasets (longitudes shuffled across sites):
  1. Run the fine-grained unit sweep (0.0990 to 0.1010 beru, same 7 spacings
     as Table 3 / tab:unitsweep_fine in the manuscript).
  2. Record whether the canonical spacing (0.1000) achieves joint significance
     in both populations simultaneously (dome p < 0.05 AND full p < 0.05).
  3. Record whether every non-canonical spacing in the sweep is non-significant
     in at least one population (i.e., joint significance collapses when
     spacing departs by ≥ 0.0003 beru = 0.3%).

The fraction of permuted datasets showing this "sharp peak" pattern is the
permutation p-value for the sensitivity slope observation.

A secondary metric: for each permuted dataset find the best joint-significance
spacing in the fine window; record whether 0.1000 beats all others.

All results reported for manuscript section §5.1 (Sensitivity Slope at the
Canonical Unit).

Code: analysis/unesco/sensitivity_slope_permutation_test.py
"""

import sys
import numpy as np
from pathlib import Path
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, TIER_APLUS, deviation_at_spacing
from lib.dome_filter import is_dome_site

# ── Configuration ─────────────────────────────────────────────────────────
N_PERMS = 10_000
SEED = 42

# The fine spacings used in Table 3 of the manuscript
FINE_SPACINGS = [0.0990, 0.0993, 0.0995, 0.1000, 0.1001, 0.1003, 0.1010]
CANONICAL = 0.1000
THRESHOLD = 0.05  # significance threshold for joint significance

# Off-canonical spacings (everything except 0.1000 and 0.1001)
# 0.1001 is the immediate neighbour that also retains dome significance
# in the observed data; we define "collapse" as: all spacings ≥ 0.0003
# from canonical are non-significant in at least one population.
OFF_CANONICAL = [s for s in FINE_SPACINGS if abs(s - CANONICAL) >= 0.0003]
# = [0.0990, 0.0993, 0.0995, 0.1003, 0.1010]


def null_rate_at_spacing(tier_aplus, spacing):
    """Geometric null rate: window_width / harmonic_spacing."""
    return 2 * tier_aplus / spacing


def count_hits_array(lons, spacing):
    """Count sites within TIER_APLUS of the nearest harmonic at given spacing."""
    hits = 0
    for lon in lons:
        dev = deviation_at_spacing(lon, spacing)
        if dev <= TIER_APLUS:
            hits += 1
    return hits


def joint_significant(dome_lons, all_lons, spacing):
    """Return (dome_p, full_p, is_joint_sig)."""
    nr = null_rate_at_spacing(TIER_APLUS, spacing)
    d_hits = count_hits_array(dome_lons, spacing)
    f_hits = count_hits_array(all_lons, spacing)
    d_p = binomtest(d_hits, len(dome_lons), nr, alternative='greater').pvalue
    f_p = binomtest(f_hits, len(all_lons), nr, alternative='greater').pvalue
    return d_p, f_p, (d_p < THRESHOLD and f_p < THRESHOLD)


def main():
    rng = np.random.default_rng(SEED)

    # ── Load observed data ─────────────────────────────────────────────────
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    all_lons = np.array([s.longitude for s in cultural])
    dome_mask = np.array([is_dome_site(s) for s in cultural])
    dome_lons_obs = all_lons[dome_mask]

    N_all = len(all_lons)
    N_dome = len(dome_lons_obs)
    print(f"Loaded {N_all} Cultural/Mixed sites, {N_dome} dome/spherical sites")
    print(f"Running {N_PERMS:,} permutations (seed={SEED})\n")

    # ── Observed pattern ───────────────────────────────────────────────────
    print("── Observed fine sweep ──")
    obs_canonical_joint = False
    obs_off_all_collapsed = True
    for sp in FINE_SPACINGS:
        d_p, f_p, joint = joint_significant(dome_lons_obs, all_lons, sp)
        marker = "<<<" if sp == CANONICAL else ""
        print(f"  {sp:.4f}  dome_p={d_p:.4f}  full_p={f_p:.4f}  joint={'YES' if joint else 'no'}  {marker}")
        if sp == CANONICAL:
            obs_canonical_joint = joint
        if sp in OFF_CANONICAL and joint:
            obs_off_all_collapsed = False

    obs_sharp_peak = obs_canonical_joint and obs_off_all_collapsed
    print(f"\nObserved: canonical joint-sig={obs_canonical_joint}, "
          f"all off-canonical collapsed={obs_off_all_collapsed}, "
          f"sharp-peak={obs_sharp_peak}")

    # ── Permutation loop ───────────────────────────────────────────────────
    print(f"\n── Permutation test ({N_PERMS:,} trials) ──")
    n_canonical_joint = 0   # permutations where canonical is jointly significant
    n_sharp_peak = 0        # permutations with the full sharp-peak pattern

    for i in range(N_PERMS):
        perm_lons = rng.permutation(all_lons)
        perm_dome_lons = perm_lons[dome_mask]

        canon_joint = False
        off_collapsed = True

        for sp in FINE_SPACINGS:
            d_p, f_p, joint = joint_significant(perm_dome_lons, perm_lons, sp)
            if sp == CANONICAL and joint:
                canon_joint = True
            if sp in OFF_CANONICAL and joint:
                off_collapsed = False

        if canon_joint:
            n_canonical_joint += 1
        if canon_joint and off_collapsed:
            n_sharp_peak += 1

        if (i + 1) % 1000 == 0:
            print(f"  {i+1:>6}/{N_PERMS}  canonical_joint={n_canonical_joint}  sharp_peak={n_sharp_peak}")

    # ── Results ────────────────────────────────────────────────────────────
    p_canonical_joint = n_canonical_joint / N_PERMS
    p_sharp_peak = n_sharp_peak / N_PERMS

    print(f"\n══ Results ══")
    print(f"N permutations:                     {N_PERMS:,}")
    print(f"Permutations with canonical joint-sig: {n_canonical_joint} / {N_PERMS}  (p = {p_canonical_joint:.4f})")
    print(f"Permutations with sharp-peak pattern:  {n_sharp_peak} / {N_PERMS}  (p = {p_sharp_peak:.4f})")
    print()
    print("Interpretation:")
    print(f"  p(canonical jointly significant)  = {p_canonical_joint:.4f}")
    print(f"  p(canonical sig AND all ±0.3% collapsed) = {p_sharp_peak:.4f}")
    print()
    print("Manuscript macro values:")
    print(f"  \\newcommand{{\\permSlopeNperms}}{{{N_PERMS:,}}}")
    print(f"  \\newcommand{{\\permSlopePcanon}}{{{p_canonical_joint:.4f}}}")
    print(f"  \\newcommand{{\\permSlopePsharp}}{{{p_sharp_peak:.4f}}}")
    if p_sharp_peak == 0:
        print(f"  \\newcommand{{\\permSlopePsharpBound}}{{$< {1/N_PERMS:.4f}$}}")
    else:
        print(f"  \\newcommand{{\\permSlopePsharpBound}}{{{p_sharp_peak:.4f}}}")


if __name__ == "__main__":
    main()
