"""
Region-conditioned permutation test for the beru-harmonic A+ enrichment.

Two complementary within-group permutation strategies are run:

  Strategy A (UNESCO region x 20-degree longitude bin):
    Sites are grouped by the cross of UNESCO administrative region and a
    20-degree longitude bin.  Longitudes are shuffled within each cell.

  Strategy B (5-degree longitude blocks, global):
    Sites are grouped into 5-degree longitude bins globally.

Macros emitted (GROUP 35): regionCondPermP/Z/Obs/Mean/Std/Nregions/Nperms/BPermP/BPermZ
"""

import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
np.random.seed(42)

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, TIER_APLUS
from lib.results_store import ResultsStore

N_PERMS      = 100_000
HARMONIC_DEG = 3.0
APLUS_DEG    = TIER_APLUS * 30.0   # 0.002 * 30 = 0.06 deg


def count_aplus_vec(lons):
    shifted = (lons - GERIZIM) % HARMONIC_DEG
    dev_deg = np.minimum(shifted, HARMONIC_DEG - shifted)
    return int(np.sum(dev_deg <= APLUS_DEG))


def run_permutation(groups, obs_ap, n_perms, label):
    keys = list(groups.keys())
    perm_counts = np.zeros(n_perms, dtype=int)
    for trial in range(n_perms):
        trial_ap = 0
        for k in keys:
            arr = groups[k].copy()
            np.random.shuffle(arr)
            trial_ap += count_aplus_vec(arr)
        perm_counts[trial] = trial_ap
        if (trial + 1) % 10_000 == 0:
            print(f"    [{label}] {trial+1:,} / {n_perms:,}", file=sys.stderr)
    perm_mean = float(perm_counts.mean())
    perm_std  = float(perm_counts.std())
    p_val     = (int(np.sum(perm_counts >= obs_ap)) + 1) / (n_perms + 1)
    z_val     = (obs_ap - perm_mean) / perm_std if perm_std > 0 else float("nan")
    return dict(perm_mean=perm_mean, perm_std=perm_std,
                p_val=p_val, z_val=z_val, n_groups=len(groups))


def main():
    print("=" * 80)
    print("  REGION-CONDITIONED PERMUTATION TEST (two strategies)")
    print(f"  Anchor: {GERIZIM}E  |  A+ threshold: {TIER_APLUS} beru  |  seed: 42")
    print("=" * 80)

    corpus = load_corpus()
    sites  = cultural_sites_with_coords(corpus)

    lons          = np.array([s.longitude for s in sites])
    region_labels = [s.regions for s in sites]
    N             = len(lons)
    obs_ap        = count_aplus_vec(lons)
    print(f"\n  Total sites: {N}  |  Observed A+: {obs_ap}  ({obs_ap/N*100:.2f}%)")

    # Strategy A: UNESCO region x 20-degree longitude bin
    LON_BIN_A = 20
    groups_a = defaultdict(list)
    for lon, reg in zip(lons, region_labels):
        bin_id = int(np.floor(lon / LON_BIN_A)) * LON_BIN_A
        groups_a[f"{reg}||{bin_id}"].append(lon)
    groups_a_arr = {k: np.array(v) for k, v in groups_a.items()}

    print(f"\n  Strategy A: UNESCO region x {LON_BIN_A}-degree longitude bin")
    print(f"    Groups: {len(groups_a_arr)}")
    res_a = run_permutation(groups_a_arr, obs_ap, N_PERMS, "A")
    sig_a = "***" if res_a["p_val"] < 0.001 else "**" if res_a["p_val"] < 0.01 else "*" if res_a["p_val"] < 0.05 else "ns"
    print(f"    mean={res_a['perm_mean']:.2f}  std={res_a['perm_std']:.2f}  Z={res_a['z_val']:.3f}  p={res_a['p_val']:.4f}  {sig_a}")

    # Strategy B: 5-degree global longitude blocks
    LON_BIN_B = 5
    groups_b = defaultdict(list)
    for lon in lons:
        bin_id = int(np.floor(lon / LON_BIN_B)) * LON_BIN_B
        groups_b[bin_id].append(lon)
    groups_b_arr = {k: np.array(v) for k, v in groups_b.items()}

    print(f"\n  Strategy B: {LON_BIN_B}-degree global longitude blocks")
    print(f"    Blocks: {len(groups_b_arr)}")
    res_b = run_permutation(groups_b_arr, obs_ap, N_PERMS, "B")
    sig_b = "***" if res_b["p_val"] < 0.001 else "**" if res_b["p_val"] < 0.01 else "*" if res_b["p_val"] < 0.05 else "ns"
    print(f"    mean={res_b['perm_mean']:.2f}  std={res_b['perm_std']:.2f}  Z={res_b['z_val']:.3f}  p={res_b['p_val']:.4f}  {sig_b}")

    def fmt_p(p):
        return "< 0.001" if p < 0.001 else f"{p:.4f}"
    def fmt_z(z):
        return "---" if (z != z) else f"{z:.2f}"

    print("\n  % LaTeX macros (GROUP 35):")
    print(f"  \\newcommand{{\\regionCondPermP}}{{{fmt_p(res_a['p_val'])}}}     % Strategy A p")
    print(f"  \\newcommand{{\\regionCondPermZ}}{{{fmt_z(res_a['z_val'])}}}      % Strategy A Z")
    print(f"  \\newcommand{{\\regionCondPermObs}}{{{obs_ap}}}              % observed A+")
    print(f"  \\newcommand{{\\regionCondPermMean}}{{{res_a['perm_mean']:.1f}}}     % Strategy A perm mean")
    print(f"  \\newcommand{{\\regionCondPermStd}}{{{res_a['perm_std']:.1f}}}      % Strategy A perm std")
    print(f"  \\newcommand{{\\regionCondNregions}}{{{res_a['n_groups']}}}          % Strategy A groups")
    print(f"  \\newcommand{{\\regionCondNperms}}{{{N_PERMS:,}}}        % permutations")
    print(f"  \\newcommand{{\\regionCondBPermP}}{{{fmt_p(res_b['p_val'])}}}     % Strategy B p")
    print(f"  \\newcommand{{\\regionCondBPermZ}}{{{fmt_z(res_b['z_val'])}}}      % Strategy B Z")
    print(f"  \\newcommand{{\\regionCondBNblocks}}{{{res_b['n_groups']}}}         % Strategy B blocks")

    ResultsStore().write_many({
        "regionCondPermP":    res_a["p_val"],
        "regionCondPermZ":    res_a["z_val"],
        "regionCondPermObs":  obs_ap,
        "regionCondPermMean": res_a["perm_mean"],
        "regionCondPermStd":  res_a["perm_std"],
        "regionCondNregions": res_a["n_groups"],
        "regionCondNperms":   N_PERMS,
        "regionCondBPermP":   res_b["p_val"],
        "regionCondBPermZ":   res_b["z_val"],
        "regionCondBNblocks": res_b["n_groups"],
    })
    print("\n  Results written to store.")


if __name__ == "__main__":
    main()
