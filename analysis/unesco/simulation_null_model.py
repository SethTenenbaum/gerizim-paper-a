# VERY IMPORTANT: OUTPUT VALUES MUST WRITE TO DISC SINCE IT'S SO SLOW TO RUN!!!!!
# RERUNNING THIS TAKES A LONG TIME PLEASE WRITE VALUES TO DISC
import re
import sys
from pathlib import Path
from collections import Counter

import numpy as np
from scipy.stats import gaussian_kde, binomtest

np.random.seed(42)

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, TIER_APLUS, TIER_B_MAX, TIER_APP, P_NULL_AP
from lib.beru import deviation as _beru_deviation
N_PERMS    = 100_000  # number of permutation trials

# (from dome_filter.py, shared with spherical_monument_test.py)
from lib.dome_filter import is_dome_site
from lib.results_store import ResultsStore

DOME_EXCLUDE = set()  # no exclusions (see spherical_monument_test.py for rationale)

def beru_deviation(lon: float) -> float:
    return _beru_deviation(lon)

def is_aplus(dev: float) -> bool:
    return dev <= TIER_APLUS

def count_aplus(lons):
    return sum(1 for lon in lons if is_aplus(beru_deviation(lon)))

def count_app(lons):
    return sum(1 for lon in lons if beru_deviation(lon) <= TIER_APP)

# ── Load UNESCO data ────────────────────────────────────────────────────────
def load_unesco():
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)

    sites = []
    for s in cultural:
        # consistent with spherical_monument_test.py
        is_dome = is_dome_site(s) if s.site not in DOME_EXCLUDE else False

        sites.append({
            "name": s.site,
            "lon": s.longitude,
            "lat": s.latitude,
            "year": s.year,
            "is_dome": is_dome,
        })

    return sites

# ── Main ────────────────────────────────────────────────────────────────────
def main():
    print("=" * 90)
    print("  SIMULATION-BASED NULL MODEL FOR BERU-GRID HYPOTHESIS")
    print(f"  Anchor: {GERIZIM}°E  |  BERU: {BERU}°  |  N_perms: {N_PERMS:,}")
    print("=" * 90)

    sites = load_unesco()
    N = len(sites)
    lons = np.array([s["lon"] for s in sites])

    # Observed values
    obs_ap = count_aplus(lons)
    obs_app = count_app(lons)
    obs_rate = obs_ap / N * 100

    # Dome sites
    dome_sites = [s for s in sites if s["is_dome"]]
    dome_lons = np.array([s["lon"] for s in dome_sites])
    N_dome = len(dome_sites)
    obs_dome_ap = count_aplus(dome_lons)

    canon_sites = [s for s in sites if s["year"] is not None and 1978 <= s["year"] <= 1984]
    canon_lons = np.array([s["lon"] for s in canon_sites])
    N_canon = len(canon_sites)
    obs_canon_ap = count_aplus(canon_lons)

    # Pre-2000 sites
    pre2k_sites = [s for s in sites if s["year"] is not None and s["year"] < 2000]
    pre2k_lons = np.array([s["lon"] for s in pre2k_sites])
    N_pre2k = len(pre2k_sites)
    obs_pre2k_ap = count_aplus(pre2k_lons)

    # Post-2000 sites
    post2k_sites = [s for s in sites if s["year"] is not None and s["year"] >= 2000]
    post2k_lons = np.array([s["lon"] for s in post2k_sites])
    N_post2k = len(post2k_sites)
    obs_post2k_ap = count_aplus(post2k_lons)

    print(f"\n  Full corpus: N = {N},  A+ = {obs_ap} ({obs_rate:.1f}%)")
    print(f"  Dome/spherical: N = {N_dome},  A+ = {obs_dome_ap} ({obs_dome_ap/N_dome*100:.1f}%)")
    print(f"  Canon (1978-84): N = {N_canon},  A+ = {obs_canon_ap} ({obs_canon_ap/N_canon*100:.1f}%)")
    print(f"  Pre-2000: N = {N_pre2k},  A+ = {obs_pre2k_ap} ({obs_pre2k_ap/N_pre2k*100:.1f}%)")
    print(f"  Post-2000: N = {N_post2k},  A+ = {obs_post2k_ap} ({obs_post2k_ap/N_post2k*100:.1f}%)")

    print("\n" + "─" * 90)
    print("  APPROACH 1: LONGITUDE PERMUTATION (shuffle longitudes among sites)")
    print("─" * 90)

    perm_ap_counts = np.zeros(N_PERMS, dtype=int)
    perm_dome_ap = np.zeros(N_PERMS, dtype=int)
    perm_canon_ap = np.zeros(N_PERMS, dtype=int)
    perm_pre2k_ap = np.zeros(N_PERMS, dtype=int)
    perm_post2k_ap = np.zeros(N_PERMS, dtype=int)

    dome_idxs = [i for i, s in enumerate(sites) if s["is_dome"]]
    canon_idxs = [i for i, s in enumerate(sites) if s["year"] is not None and 1978 <= s["year"] <= 1984]
    pre2k_idxs = [i for i, s in enumerate(sites) if s["year"] is not None and s["year"] < 2000]
    post2k_idxs = [i for i, s in enumerate(sites) if s["year"] is not None and s["year"] >= 2000]

    for trial in range(N_PERMS):
        shuffled = np.random.permutation(lons)
        ap_mask = np.array([is_aplus(beru_deviation(l)) for l in shuffled])
        perm_ap_counts[trial] = ap_mask.sum()
        perm_dome_ap[trial] = ap_mask[dome_idxs].sum()
        perm_canon_ap[trial] = ap_mask[canon_idxs].sum()
        perm_pre2k_ap[trial] = ap_mask[pre2k_idxs].sum()
        perm_post2k_ap[trial] = ap_mask[post2k_idxs].sum()

        if (trial + 1) % 10000 == 0:
            print(f"    ... {trial+1:,} / {N_PERMS:,} permutations", file=sys.stderr)

    p_perm_full = np.mean(perm_ap_counts >= obs_ap)
    mean_perm_ap = perm_ap_counts.mean()
    std_perm_ap = perm_ap_counts.std()
    z_perm = (obs_ap - mean_perm_ap) / std_perm_ap if std_perm_ap > 0 else 0

    print(f"\n  Full corpus (N={N}):")
    print(f"    Observed A+ = {obs_ap}")
    print(f"    Permutation mean A+ = {mean_perm_ap:.1f} ± {std_perm_ap:.1f}")
    print(f"    Permutation Z = {z_perm:.2f}")
    print(f"    p_perm (≥ {obs_ap}) = {p_perm_full:.6f}")
    print(f"    → {'SIGNIFICANT' if p_perm_full < 0.05 else 'NOT significant'} at α=0.05")

    if std_perm_ap < 0.01:
        print(f"\n    *** NOTE: Full-corpus A+ count is invariant under permutation")
        print(f"    *** (shuffling longitudes doesn't change total A+ count)")
        print(f"    *** The meaningful tests are for SUB-POPULATIONS below:")

    # Dome sub-population
    p_perm_dome = np.mean(perm_dome_ap >= obs_dome_ap)
    mean_dome = perm_dome_ap.mean()
    std_dome = perm_dome_ap.std()
    z_dome = (obs_dome_ap - mean_dome) / std_dome if std_dome > 0 else 0

    print(f"\n  Dome/spherical (N={N_dome}):")
    print(f"    Observed A+ = {obs_dome_ap}")
    print(f"    Permutation mean A+ = {mean_dome:.2f} ± {std_dome:.2f}")
    print(f"    Permutation Z = {z_dome:.2f}")
    print(f"    p_perm (≥ {obs_dome_ap}) = {p_perm_dome:.6f}")
    print(f"    → {'SIGNIFICANT' if p_perm_dome < 0.05 else 'NOT significant'} at α=0.05")

    # Canon sub-population
    p_perm_canon = np.mean(perm_canon_ap >= obs_canon_ap)
    mean_canon = perm_canon_ap.mean()
    std_canon = perm_canon_ap.std()
    z_canon = (obs_canon_ap - mean_canon) / std_canon if std_canon > 0 else 0

    print(f"\n  Founding canon (1978-84, N={N_canon}):")
    print(f"    Observed A+ = {obs_canon_ap}")
    print(f"    Permutation mean A+ = {mean_canon:.2f} ± {std_canon:.2f}")
    print(f"    Permutation Z = {z_canon:.2f}")
    print(f"    p_perm (≥ {obs_canon_ap}) = {p_perm_canon:.6f}")
    print(f"    → {'SIGNIFICANT' if p_perm_canon < 0.05 else 'NOT significant'} at α=0.05")

    # Pre-2000 sub-population
    p_perm_pre2k = np.mean(perm_pre2k_ap >= obs_pre2k_ap)
    mean_pre2k = perm_pre2k_ap.mean()
    std_pre2k = perm_pre2k_ap.std()
    z_pre2k = (obs_pre2k_ap - mean_pre2k) / std_pre2k if std_pre2k > 0 else 0

    print(f"\n  Pre-2000 (N={N_pre2k}):")
    print(f"    Observed A+ = {obs_pre2k_ap}")
    print(f"    Permutation mean A+ = {mean_pre2k:.2f} ± {std_pre2k:.2f}")
    print(f"    Permutation Z = {z_pre2k:.2f}")
    print(f"    p_perm (≥ {obs_pre2k_ap}) = {p_perm_pre2k:.6f}")

    # Post-2000 sub-population
    p_perm_post2k = np.mean(perm_post2k_ap >= obs_post2k_ap)
    mean_post2k = perm_post2k_ap.mean()
    std_post2k = perm_post2k_ap.std()
    z_post2k = (obs_post2k_ap - mean_post2k) / std_post2k if std_post2k > 0 else 0

    print(f"\n  Post-2000 (N={N_post2k}):")
    print(f"    Observed A+ = {obs_post2k_ap}")
    print(f"    Permutation mean A+ = {mean_post2k:.2f} ± {std_post2k:.2f}")
    print(f"    Permutation Z = {z_post2k:.2f}")
    print(f"    p_perm (≥ {obs_post2k_ap}) = {p_perm_post2k:.6f}")

    print("\n" + "─" * 90)
    print("  APPROACH 2: BOOTSTRAP FROM EMPIRICAL LONGITUDE DISTRIBUTION")
    print("  Draw N longitudes with replacement from observed; count A+")
    print("─" * 90)

    boot_ap_counts = np.zeros(N_PERMS, dtype=int)
    for trial in range(N_PERMS):
        sample = np.random.choice(lons, size=N, replace=True)
        boot_ap_counts[trial] = count_aplus(sample)
        if (trial + 1) % 10000 == 0:
            print(f"    ... {trial+1:,} / {N_PERMS:,} bootstrap samples", file=sys.stderr)

    p_boot = np.mean(boot_ap_counts >= obs_ap)
    boot_mean = boot_ap_counts.mean()
    boot_std = boot_ap_counts.std()
    z_boot = (obs_ap - boot_mean) / boot_std if boot_std > 0 else 0

    print(f"\n  Bootstrap (N={N}, {N_PERMS:,} samples):")
    print(f"    Observed A+ = {obs_ap}")
    print(f"    Bootstrap mean A+ = {boot_mean:.1f} ± {boot_std:.1f}")
    print(f"    Bootstrap Z = {z_boot:.2f}")
    print(f"    p_boot (≥ {obs_ap}) = {p_boot:.6f}")

    print("\n" + "─" * 90)
    print("  APPROACH 3: KDE-SMOOTHED NULL")
    print("  Fit KDE to observed longitudes; draw N samples; count A+")
    print("─" * 90)

    kde = gaussian_kde(lons, bw_method=0.05)  # bandwidth ~1.8° at this data range
    kde_ap_counts = np.zeros(N_PERMS, dtype=int)
    for trial in range(N_PERMS):
        sample = kde.resample(N).flatten()
        kde_ap_counts[trial] = count_aplus(sample)
        if (trial + 1) % 10000 == 0:
            print(f"    ... {trial+1:,} / {N_PERMS:,} KDE samples", file=sys.stderr)

    p_kde = np.mean(kde_ap_counts >= obs_ap)
    kde_mean = kde_ap_counts.mean()
    kde_std = kde_ap_counts.std()
    z_kde = (obs_ap - kde_mean) / kde_std if kde_std > 0 else 0

    print(f"\n  KDE-smoothed (N={N}, {N_PERMS:,} samples, bw=0.05):")
    print(f"    Observed A+ = {obs_ap}")
    print(f"    KDE mean A+ = {kde_mean:.1f} ± {kde_std:.1f}")
    print(f"    KDE Z = {z_kde:.2f}")
    print(f"    p_kde (≥ {obs_ap}) = {p_kde:.6f}")

    print("\n" + "─" * 90)
    print("  APPROACH 2b: BOOTSTRAP FOR DOME/SPHERICAL SUB-POPULATION")
    print(f"  Draw {N_dome} longitudes from the full corpus; count A+")
    print("─" * 90)

    boot_dome_counts = np.zeros(N_PERMS, dtype=int)
    for trial in range(N_PERMS):
        sample = np.random.choice(lons, size=N_dome, replace=True)
        boot_dome_counts[trial] = count_aplus(sample)
        if (trial + 1) % 10000 == 0:
            print(f"    ... {trial+1:,} / {N_PERMS:,} dome bootstrap", file=sys.stderr)

    p_boot_dome = np.mean(boot_dome_counts >= obs_dome_ap)
    boot_dome_mean = boot_dome_counts.mean()
    boot_dome_std = boot_dome_counts.std()
    z_boot_dome = (obs_dome_ap - boot_dome_mean) / boot_dome_std if boot_dome_std > 0 else 0

    print(f"\n  Dome/spherical bootstrap (N={N_dome}, {N_PERMS:,} samples):")
    print(f"    Observed A+ = {obs_dome_ap}")
    print(f"    Bootstrap mean A+ = {boot_dome_mean:.2f} ± {boot_dome_std:.2f}")
    print(f"    Bootstrap Z = {z_boot_dome:.2f}")
    print(f"    p_boot_dome (≥ {obs_dome_ap}) = {p_boot_dome:.6f}")

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 90)
    print("  SUMMARY: SIMULATION-BASED NULL MODEL RESULTS")
    print("=" * 90)
    print(f"  {'Test':<40} {'Obs':>5} {'Sim mean':>10} {'Z':>7} {'p_sim':>10} {'Sig':>5}")
    print(f"  {'-'*40} {'-'*5:>5} {'-'*10:>10} {'-'*7:>7} {'-'*10:>10} {'-'*5:>5}")

    results = [
        ("Full corpus (permutation)", obs_ap, mean_perm_ap, z_perm, p_perm_full),
        ("Full corpus (bootstrap)", obs_ap, boot_mean, z_boot, p_boot),
        ("Full corpus (KDE)", obs_ap, kde_mean, z_kde, p_kde),
        ("Dome/spher. (permutation)", obs_dome_ap, mean_dome, z_dome, p_perm_dome),
        ("Dome/spher. (bootstrap)", obs_dome_ap, boot_dome_mean, z_boot_dome, p_boot_dome),
        ("Canon 1978-84 (permutation)", obs_canon_ap, mean_canon, z_canon, p_perm_canon),
        ("Pre-2000 (permutation)", obs_pre2k_ap, mean_pre2k, z_pre2k, p_perm_pre2k),
        ("Post-2000 (permutation)", obs_post2k_ap, mean_post2k, z_post2k, p_perm_post2k),
    ]

    for name, obs, mean, z, p in results:
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "~" if p < 0.10 else "ns"
        print(f"  {name:<40} {obs:5d} {mean:10.2f} {z:7.2f} {p:10.6f} {sig:>5}")

    # ── Interpretation ───────────────────────────────────────────────────────
    print("\n" + "─" * 90)
    print("─" * 90)
    print("""
  (including the Europe-heavy clustering) but DESTROYS any association
  between a site's cultural identity and its beru-grid position.

  permutation (shuffling longitudes doesn't change the total count).
  The meaningful tests are for SUB-POPULATIONS:

  - If dome/spherical monuments have more A+ than expected under permutation,
    it means dome sites are disproportionately at harmonic longitudes —

  - If the founding canon (1978-84) has more A+ than expected under
    permutation, it means early-inscribed sites are disproportionately
    at harmonic longitudes — the temporal gradient survives the
    non-uniform null.

  is anomalous given the empirical longitude clumping.
""")

    # ── LaTeX macros (GROUP 25) ──────────────────────────────────────────────
    print("  % LaTeX macros (GROUP 25):")
    _nperms_fmt = f"{N_PERMS:,}".replace(",", "{,}")
    print(f"  \\newcommand{{\\simNperms}}{{{_nperms_fmt}}}        % number of permutations in null model")
    print(f"  \\newcommand{{\\simDomePermZ}}{{{z_dome:.2f}}}         % permutation Z-score, dome corpus")
    print(f"  \\newcommand{{\\simDomePermP}}{{{p_perm_dome:.3f}}}          % permutation p, dome corpus")
    print(f"  \\newcommand{{\\simDomeBootP}}{{{p_boot_dome:.3f}}}          % bootstrap p, dome corpus")
    print(f"  \\newcommand{{\\simCanonPermP}}{{{p_perm_canon:.3f}}}         % permutation p, canonical sites")
    print(f"  \\newcommand{{\\simPreTwoKPermP}}{{{p_perm_pre2k:.3f}}}        % permutation p, pre-2000 BCE sites")
    print(f"  \\newcommand{{\\simPostTwoKPermP}}{{{p_perm_post2k:.3f}}}       % permutation p, post-2000 BCE sites")
    print(f"  \\newcommand{{\\simKDEZ}}{{{z_kde:.2f}}}           % KDE-based Z-score, null model")
    print(f"  \\newcommand{{\\simKDEP}}{{{p_kde:.3f}}}           % KDE-based p-value, null model")
    print(f"  \\newcommand{{\\simBootFullMean}}{{{boot_mean:.1f}}}         % bootstrap mean A+, full corpus")
    print(f"  \\newcommand{{\\simBootFullStd}}{{{boot_std:.1f}}}          % bootstrap std A+, full corpus")
    print(f"  \\newcommand{{\\simKDEMean}}{{{kde_mean:.1f}}}           % KDE null-model mean A+")
    print(f"  \\newcommand{{\\simKDEStd}}{{{kde_std:.1f}}}            % KDE null-model std A+")

    # ── Write to results store ────────────────────────────────────────────────
    ResultsStore().write_many({
        "simDomePermP":      p_perm_dome,    # permutation p, dome corpus
        "simDomeBootP":      p_boot_dome,    # bootstrap p, dome corpus
        "simKDEP":           p_kde,          # KDE-based p-value
        "simCanonPermP":     p_perm_canon,   # permutation p, canonical
        "simPreTwoKPermP":   p_perm_pre2k,   # permutation p, pre-2000
        "simPostTwoKPermP":  p_perm_post2k,  # permutation p, post-2000
    })

if __name__ == "__main__":
    main()
