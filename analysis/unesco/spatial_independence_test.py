import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
from scipy.stats import binomtest

np.random.seed(42)

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, TIER_APLUS, P_NULL_AP
from lib.beru import deviation as _beru_deviation
from lib.results_store import ResultsStore

N_BOOT     = 100_000

def beru_deviation(lon: float) -> float:
    return _beru_deviation(lon)

def is_aplus(lon: float) -> bool:
    return beru_deviation(lon) <= TIER_APLUS

# ── Load UNESCO data ────────────────────────────────────────────────────────
def load_sites():
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)

    sites = []
    for s in cultural:
        sites.append({
            "name": s.site,
            "lon": s.longitude,
            "lat": s.latitude,
            "year": s.year,
            "is_ap": is_aplus(s.longitude),
        })
    return sites

def cluster_effective_n(sites, cluster_radius=0.5):
    sorted_sites = sorted(sites, key=lambda s: s["lon"])
    clusters = []
    current_cluster = [sorted_sites[0]]

    for s in sorted_sites[1:]:
        if abs(s["lon"] - current_cluster[0]["lon"]) <= cluster_radius * 2:
            current_cluster.append(s)
        else:
            clusters.append(current_cluster)
            current_cluster = [s]
    clusters.append(current_cluster)

    n_clusters = len(clusters)
    n_ap_clusters = sum(1 for c in clusters if any(s["is_ap"] for s in c))

    n_ap_majority = sum(
        1 for c in clusters
        if sum(1 for s in c if s["is_ap"]) > len(c) / 2
    )

    weighted_ap = sum(
        sum(1 for s in c if s["is_ap"]) / len(c)
        for c in clusters
    )

    return {
        "n_clusters": n_clusters,
        "n_ap_any": n_ap_clusters,
        "n_ap_majority": n_ap_majority,
        "weighted_ap": weighted_ap,
        "clusters": clusters,
    }

def design_effect_correction(sites, bandwidth=0.5):
    # Form clusters
    sorted_sites = sorted(sites, key=lambda s: s["lon"])
    clusters = []
    current = [sorted_sites[0]]
    for s in sorted_sites[1:]:
        if abs(s["lon"] - current[0]["lon"]) <= bandwidth * 2:
            current.append(s)
        else:
            clusters.append(current)
            current = [s]
    clusters.append(current)

    N = len(sites)
    K = len(clusters)  # number of clusters
    m_sizes = [len(c) for c in clusters]
    m_bar = N / K  # average cluster size

    p_hat = sum(s["is_ap"] for s in sites) / N

    cluster_means = []
    ssw = 0.0
    for c in clusters:
        c_mean = sum(1 for s in c if s["is_ap"]) / len(c)
        cluster_means.append(c_mean)
        for s in c:
            ssw += (int(s["is_ap"]) - c_mean) ** 2

    ssb = sum(len(c) * (cm - p_hat) ** 2 for c, cm in zip(clusters, cluster_means))

    MSB = ssb / max(K - 1, 1)
    MSW = ssw / max(N - K, 1)

    # ICC (rho)
    m0 = (N - sum(m ** 2 for m in m_sizes) / N) / max(K - 1, 1)
    rho = max(0, (MSB - MSW) / (MSB + (m0 - 1) * MSW)) if (MSB + (m0 - 1) * MSW) > 0 else 0

    # Design effect
    DEFF = 1 + (m_bar - 1) * rho
    N_eff = N / DEFF

    return {
        "N": N,
        "K": K,
        "m_bar": m_bar,
        "rho": rho,
        "DEFF": DEFF,
        "N_eff": N_eff,
    }

def spatial_block_bootstrap(sites, block_width=3.0, n_boot=N_BOOT):
    blocks = defaultdict(list)
    for s in sites:
        block_id = int(np.floor(s["lon"] / block_width))
        blocks[block_id].append(s)

    block_ids = list(blocks.keys())
    n_blocks = len(block_ids)

    # Observed
    obs_ap = sum(1 for s in sites if s["is_ap"])
    obs_N = len(sites)

    # Bootstrap
    boot_ap_rates = np.zeros(n_boot)
    for trial in range(n_boot):
        sampled_ids = np.random.choice(block_ids, size=n_blocks, replace=True)
        sample_ap = 0
        sample_n = 0
        for bid in sampled_ids:
            block = blocks[bid]
            sample_ap += sum(1 for s in block if s["is_ap"])
            sample_n += len(block)
        boot_ap_rates[trial] = sample_ap / sample_n if sample_n > 0 else 0

        if (trial + 1) % 10000 == 0:
            print(f"    ... {trial+1:,} / {n_boot:,} block bootstrap", file=sys.stderr)

    obs_rate = obs_ap / obs_N
    p_block = np.mean(boot_ap_rates >= obs_rate)

    ci_lo = np.percentile(boot_ap_rates, 2.5)
    ci_hi = np.percentile(boot_ap_rates, 97.5)
    null_excluded = ci_lo > P_NULL_AP

    return {
        "n_blocks": n_blocks,
        "block_width": block_width,
        "obs_rate": obs_rate,
        "boot_mean": boot_ap_rates.mean(),
        "boot_std": boot_ap_rates.std(),
        "p_block": p_block,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "null_excluded": null_excluded,
        "boot_rates": boot_ap_rates,
    }

# ── Main ────────────────────────────────────────────────────────────────────
def main():
    print("=" * 90)
    print("  SPATIAL INDEPENDENCE TEST: EFFECTIVE-N AND BLOCK BOOTSTRAP")
    print(f"  Anchor: {GERIZIM}°E  |  BERU: {BERU}°  |  A+ threshold: {TIER_APLUS}")
    print("=" * 90)

    sites = load_sites()
    N = len(sites)
    obs_ap = sum(1 for s in sites if s["is_ap"])
    obs_rate = obs_ap / N

    print(f"\n  Raw: N = {N}, A+ = {obs_ap} ({obs_rate*100:.1f}%)")
    print(f"  Binomial p (raw, assuming independence) = "
          f"{binomtest(obs_ap, N, P_NULL_AP, alternative='greater').pvalue:.6f}")

    # ── Approach 1: Cluster-level ────────────────────────────────────────────
    print("\n" + "─" * 90)
    print("  APPROACH 1: CLUSTER-LEVEL EFFECTIVE-N")
    print("─" * 90)

    for radius in [0.25, 0.50, 1.00]:
        c = cluster_effective_n(sites, cluster_radius=radius)
        p_any = binomtest(c["n_ap_any"], c["n_clusters"], P_NULL_AP,
                          alternative='greater').pvalue
        p_maj = binomtest(c["n_ap_majority"], c["n_clusters"], P_NULL_AP,
                          alternative='greater').pvalue

        print(f"\n  Cluster radius: ±{radius}°")
        print(f"    N_clusters = {c['n_clusters']}")
        print(f"    A+ clusters (any member): {c['n_ap_any']} "
              f"({c['n_ap_any']/c['n_clusters']*100:.1f}%)")
        print(f"    p (any-member): {p_any:.6f}"
              f"  {'*' if p_any < 0.05 else 'ns'}")
        print(f"    A+ clusters (majority): {c['n_ap_majority']} "
              f"({c['n_ap_majority']/c['n_clusters']*100:.1f}%)")
        print(f"    p (majority): {p_maj:.6f}"
              f"  {'*' if p_maj < 0.05 else 'ns'}")
        print(f"    Weighted A+ (fractional): {c['weighted_ap']:.1f} "
              f"({c['weighted_ap']/c['n_clusters']*100:.1f}%)")

    print("\n" + "─" * 90)
    print("  APPROACH 2: DESIGN EFFECT / ICC CORRECTION")
    print("─" * 90)

    for bw in [0.25, 0.50, 1.00]:
        d = design_effect_correction(sites, bandwidth=bw)
        n_eff_int = max(1, int(round(d["N_eff"])))
        ap_eff = max(1, int(round(obs_ap * d["N_eff"] / d["N"])))
        p_eff = binomtest(ap_eff, n_eff_int, P_NULL_AP,
                          alternative='greater').pvalue

        from scipy.stats import norm as _norm
        z_eff = (obs_rate - P_NULL_AP) / np.sqrt(P_NULL_AP * (1 - P_NULL_AP) / d["N_eff"])
        p_z   = 1 - _norm.cdf(z_eff)

        print(f"\n  Bandwidth: ±{bw}°")
        print(f"    K = {d['K']} clusters, m̄ = {d['m_bar']:.1f}")
        print(f"    ICC (ρ) = {d['rho']:.4f}")
        print(f"    DEFF = {d['DEFF']:.3f}")
        print(f"    N_eff = {d['N_eff']:.0f}")
        print(f"    p (N_eff-corrected, binomial rounded) = {p_eff:.6f}"
              f"  {'*' if p_eff < 0.05 else 'ns'}")
        print(f"    p (N_eff-corrected, Z-approx exact)   = {p_z:.6f}"
              f"  Z = {z_eff:.2f}  {'*' if p_z < 0.05 else 'ns'}")

    print("\n" + "─" * 90)
    print("  APPROACH 3: SPATIAL BLOCK BOOTSTRAP")
    print("─" * 90)

    for block_w in [3.0, 5.0, 10.0]:
        b = spatial_block_bootstrap(sites, block_width=block_w, n_boot=N_BOOT)
        z = (b["obs_rate"] - P_NULL_AP) / b["boot_std"] if b["boot_std"] > 0 else 0
        from scipy.stats import norm as _norm
        p_null_z = 1 - _norm.cdf(z)

        print(f"\n  Block width: {block_w}°")
        print(f"    N_blocks = {b['n_blocks']}")
        print(f"    Observed rate = {b['obs_rate']*100:.2f}%")
        print(f"    Bootstrap mean rate = {b['boot_mean']*100:.2f}% ± {b['boot_std']*100:.2f}%")
        print(f"    95% CI = [{b['ci_lo']*100:.2f}%, {b['ci_hi']*100:.2f}%]")
        print(f"    {P_NULL_AP:.1%} null excluded from CI: {'YES' if b['null_excluded'] else 'NO'}")
        print(f"    Z (observed rate vs {P_NULL_AP:.1%} null, SE from bootstrap) = {z:.2f}")
        print(f"    p (Z-test vs {P_NULL_AP:.1%} null) = {p_null_z:.6f}"
              f"  {'*' if p_null_z < 0.05 else 'ns'}")
        sig = "YES ★" if b['null_excluded'] else "NO"
        print(f"    → Spatial-dependence-adjusted significance: {sig}")

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 90)
    print("  SUMMARY")
    print("=" * 90)
    print("""
  The raw binomial test (N=1011) assumes independence; this script
  tests whether spatial clustering inflates the result.

  Three corrections:
  1. Cluster-level: collapse nearby sites into clusters and test
     clusters, not individual sites.
  2. Design effect (DEFF): compute effective-N via ICC of A+ status
     within longitude clusters.
  3. Spatial block bootstrap: resample longitude blocks with
     replacement, preserving within-block dependence. The 95% CI
     on the A+ rate accounts for spatial autocorrelation; if the
     lower bound exceeds the geometric null (P_NULL_AP), the enrichment is confirmed even
     under dependence.

  If the independence assumption were driving the result, we would
  expect the A+ rate to be consistent across space. The methods
  above test and account for spatial structure. If the 95% CI
  excludes the geometric null (P_NULL_AP) from the CI, the independence assumption is not driving
  the result.
""")

    # ── LaTeX macros (GROUP 17) ──────────────────────────────────────────────
    # DEFF at ±0.5°
    d_half = design_effect_correction(sites, bandwidth=0.50)
    n_eff_half = max(1, int(round(d_half["N_eff"])))
    ap_eff_half = max(1, int(round(obs_ap * d_half["N_eff"] / d_half["N"])))
    p_eff_half = binomtest(ap_eff_half, n_eff_half, P_NULL_AP, alternative='greater').pvalue

    # DEFF at ±0.25°
    d_qtr = design_effect_correction(sites, bandwidth=0.25)
    n_eff_qtr = max(1, int(round(d_qtr["N_eff"])))
    ap_eff_qtr = max(1, int(round(obs_ap * d_qtr["N_eff"] / d_qtr["N"])))
    p_eff_qtr = binomtest(ap_eff_qtr, n_eff_qtr, P_NULL_AP, alternative='greater').pvalue

    # Block bootstrap at 5.0°
    b5 = spatial_block_bootstrap(sites, block_width=5.0, n_boot=N_BOOT)
    z5 = (b5["obs_rate"] - P_NULL_AP) / b5["boot_std"] if b5["boot_std"] > 0 else 0

    print("  % LaTeX macros (GROUP 17):")
    print(f"  \\newcommand{{\\DEFFhalf}}{{{d_half['DEFF']:.3f}}}          % design effect (DEFF), half-degree blocks")
    print(f"  \\newcommand{{\\NeffHalf}}{{{n_eff_half}}}            % effective N, half-degree blocks")
    print(f"  \\newcommand{{\\pNeffHalf}}{{{p_eff_half:.3f}}}           % p-value, A+ test with effective N (half-°)")
    print(f"  \\newcommand{{\\ICChalf}}{{{d_half['rho']:.3f}}}           % intracluster correlation, half-degree blocks")
    print(f"  \\newcommand{{\\DEFFquarter}}{{{d_qtr['DEFF']:.3f}}}         % design effect (DEFF), quarter-degree blocks")
    print(f"  \\newcommand{{\\NeffQuarter}}{{{n_eff_qtr}}}          % effective N, quarter-degree blocks")
    print(f"  \\newcommand{{\\pNeffQuarter}}{{{p_eff_qtr:.3f}}}         % p-value, A+ test with effective N (quarter-°)")
    print(f"  \\newcommand{{\\blockBootCIlo}}{{{b5['ci_lo']*100:.2f}}}        % block-bootstrap CI lower bound, A+ rate (%)")
    print(f"  \\newcommand{{\\blockBootCIhi}}{{{b5['ci_hi']*100:.2f}}}        % block-bootstrap CI upper bound, A+ rate (%)")
    print(f"  \\newcommand{{\\blockBootZ}}{{{z5:.2f}}}           % block-bootstrap Z-score, A+")

    # ── Write to results store ────────────────────────────────────────────────
    from scipy.stats import norm as _norm
    blockBootZ_p = _norm.sf(z5)   # one-sided p from Z-score
    ResultsStore().write_many({
        "pNeffQuarter":  p_eff_qtr,   # DEFF-corrected binomial p (quarter-degree)
        "pNeffHalf":     p_eff_half,  # DEFF-corrected binomial p (half-degree)
        "blockBootZ_p":  blockBootZ_p, # one-sided p from block-bootstrap Z
    })

if __name__ == "__main__":
    main()
