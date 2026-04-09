
import numpy as np

def benjamini_hochberg(pvals, alpha=0.05):
    m = len(pvals)
    sorted_idx = np.argsort(pvals)
    sorted_pvals = np.array(pvals)[sorted_idx]

    # Adjusted p-values
    q_values = np.zeros(m)
    q_values[m - 1] = sorted_pvals[m - 1]
    for i in range(m - 2, -1, -1):
        q_values[i] = min(q_values[i + 1], sorted_pvals[i] * m / (i + 1))

    q_out = np.zeros(m)
    for i, idx in enumerate(sorted_idx):
        q_out[idx] = q_values[i]

    return q_out

def holm_correction(pvals):
    m = len(pvals)
    sorted_idx = np.argsort(pvals)
    sorted_pvals = np.array(pvals)[sorted_idx]

    adj = np.zeros(m)
    adj[0] = sorted_pvals[0] * m
    for i in range(1, m):
        adj[i] = max(adj[i - 1], sorted_pvals[i] * (m - i))
    adj = np.minimum(adj, 1.0)

    out = np.zeros(m)
    for i, idx in enumerate(sorted_idx):
        out[idx] = adj[i]
    return out

def main():
    # Format: (name, p_value, category)

    tests = [
        # are listed separately below.
        ("Test 1: Dome/spherical A+ (binomial)", 0.0009, "C"),
        ("Test 1: Dome/spherical A (binomial)", 0.0043, "C"),
        ("Test 1: Dome χ²-uniform", 0.1266, "C"),
        ("Test 2/4: Full corpus A+ (binomial)", 0.0103, "C"),
        ("Test 2: Cluster MWU (one-tailed)", 0.0182, "C"),
        ("Test 2: Cluster permutation", 0.0033, "C"),
        ("Test 2: Harmonic-level MWU", 0.0001, "C"),
        ("Test 3: Buddhist A+ (binomial)", 0.5299, "C"),
        ("Test 3: Buddhist A (binomial)", 0.1256, "C"),

        # ── Temporal gradient tests ──
        ("Cochran-Armitage 3-cohort", 0.0481, "C"),
        ("Cochran-Armitage 5-cohort", 0.0489, "C"),
        ("Canon (1978-84) A+ binomial", 0.0232, "C"),
        ("Pre-2000 A+ binomial", 0.0041, "C"),
        ("Modern (2000+) A+ binomial", 0.3081, "C"),
        ("Fisher canon vs modern", 0.0844, "C"),
        ("First-half A+ binomial", 0.0046, "C"),
        ("Second-half A+ binomial", 0.2952, "C"),
        ("Fisher first-half vs second-half", 0.0849, "C"),

        # ── Sequential significance (peak) ──
        ("Sequential peak (N=605)", 0.00026, "E"),  # selected as best

        # ── 5-cohort individual p-values ──
        ("Cohort 1985-91 binomial", 0.205, "E"),
        ("Cohort 1992-99 binomial", 0.035, "E"),
        ("Cohort 2000-09 binomial", 0.338, "E"),
        ("Cohort 2010-25 binomial", 0.419, "E"),

        # ── Religion keyword tests ──
        ("Any-religion A+ binomial", 0.0369, "C"),
        ("Christian A+ binomial", 0.0312, "E"),

        # ── Founding-era dating ──
        ("Pre-CE A+ binomial", 0.2577, "E"),
        ("Early modern A+ binomial", 0.0303, "E"),
        ("Canon × post-CE interaction", 0.0012, "E"),

        # ── Founding keyword enrichment ──
        ("Founding keyword Fisher exact", 0.154, "C"),

        ("Sacred Origin: A+ binomial (geo null)", 6.55e-12, "C"),
        ("Sacred Origin: A++ binomial (geo null)", 3.6e-6, "C"),
        ("Sacred Origin: A+ Fisher (vs rest)", 1.3e-7, "C"),
        ("Sacred Origin: A++ Fisher (vs rest)", 6e-6, "C"),
        ("Sacred Origin: A++ hypergeometric", 6e-6, "C"),
        ("Sacred Origin: A+ permutation (empirical)", 0.0, "C"),
        ("Sacred Origin: A++ permutation (empirical)", 1.0e-5, "C"),
        ("Sacred Origin: mean-delta permutation", 0.0, "C"),
        ("Sacred Origin: A++ within-A+ permutation", 0.0313, "E"),

        # ── Sliding-window blocks ──
        ("Block 1 (sites 1-100)", 0.0190, "E"),
        ("Block 3 (sites 201-300)", 0.0022, "E"),
        ("Block 5 (sites 401-500)", 0.0475, "E"),

        ("Spatial: DEFF-corrected binomial (0.25° band)", 0.014, "C"),
        ("Spatial: DEFF-corrected binomial (0.50° band)", 0.015, "C"),
        ("Spatial: block-bootstrap Z-test (3.0° blocks)", 0.012, "C"),

        ("CMH stratified test (early vs late, by region)", 0.082, "C"),

        ("Deep temporal: Spearman (founding date vs A+)", 0.941, "E"),
        ("Deep temporal: Fisher (pre-CE vs post-CE rate)", 0.697, "E"),

        # ── Simulation null model ──
        ("Sim: Dome permutation", 0.004, "C"),
        ("Sim: Dome bootstrap", 0.006, "C"),
        ("Sim: KDE full corpus", 0.011, "C"),
        ("Sim: Canon permutation", 0.127, "C"),
        ("Sim: Pre-2000 permutation", 0.094, "C"),
        ("Sim: Post-2000 permutation", 0.944, "C"),

        ("Global sweep: S1 megalith A+ binomial (N=16)", 0.0242, "E"),

        ("Regional: Asia & Pacific A+ binomial", 0.0477, "E"),
        ("Regional: Europe & N. America A+ binomial", 0.0234, "E"),
        ("Regional: Arab States Cochran-Armitage trend", 0.0347, "E"),

        # ── Anchor scan tests (analysis/anchor_uniqueness_audit.py,
        #                       analysis/landmark_anchor_ranking.py) ──
        ("Anchor scan: Gerizim rank 41/1011 (landmark scan)", 41/1011, "E"),
        ("Anchor scan: Gerizim ≥56 A+, 41/1011 (conservative)", 41/1011, "E"),

        ("Corridor: Occupancy binomial (17/17)", 1.72e-24, "C"),
        ("Corridor: Fisher combined p", 1.49e-19, "C"),
        ("Corridor: A++ enrichment (Fisher)", 0.0022, "C"),
        ("Corridor: Precision MW (corridor < global)", 0.016, "C"),

        ("Test 5: Evolution combined A+ (binomial)", 0.0009, "C"),
        ("Test 5: Stupa stage A+ (binomial)", 0.0050, "C"),
        ("Test 5: Dome stage A+ (binomial)", 0.0102, "C"),
        ("Test 5: Mound stage A+ (binomial)", 0.3012, "C"),  # pre-specified, ns as predicted
        ("Test 5: Fisher mound vs dome/stupa", 0.3826, "C"),  # two-sided, ns as predicted

        ("Test 6: Stupa A+ enrichment (binomial)", 0.5299, "C"),   # p_S ~ null rate; ns
        ("Test 6: Non-stupa A+ depletion (binomial)", 0.8306, "C"),# p_NS not depleted; ns
        ("Test 6: NS depletion Z-test (S>NS differential)", 0.6554, "C"),  # Z=-0.40; ns
        ("Test 6: Cluster enrichment (matched-size permutation)", 0.0, "E"),  # p=0.0000 *** (R=33.3km, beru-derived)
        ("Test 6: Cluster enrichment (rate-shuffle permutation)", 0.0462, "E"),  # p=0.046 * (R=33.3km)
        ("Test 6: Spacing uniqueness σ=0.101 diff Z", 0.0, "E"),   # Z=+3.97; p≈0; flagged exploratory

        # ── Anti-node tests (analysis/unesco_antinode_test.py,
        #                    analysis/full_antinode_analysis.py) ────────────
        ("Anti-node: enrichment binomial (anti-node A+>4%)", 0.0000, "C"),  # 49/500=9.8%; p≈0 ***
        ("Anti-node: depletion binomial (anti-node A+<4%)", 1.0000, "C"),   # not depleted; p=1.0 ns
        ("Anti-node: node>anti z-test (differential)",       0.2411, "C"),  # z=0.70; ns
        ("Anti-node: Fisher exact (node>anti, 2×2)",         0.2743, "C"),  # OR=1.16; ns
        ("Anti-node: KS distribution test",                  0.1406, "C"),  # D=0.07; ns
        ("Anti-node: joint node+anti enrichment (>7.84%)",   0.0016, "C"),  # 106/1011=10.5%; p=0.0016 **
        ("Anti-node: A++ enrichment binomial (>0.4%)",       0.0044, "C"),  # p=0.0044 **
        # Node founding=43.9%, anti founding=28.6%; diff=+15.3%.
        ("Anti-node: founding differential permutation",     0.0793, "E"),  # p=0.079; marginal
        ("Anti-node: block-bootstrap longitude correction",  1.0000, "C"),  # p=1.000; ns — artifact confirmed
        ("Anti-node: node>anti z-test (Tier-A, 0.010 beru)", 0.0875, "C"),  # z=1.36; p=0.088 ~
    ]

    names = [t[0] for t in tests]
    pvals = np.array([t[1] for t in tests])
    cats = [t[2] for t in tests]

    m = len(tests)

    # ── Corrections ──────────────────────────────────────────────────────────
    bonf = np.minimum(pvals * m, 1.0)
    holm = holm_correction(pvals)
    bh_q = benjamini_hochberg(pvals)

    conf_idx = [i for i, c in enumerate(cats) if c == "C"]
    expl_idx = [i for i, c in enumerate(cats) if c == "E"]

    conf_pvals = pvals[conf_idx]
    bh_conf = benjamini_hochberg(conf_pvals)

    print("=" * 110)
    print("  MULTIPLE COMPARISONS: FDR AND FAMILY-WISE CORRECTIONS")
    print(f"  Total tests: {m}  |  Confirmatory: {len(conf_idx)}  |  Exploratory: {len(expl_idx)}")
    print("=" * 110)

    # ── Table ────────────────────────────────────────────────────────────────
    header = f"  {'Test':<45} {'Cat':>3} {'p_raw':>9} {'Bonf':>9} {'Holm':>9} {'BH q':>9} {'Sig':>5}"
    print(f"\n{header}")
    print(f"  {'-'*45} {'-'*3:>3} {'-'*9:>9} {'-'*9:>9} {'-'*9:>9} {'-'*9:>9} {'-'*5:>5}")

    for i in range(m):
        sig = "***" if bh_q[i] < 0.001 else "**" if bh_q[i] < 0.01 else "*" if bh_q[i] < 0.05 else "~" if bh_q[i] < 0.10 else "ns"
        print(f"  {names[i]:<45} {cats[i]:>3} {pvals[i]:9.4f} {bonf[i]:9.4f} {holm[i]:9.4f} {bh_q[i]:9.4f} {sig:>5}")

    # ── Summary counts ───────────────────────────────────────────────────────
    print(f"\n" + "─" * 110)
    print(f"  SUMMARY: How many tests survive each correction at α = 0.05?")
    print(f"  {'Correction':<35} {'Significant':>12} {'of':>4} {m}")
    print(f"  {'-'*35} {'-'*12:>12} {'-'*4:>4}")
    print(f"  {'Raw (uncorrected)':<35} {sum(pvals < 0.05):>12}")
    print(f"  {'Bonferroni (k={})'.format(m):<35} {sum(bonf < 0.05):>12}")
    print(f"  {'Holm step-down (k={})'.format(m):<35} {sum(holm < 0.05):>12}")
    print(f"  {'Benjamini-Hochberg FDR (k={})'.format(m):<35} {sum(bh_q < 0.05):>12}")

    # ── Confirmatory subset ──────────────────────────────────────────────────
    print(f"\n" + "─" * 110)
    print(f"  CONFIRMATORY TESTS ONLY (pre-specified, N={len(conf_idx)})")
    print(f"  {'Correction':<35} {'Significant':>12} {'of':>4} {len(conf_idx)}")
    print(f"  {'-'*35} {'-'*12:>12} {'-'*4:>4}")
    bonf_conf = np.minimum(conf_pvals * len(conf_idx), 1.0)
    holm_conf = holm_correction(conf_pvals)
    print(f"  {'Raw (uncorrected)':<35} {sum(conf_pvals < 0.05):>12}")
    print(f"  {'Bonferroni (k={})'.format(len(conf_idx)):<35} {sum(bonf_conf < 0.05):>12}")
    print(f"  {'Holm step-down (k={})'.format(len(conf_idx)):<35} {sum(holm_conf < 0.05):>12}")
    print(f"  {'Benjamini-Hochberg FDR (k={})'.format(len(conf_idx)):<35} {sum(bh_conf < 0.05):>12}")

    print(f"\n" + "─" * 110)
    print(f"  PAPER'S BONFERRONI (k=5) FOR REFERENCE")
    print(f"  Test 1  full corpus A+:    p_raw = 0.0103, p_adj = {0.0103*5:.4f}  {'*' if 0.0103*5 < 0.05 else 'ns'}")
    print(f"  Test 2  dome/sph. raw A+:  p_raw = 0.0009, p_adj = {0.0009*5:.4f}  {'*' if 0.0009*5 < 0.05 else 'ns'}")
    print(f"  Test 2b hemiph. evo. A+:   p_raw = 0.0009, p_adj = {0.0009*5:.4f}  {'*' if 0.0009*5 < 0.05 else 'ns'}")
    print(f"  Test 3  cluster asymmetry: p_raw = 0.0033, p_adj = {0.0033*5:.4f}  {'*' if 0.0033*5 < 0.05 else 'ns'}")
    print(f"  Test 4  temporal gradient: p_raw = 0.0481, p_adj = {0.0481*5:.4f}  {'*' if 0.0481*5 < 0.05 else 'ns'}")

    # ── Interpretation ───────────────────────────────────────────────────────
    print(f"\n" + "=" * 110)
    print("=" * 110)
    n_bh_sig = sum(bh_q < 0.05)
    n_raw_sig = sum(pvals < 0.05)
    print(f"""
  {n_raw_sig} tests are significant at α = 0.05 (uncorrected).
  {n_bh_sig} survive Benjamini-Hochberg FDR correction across all {m} tests.

  The paper's Bonferroni k=5 correction covers the five confirmatory tests;
  this analysis extends to ALL {m} hypothesis tests actually performed.

  Key results surviving full FDR correction (q < 0.05):""")

    for i in range(m):
        if bh_q[i] < 0.05:
            print(f"    • {names[i]:45s}  p = {pvals[i]:.4f}  q = {bh_q[i]:.4f}")

    print()
    return 0

if __name__ == "__main__":
    main()
