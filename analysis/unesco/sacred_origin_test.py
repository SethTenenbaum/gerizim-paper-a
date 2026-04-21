"""
sacred_origin_test.py
=====================
Enrichment analysis: are UNESCO sites matching the "sacred_origin"
keyword category disproportionately close to beru harmonics?

DESIGN (non-cherry-picked)
──────────────────────────
1. Apply the context-aware founding_filter (lib/founding_filter.py)
   to ALL cultural/mixed UNESCO sites. This uses:
     - Unambiguous keywords (accepted on sight)
     - Ambiguous keywords with sentence-level context validation
       (require religious/sacred vocabulary in the same sentence)
   All keywords come from config.json; no site names are hardcoded.
2. Split corpus into:
     • Sacred-Origin (SO)  — sites classified as "S" by the filter
     • Remainder           — all other sites
3. Test: is the A+ rate in SO higher than in Remainder?

This is a proper 2×2 enrichment test — no site names are hardcoded,
no threshold is tuned post-hoc. The keyword list and context validation
rules are applied blinded to the beru outcome.

Tests performed:
  1. Fisher exact (SO A+ rate vs Remainder A+ rate)
  2. Binomial vs geometric null (4%)
  3. Permutation test (label-shuffling)
  4. Benjamini-Hochberg + Bonferroni correction
  5. Anchor comparison: Gerizim vs Jerusalem vs Megiddo
  6. Phase-independence check (x.18° artifact)

USAGE
-----
  python3 analysis/sacred_origin_test.py
"""

import sys
from pathlib import Path
import numpy as np
from scipy.stats import binomtest, fisher_exact, hypergeom

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APP, TIER_APLUS, TIER_B_MAX, P_NULL_AP, P_NULL_APP,
    deviation as beru_deviation, tier_label,
    CONFIG,
)
from lib.founding_filter import classify_site as _classify_site_ctx

# ── Configuration ─────────────────────────────────────────────────────────────
JERUSALEM = CONFIG["anchors"]["jerusalem"]["longitude"]
MEGIDDO   = CONFIG["anchors"]["megiddo"]["longitude"]

TIER_AP = TIER_APLUS

P0_APP = P_NULL_APP   # geometric null for A++ — from lib.beru (config-driven)
P0_AP  = P_NULL_AP    # geometric null for A+  — from lib.beru (config-driven)

N_PERM = CONFIG["simulation"]["n_permutations"]   # 100,000 from config


def _matches_so(site) -> bool:
    """True if the context-aware founding filter classifies the site as Sacred Origin (S)."""
    cats = _classify_site_ctx(site)
    return "S" in cats


def dev_from(lon, anchor):
    arc = abs(lon - anchor)
    b = arc / BERU
    return abs(b - round(b / 0.1) * 0.1)


def main():
    corpus = load_corpus()
    csites = cultural_sites_with_coords(corpus)
    N = len(csites)

    # ── Partition corpus ──────────────────────────────────────────────────────
    so_sites  = [s for s in csites if _matches_so(s)]
    rem_sites = [s for s in csites if not _matches_so(s)]

    n_so  = len(so_sites)
    n_rem = len(rem_sites)

    all_devs  = np.array([beru_deviation(s.longitude) for s in csites])
    so_devs   = np.array([beru_deviation(s.longitude) for s in so_sites])
    rem_devs  = np.array([beru_deviation(s.longitude) for s in rem_sites])

    n_ap_total  = int(np.sum(all_devs  <= TIER_AP))
    n_app_total = int(np.sum(all_devs  <= TIER_APP))
    so_ap       = int(np.sum(so_devs   <= TIER_AP))
    so_app      = int(np.sum(so_devs   <= TIER_APP))
    rem_ap      = int(np.sum(rem_devs  <= TIER_AP))
    rem_app     = int(np.sum(rem_devs  <= TIER_APP))

    emp_ap  = n_ap_total  / N
    emp_app = n_app_total / N

    print("=" * 100)
    print("  SACRED ORIGIN KEYWORD ENRICHMENT ANALYSIS")
    print(f"  Classification: context-aware filter (lib/founding_filter.py)")
    print("=" * 100)
    print(f"""
  Corpus:              N = {N}
  Sacred-Origin (SO):  n = {n_so}  ({100*n_so/N:.1f}% of corpus)
  Remainder:           n = {n_rem}

  SO  — A+ : {so_ap}/{n_so} = {100*so_ap/n_so:.1f}%
  SO  — A++: {so_app}/{n_so} = {100*so_app/n_so:.1f}%
  Rem — A+ : {rem_ap}/{n_rem} = {100*rem_ap/n_rem:.1f}%
  Rem — A++: {rem_app}/{n_rem} = {100*rem_app/n_rem:.1f}%
  All — A+ : {n_ap_total}/{N} = {100*emp_ap:.1f}%   (empirical rate)
  Null rate: {100*P0_AP:.1f}%  (geometric)
""")

    # ── Individual SO site table ──────────────────────────────────────────────
    print(f"  {'Site':<55s}  {'lon':>8s}  {'δ(beru)':>10s}  {'km':>7s}  {'Tier':>5s}")
    print("  " + "─" * 95)
    for s, d in sorted(zip(so_sites, so_devs), key=lambda x: x[1]):
        km   = d * BERU * 111.32
        tier = tier_label(d)
        mark = " ★" if d <= TIER_APP else ""
        print(f"  {s.site[:53]:53s}  {s.longitude:8.4f}  {d:10.6f}  {km:7.2f}  {tier:>5s}{mark}")

    # ── 1. Fisher exact: SO vs Remainder ────────────────────────────────────
    print("\n" + "─" * 100)
    print("  FISHER EXACT: Sacred-Origin sites vs Remainder")
    print("─" * 100)

    table_ap = [[so_ap,  n_so  - so_ap ],
                [rem_ap, n_rem - rem_ap]]
    or_ap, p_fisher_ap = fisher_exact(table_ap, alternative='greater')

    table_app = [[so_app,  n_so  - so_app ],
                 [rem_app, n_rem - rem_app]]
    or_app, p_fisher_app = fisher_exact(table_app, alternative='greater')

    def sig(p):
        return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

    print(f"""
  A+ Fisher exact (SO rate > Remainder rate):
    SO:  {so_ap}/{n_so}  |  Rem: {rem_ap}/{n_rem}
    OR = {or_ap:.2f},  p = {p_fisher_ap:.6f}  {sig(p_fisher_ap)}

  A++ Fisher exact:
    SO:  {so_app}/{n_so}  |  Rem: {rem_app}/{n_rem}
    OR = {or_app:.2f},  p = {p_fisher_app:.6f}  {sig(p_fisher_app)}
""")

    # ── 2. Binomial vs geometric null ────────────────────────────────────────
    print("─" * 100)
    print("  BINOMIAL TESTS (geometric null)")
    print("─" * 100)

    bt_ap  = binomtest(so_ap,  n_so, P0_AP,  alternative='greater')
    bt_app = binomtest(so_app, n_so, P0_APP, alternative='greater')

    print(f"""
  A+:  {so_ap}/{n_so} = {100*so_ap/n_so:.1f}%  vs null {100*P0_AP:.1f}%
    Enrichment = {(so_ap/n_so)/P0_AP:.2f}×
    p = {bt_ap.pvalue:.4f}  {sig(bt_ap.pvalue)}

  A++: {so_app}/{n_so} = {100*so_app/n_so:.1f}%  vs null {100*P0_APP:.1f}%
    Enrichment = {(so_app/n_so)/P0_APP:.2f}×
    p = {bt_app.pvalue:.4f}  {sig(bt_app.pvalue)}
""")

    # ── 3. Binomial vs empirical rate ────────────────────────────────────────
    print("─" * 100)
    print("  BINOMIAL TESTS (empirical corpus rate)")
    print("─" * 100)

    bt_ap_emp  = binomtest(so_ap,  n_so, emp_ap,  alternative='greater')
    bt_app_emp = binomtest(so_app, n_so, emp_app, alternative='greater')

    print(f"""
  A+ vs empirical {100*emp_ap:.1f}%:  p = {bt_ap_emp.pvalue:.4f}  {sig(bt_ap_emp.pvalue)}
  A++ vs empirical {100*emp_app:.2f}%:  p = {bt_app_emp.pvalue:.4f}  {sig(bt_app_emp.pvalue)}
""")

    # ── 4. Hypergeometric ────────────────────────────────────────────────────
    print("─" * 100)
    print("  HYPERGEOMETRIC TEST (A+ concentration in SO)")
    print("─" * 100)

    p_hyper_ap  = hypergeom.sf(so_ap  - 1, N, n_ap_total,  n_so)
    p_hyper_app = hypergeom.sf(so_app - 1, N, n_app_total, n_so)

    print(f"""
  A+:  P(X ≥ {so_ap} | N={N}, K={n_ap_total}, n={n_so}) = {p_hyper_ap:.6f}  {sig(p_hyper_ap)}
  A++: P(X ≥ {so_app} | N={N}, K={n_app_total}, n={n_so}) = {p_hyper_app:.6f}  {sig(p_hyper_app)}
""")

    # ── 5. Permutation test (label shuffle) ─────────────────────────────────
    print("─" * 100)
    print(f"  PERMUTATION TEST ({N_PERM:,} label shuffles)")
    print("─" * 100)

    rng = np.random.default_rng(CONFIG["simulation"]["random_seed"])
    perm_ap  = np.zeros(N_PERM, dtype=int)
    perm_app = np.zeros(N_PERM, dtype=int)
    for i in range(N_PERM):
        idx = rng.choice(N, size=n_so, replace=False)
        d   = all_devs[idx]
        perm_ap[i]  = np.sum(d <= TIER_AP)
        perm_app[i] = np.sum(d <= TIER_APP)

    p_perm_ap  = float(np.mean(perm_ap  >= so_ap))
    p_perm_app = float(np.mean(perm_app >= so_app))

    print(f"""
  A+:  P(random {n_so} sites get ≥{so_ap} A+) = {p_perm_ap:.6f}  {sig(p_perm_ap)}
  A++: P(random {n_so} sites get ≥{so_app} A++) = {p_perm_app:.6f}  {sig(p_perm_app)}
""")

    # ── 6. Anchor comparison ─────────────────────────────────────────────────
    print("─" * 100)
    print("  ANCHOR COMPARISON: GERIZIM vs JERUSALEM vs MEGIDDO")
    print("─" * 100)

    for aname, alon in [("Gerizim", GERIZIM), ("Jerusalem", JERUSALEM), ("Megiddo", MEGIDDO)]:
        devs_a = np.array([dev_from(s.longitude, alon) for s in so_sites])
        n_ap_a  = int(np.sum(devs_a <= TIER_AP))
        n_app_a = int(np.sum(devs_a <= TIER_APP))
        bt_a    = binomtest(n_ap_a,  n_so, P0_AP,  alternative='greater')
        bt_app_a = binomtest(n_app_a, n_so, P0_APP, alternative='greater')
        print(f"\n  {aname} ({alon}°E):")
        print(f"    A+ = {n_ap_a}/{n_so}  p(binom) = {bt_a.pvalue:.4f}  {sig(bt_a.pvalue)}")
        print(f"    A++ = {n_app_a}/{n_so}  p(binom) = {bt_app_a.pvalue:.4f}  {sig(bt_app_a.pvalue)}")
        for s, d in sorted(zip(so_sites, devs_a), key=lambda x: x[1])[:10]:
            km = d * BERU * 111.32
            self_mark = " [self]" if aname == "Jerusalem" and "jerusalem" in s.site.lower() else ""
            print(f"      {s.site[:52]:52s}  {d:.6f}  {km:6.2f} km  {tier_label(d)}{self_mark}")

    # ── 7. x.18° phase independence check ────────────────────────────────────
    print("\n" + "─" * 100)
    print("  x.18° PHASE INDEPENDENCE CHECK")
    print("─" * 100)
    print()
    n_near_artifact = 0
    for s, d in zip(so_sites, so_devs):
        frac = abs(s.longitude) % 1.0
        dist_18 = min(abs(frac - 0.18), abs(frac - 0.82),
                      1 - abs(frac - 0.18), 1 - abs(frac - 0.82))
        phase_flag = " ← NEAR x.18°" if dist_18 < 0.02 else ""
        if dist_18 < 0.02:
            n_near_artifact += 1
        km = d * BERU * 111.32
        print(f"  {s.site[:45]:45s}  frac={frac:.4f}  dist_from_.18={dist_18:.4f}{phase_flag}")
    print(f"\n  SO sites near x.18° artifact: {n_near_artifact}/{n_so}")

    # ── 8. Clopper-Pearson CIs ───────────────────────────────────────────────
    print("\n" + "─" * 100)
    print("  CLOPPER-PEARSON 95% CONFIDENCE INTERVALS")
    print("─" * 100)
    from scipy.stats import beta as beta_dist
    alpha = 0.05
    lo_ap = beta_dist.ppf(alpha/2, so_ap, n_so - so_ap + 1) if so_ap > 0 else 0.0
    hi_ap = beta_dist.ppf(1 - alpha/2, so_ap + 1, n_so - so_ap) if so_ap < n_so else 1.0
    print(f"\n  A+ rate: {100*so_ap/n_so:.1f}%  95% CI: [{100*lo_ap:.1f}%, {100*hi_ap:.1f}%]")
    print(f"    (Null = {100*P0_AP:.1f}%)")

    lo_app = beta_dist.ppf(alpha/2, so_app, n_so - so_app + 1) if so_app > 0 else 0.0
    hi_app = beta_dist.ppf(1 - alpha/2, so_app + 1, n_so - so_app) if so_app < n_so else 1.0
    print(f"  A++ rate: {100*so_app/n_so:.1f}%  95% CI: [{100*lo_app:.1f}%, {100*hi_app:.1f}%]")
    print(f"    (Null = {100*P0_APP:.1f}%)")

    # ── 9. Summary and multiple-comparisons correction ───────────────────────
    print("\n" + "=" * 100)
    print("  SUMMARY: ALL P-VALUES AND CORRECTIONS")
    print("=" * 100)

    all_tests = [
        ("A+ binomial (geometric null)",       bt_ap.pvalue),
        ("A+ binomial (empirical null)",        bt_ap_emp.pvalue),
        ("A++ binomial (geometric null)",       bt_app.pvalue),
        ("A++ binomial (empirical null)",       bt_app_emp.pvalue),
        ("A+ Fisher (SO vs Remainder)",         p_fisher_ap),
        ("A++ Fisher (SO vs Remainder)",        p_fisher_app),
        ("A+ hypergeometric",                   p_hyper_ap),
        ("A++ hypergeometric",                  p_hyper_app),
        ("A+ permutation",                      p_perm_ap),
        ("A++ permutation",                     p_perm_app),
    ]

    pvals = np.array([p for _, p in all_tests])
    m = len(pvals)
    bonf = np.minimum(pvals * m, 1.0)

    ranks    = np.argsort(np.argsort(pvals)) + 1
    bh_q     = np.minimum(pvals * m / ranks, 1.0)
    sort_idx = np.argsort(pvals)
    bh_adj   = np.zeros(m)
    bh_adj[sort_idx[-1]] = min(pvals[sort_idx[-1]], 1.0)
    for i in range(m - 2, -1, -1):
        bh_adj[sort_idx[i]] = min(bh_adj[sort_idx[i+1]],
                                  pvals[sort_idx[i]] * m / (i + 1))

    print(f"\n  {'Test':<45s}  {'p_raw':>10s}  {'Bonf':>10s}  {'BH q':>10s}  {'Sig':>5s}")
    print(f"  {'-'*45}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*5}")
    for i, (name, p) in enumerate(all_tests):
        s = "***" if bh_adj[i] < 0.001 else "**" if bh_adj[i] < 0.01 else "*" if bh_adj[i] < 0.05 else "ns"
        print(f"  {name:<45s}  {p:10.6f}  {bonf[i]:10.6f}  {bh_adj[i]:10.6f}  {s:>5s}")

    print(f"\n  Survive Bonferroni: {int(np.sum(bonf < 0.05))}/{m}")
    print(f"  Survive BH FDR:     {int(np.sum(bh_adj < 0.05))}/{m}")

    print("\n" + "=" * 100)
    print(f"  INTERPRETATION")
    print(f"  ─────────────────────────────────────────────────────────────────────")
    print(f"  The Sacred-Origin keyword category contains {n_so} sites ({100*n_so/N:.1f}% of corpus).")
    print(f"  Of these, {so_ap} are A+ ({100*so_ap/n_so:.1f}%) vs {100*emp_ap:.1f}% corpus rate.")
    if p_fisher_ap < 0.05:
        print(f"  The enrichment is STATISTICALLY SIGNIFICANT (Fisher p = {p_fisher_ap:.4f}).")
    else:
        print(f"  The enrichment is NOT significant (Fisher p = {p_fisher_ap:.4f}).")
        print(f"  The keyword list may be too broad, or the A+ signal is not")
        print(f"  concentrated specifically in sacred-origin sites vs. the corpus.")
    print("=" * 100)

    # ── LaTeX macros (GROUP 7) ────────────────────────────────────────────────
    # Sort sacred-origin A+ sites by deviation for top-3 naming.
    # NOTE: these macros are named \sacredHitOne/Two/Three to avoid collision
    # with \topHitOne/Two/Three, which are corridor-precision macros produced
    # by corridor_precision_test.py and anchor_site_comparison.py.
    so_ap_sites = sorted(
        [(s, beru_deviation(s.longitude)) for s in so_sites
         if beru_deviation(s.longitude) <= TIER_AP],
        key=lambda x: x[1],
    )
    top_names = [s.site for s, _ in so_ap_sites[:3]]
    print("  % LaTeX macros (GROUP 7):")
    if len(top_names) >= 1:
        print(f"  \\newcommand{{\\sacredHitOne}}{{{top_names[0]}}}  % top-ranked sacred-origin A+ site (by beru deviation)")
    if len(top_names) >= 2:
        print(f"  \\newcommand{{\\sacredHitTwo}}{{{top_names[1]}}}  % 2nd-ranked sacred-origin A+ site")
    if len(top_names) >= 3:
        print(f"  \\newcommand{{\\sacredHitThree}}{{{top_names[2]}}}  % 3rd-ranked sacred-origin A+ site")


if __name__ == "__main__":
    main()
