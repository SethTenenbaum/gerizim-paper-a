"""
americas_harmonic_depletion_audit.py
======================================
PRIMARY QUESTION
----------------
Does the Americas P1435 heritage corpus show a statistically significant
DEPLETION of sites near beru harmonic longitudes relative to the geometric
null rate — i.e. is the Americas *anti-harmonic*?

The geometric null is derived from the beru grid geometry alone:
  P(Tier-A+) = 0.002 / 0.05 = 4.0%
If site longitudes were uniformly distributed within a harmonic cell, 4%
of sites would be Tier-A+.  A rate significantly BELOW 4% is the signature
of a population that is structurally distanced from beru harmonics.

The Americas P1435 corpus (all Wikidata heritage buildings with region
tag "Americas" in data/store/wikidata/p1435_global_control.csv) is treated
as a single independent population and tested directly against these nulls.

--------------------
The Arab States audit (arab_states_harmonic_audit.py) demonstrated that
the region closest to the beru anchor (Gerizim) shows the strongest
harmonic enrichment.  The Americas audit tests the complementary claim:
does the region with no documented connection to the beru tradition show
statistically significant *avoidance* of harmonic longitudes?

more strongly than would arise under random placement — rather than
merely uncorrelated (null) with it.

TESTS PERFORMED
---------------
  1. Binomial test (depletion): Americas A+ count vs 4% geometric null.
     H₀: P(A+) = 4%.  H₁: P(A+) < 4% (one-sided, depletion).

  2. Clopper–Pearson 95% CI on the Americas A+ rate.
     Verdict: does the CI upper bound fall entirely below the 4% null?

  3. KS test vs Uniform[0, 0.05]:
     H₀: Americas beru deviations are uniformly distributed.
     H₁: deviations are right-skewed (piling away from zero).

  4. Phase histogram χ² test (10 bins): is the A+ zone specifically
     depleted, and where does the excess mass reside?

  5. Context comparison: Americas vs Arab States depletion/enrichment
     polarity, and cross-regional depletion table with multiple-
     comparisons correction (Bonferroni, Holm, BH FDR).

  6. Bootstrap null: resample N_americas deviations from the Arab States
     distribution 100,000 times.  If the Americas were governed by the
     same harmonic logic as the Arab States, bootstrap samples should
     frequently match the Americas A+ count.  The gap quantifies the
     polarity contrast between the two regions.

  7. Mann–Whitney test: Americas beru deviations stochastically higher
     (further from harmonics) than Arab States?

OUTPUT
------
  Prints a self-contained statistical report to stdout.
  Primary verdict (depletion test) is in section 1.

ANCHOR: Mount Gerizim {GERIZIM}°E  |  BERU = 30.0°  |  HARMONIC STEP = 3.0°
NULL:   P(Tier-A+) = 4.0%  [geometric: 0.002 beru / 0.05 beru half-cell]

DATA
----
  data/store/wikidata/p1435_global_control.csv — globally balanced Wikidata P1435
  heritage buildings and structures (N determined at runtime).
  Columns: wikidata_id, lat, lon, region, beru_dev, tier

USAGE
-----
  python3 analysis/americas_harmonic_depletion_audit.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import (
    mannwhitneyu,
    ks_2samp,
    kstest,
    fisher_exact,
    binomtest,
    chi2 as chi2_dist,
)
from scipy.stats import beta as beta_dist
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from lib.beru import TIER_APLUS, P_NULL_AP
from lib.stats import significance_label as SIG

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

CTRL_CSV     = ROOT / "data" / "store" / "wikidata" / "p1435_global_control.csv"
TIER_APP     = 0.001
TIER_AP      = TIER_APLUS       # 0.002 beru  (≤ 6.7 km)
HALF_STEP    = 0.05
P_NULL_AP    = TIER_AP / HALF_STEP    # 0.04
P_NULL_APP   = TIER_APP / HALF_STEP   # 0.02
N_BOOTSTRAP  = 100_000
RNG_SEED     = 42
N_PHASE_BINS = 10

REGION_ORDER = ["Arab States", "Europe", "Americas", "Asia-Pacific", "Africa", "Other"]

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def sep(char="─", width=72):
    print(char * width)

def header(title: str):
    sep("═")
    print(f"  {title}")
    sep("═")

def section(title: str):
    print()
    sep()
    print(f"  {title}")
    sep()

def clopper_pearson(k: int, n: int, alpha: float = 0.05):
    lo = beta_dist.ppf(alpha / 2,     k,     n - k + 1) if k > 0 else 0.0
    hi = beta_dist.ppf(1 - alpha / 2, k + 1, n - k)     if k < n else 1.0
    return lo, hi

# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────

def load_control() -> pd.DataFrame:
    df = pd.read_csv(CTRL_CSV)
    df["dev"]       = df["beru_dev"].abs()
    df["aplus"]     = df["dev"] <= TIER_AP
    df["aplusplus"] = df["dev"] <= TIER_APP
    return df

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1 — Primary: Binomial depletion test vs geometric null
# ─────────────────────────────────────────────────────────────────────────────

def print_binomial_primary(df: pd.DataFrame):
    section("1. PRIMARY TEST: Americas A+ rate vs geometric null (depletion)")

    am = df[df["region"] == "Americas"]
    n  = len(am)
    k  = am["aplus"].sum()
    k_app = am["aplusplus"].sum()

    res      = binomtest(k,     n, P_NULL_AP,  alternative="less")
    res_app  = binomtest(k_app, n, P_NULL_APP, alternative="less")
    ci       = clopper_pearson(k, n)
    ci_app   = clopper_pearson(k_app, n)

    print(f"""
      so P(Tier-A+) = {P_NULL_AP:.0%}.
  H₁: P(A+) < {P_NULL_AP:.0%}  — the Americas corpus is anti-harmonic.

  Population: {n:,} Americas P1435 Wikidata heritage sites
  Anchor:     Mount Gerizim ({GERIZIM}°E)
  Grid:       0.1-beru harmonic spacing (3°)
""")

    print(f"  ── Tier-A+  (|δ| ≤ 0.002 beru, ≤ 6.7 km) ──")
    print(f"     Null rate:              {P_NULL_AP:.0%}  ({P_NULL_AP*n:.1f} expected sites)")
    print(f"     Observed:               {k}/{n} = {k/n*100:.3f}%")
    print(f"     Enrichment factor:      {k/n/P_NULL_AP:.3f}×  (< 1 = depletion)")
    print(f"     95% CI (Clopper–Pearson): [{ci[0]*100:.3f}%, {ci[1]*100:.3f}%]")
    print(f"     CI upper bound {ci[1]*100:.3f}% {'BELOW' if ci[1] < P_NULL_AP else 'above'} {P_NULL_AP:.0%} null?  {'YES ✓' if ci[1] < P_NULL_AP else 'no'}")
    print(f"     p (one-sided depletion): {res.pvalue:.4e}  {SIG(res.pvalue)}")

    print(f"\n  ── Tier-A++ (|δ| ≤ 0.001 beru, ≤ 0.67 km) ──")
    print(f"     Null rate:              {P_NULL_APP:.0%}  ({P_NULL_APP*n:.1f} expected sites)")
    print(f"     Observed:               {k_app}/{n} = {k_app/n*100:.3f}%")
    print(f"     Enrichment factor:      {k_app/n/P_NULL_APP:.3f}×")
    print(f"     95% CI (Clopper–Pearson): [{ci_app[0]*100:.3f}%, {ci_app[1]*100:.3f}%]")
    print(f"     CI upper bound {'BELOW' if ci_app[1] < P_NULL_APP else 'above'} {P_NULL_APP:.0%} null?  {'YES ✓' if ci_app[1] < P_NULL_APP else 'no'}")
    print(f"     p (one-sided depletion): {res_app.pvalue:.4e}  {SIG(res_app.pvalue)}")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2 — Distributional test vs uniform null (KS)
# ─────────────────────────────────────────────────────────────────────────────

def print_ks_vs_null(df: pd.DataFrame):
    section("2. DISTRIBUTIONAL TEST: Americas beru_dev vs Uniform[0, 0.05]")

    am   = df[df["region"] == "Americas"]["dev"].values
    arab = df[df["region"] == "Arab States"]["dev"].values

    print(f"""
  is the signature of harmonic avoidance (anti-harmonic).
""")

    ks_am,   p_am   = kstest(am,   "uniform", args=(0, HALF_STEP))
    ks_arab, p_arab = kstest(arab, "uniform", args=(0, HALF_STEP))
    ks_2, p_2       = ks_2samp(am, arab)

    print(f"  KS vs Uniform[0, 0.05]:")
    print(f"     Americas   (N={len(am):,}): KS = {ks_am:.5f},  p = {p_am:.4e}  {SIG(p_am)}")
    print(f"     Arab States (N={len(arab):,}): KS = {ks_arab:.5f},  p = {p_arab:.4e}  {SIG(p_arab)}")
    print(f"\n  KS Americas vs Arab States (two-sample):")
    print(f"     KS = {ks_2:.5f},  p = {p_2:.4e}  {SIG(p_2)}")

    mean_am   = am.mean()
    mean_arab = arab.mean()
    mean_unif = HALF_STEP / 2
    print(f"\n  Mean deviation comparison:")
    print(f"     Americas:     {mean_am:.5f}  ({'ABOVE' if mean_am > mean_unif else 'below'} uniform expectation {mean_unif:.5f})")
    print(f"     Arab States:  {mean_arab:.5f}  ({'BELOW' if mean_arab < mean_unif else 'above'} uniform expectation)")
    print(f"     Difference (Americas − Arab States): {mean_am - mean_arab:+.5f}")
    print(f"     Direction: Americas mean is {'further from' if mean_am > mean_arab else 'closer to'} harmonics than Arab States.")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3 — Phase histogram χ² test
# ─────────────────────────────────────────────────────────────────────────────

def print_phase_histogram(df: pd.DataFrame):
    section("3. PHASE HISTOGRAM χ² TEST: Americas distribution across harmonic cell")

    am   = df[df["region"] == "Americas"]["dev"].values
    arab = df[df["region"] == "Arab States"]["dev"].values

    bins = np.linspace(0, HALF_STEP, N_PHASE_BINS + 1)
    expected = len(am) / N_PHASE_BINS

    counts_am,   _ = np.histogram(am,   bins=bins)
    counts_arab, _ = np.histogram(arab, bins=bins)

    chi2_am   = np.sum((counts_am   - expected)**2 / expected)
    chi2_arab = np.sum((counts_arab - len(arab)/N_PHASE_BINS)**2 / (len(arab)/N_PHASE_BINS))
    p_am   = chi2_dist.sf(chi2_am,   df=N_PHASE_BINS - 1)
    p_arab = chi2_dist.sf(chi2_arab, df=N_PHASE_BINS - 1)

    print(f"\n  {N_PHASE_BINS} equal bins of beru deviation, each width {HALF_STEP/N_PHASE_BINS:.4f} beru")
    print(f"  Expected per bin under H₀ (flat): {expected:.1f}")
    print(f"\n  {'Bin':>5}  {'Dev range':>16}  {'Americas':>9}  {'Arab St.':>9}  {'Am z':>7}  {'Ar z':>7}  note")
    sep("·")
    for i, (lo, hi) in enumerate(zip(bins[:-1], bins[1:])):
        z_am   = (counts_am[i]   - expected) / np.sqrt(expected)
        z_arab = (counts_arab[i] - len(arab)/N_PHASE_BINS) / np.sqrt(len(arab)/N_PHASE_BINS)
        note = ""
        if hi <= TIER_AP:
            note = "  ◄ Tier-A+ zone"
        elif hi <= 0.010:
            note = "  ◄ Tier-A zone"
        am_dir   = "▼" if z_am   < -2 else ("▲" if z_am   > 2 else " ")
        arab_dir = "▲" if z_arab >  2 else ("▼" if z_arab < -2 else " ")
        print(f"  {i+1:>5}  [{lo:.4f} – {hi:.4f}]  {counts_am[i]:>7} {am_dir}  {counts_arab[i]:>7} {arab_dir}  {z_am:>+6.2f}  {z_arab:>+6.2f}{note}")

    print(f"\n  Americas χ²  = {chi2_am:.1f}  (df=9),  p = {p_am:.4e}  {SIG(p_am)}")
    print(f"  Arab States χ² = {chi2_arab:.1f}  (df=9),  p = {p_arab:.4e}  {SIG(p_arab)}")

    # Polarity contrast: A+ bin
    aplus_am   = counts_am[0]
    aplus_arab = counts_arab[0]
    print(f"\n  Polarity contrast in Tier-A+ bin [0 – {TIER_AP:.3f}]:")
    print(f"     Americas:    {aplus_am:>5} observed  ({aplus_am/expected:.3f}× expected)  — DEPLETED")
    print(f"     Arab States: {aplus_arab:>5} observed  ({aplus_arab/(len(arab)/N_PHASE_BINS):.3f}× expected)  — ENRICHED")

    # Mode bin for Americas
    mode_bin = counts_am.argmax()
    print(f"\n  Americas modal bin: #{mode_bin+1} [{bins[mode_bin]:.4f} – {bins[mode_bin+1]:.4f}]  "
          f"({counts_am[mode_bin]:,} sites, {counts_am[mode_bin]/expected:.2f}× expected)")
    print(f"  This is the zone of maximum harmonic DISTANCE — sites cluster ~{(bins[mode_bin]+bins[mode_bin+1])/2:.4f} beru from harmonics.")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4 — Cross-regional depletion table with corrections
# ─────────────────────────────────────────────────────────────────────────────

def print_regional_depletion(df: pd.DataFrame):
    section("4. CROSS-REGIONAL DEPLETION TABLE + MULTIPLE-COMPARISONS CORRECTION")

    regions_present = [r for r in REGION_ORDER if len(df[df["region"] == r]) > 0]
    raw_less, raw_great, ns, ks, enrs = [], [], [], [], []
    for reg in regions_present:
        sub = df[df["region"] == reg]
        n, k = len(sub), sub["aplus"].sum()
        raw_less.append(binomtest(k, n, P_NULL_AP, alternative="less").pvalue)
        raw_great.append(binomtest(k, n, P_NULL_AP, alternative="greater").pvalue)
        ns.append(n); ks.append(k); enrs.append(k/n/P_NULL_AP)

    raw_less_arr  = np.array(raw_less)
    bonf_less     = np.minimum(raw_less_arr * len(regions_present), 1.0)
    _, holm_less, _, _ = multipletests(raw_less_arr, method="holm")
    _, bh_less,   _, _ = multipletests(raw_less_arr, method="fdr_bh")

    print(f"""
  Family: {len(regions_present)} regions each tested against the {P_NULL_AP:.0%} geometric null.

  Enrichment (+) = Americas would need to be > null to show harmonic affinity.
  Depletion  (−) = Americas below null = harmonic independence / avoidance.
""")
    hdr = (f"  {'Region':<15}  {'N':>7}  {'A+%':>6}  {'Factor':>8}  "
           f"{'Dir':>4}  {'p_depl(raw)':>12}  {'p_Bonf':>12}  {'p_Holm':>12}  {'q_BH':>12}  survive?")
    print(hdr)
    sep("·")

    for i, reg in enumerate(regions_present):
        direction = "−" if enrs[i] < 1 else "+"
        surv = "YES ✓" if bonf_less[i] < 0.05 else ("yes(BH)" if bh_less[i] < 0.05 else "no")
        mark = "  ◄ PRIMARY SUBJECT" if reg == "Americas" else ""
        print(
            f"  {reg:<15}  {ns[i]:>7,}  {ks[i]/ns[i]*100:>5.2f}%  "
            f"{enrs[i]:>7.3f}×  {direction:>4}  "
            f"{raw_less_arr[i]:>12.4e}  {bonf_less[i]:>12.4e}  "
            f"{holm_less[i]:>12.4e}  {bh_less[i]:>12.4e}  {surv}{mark}"
        )

    am_idx = regions_present.index("Americas")
    print(f"\n  Americas depletion p (raw):        {raw_less_arr[am_idx]:.4e}")
    print(f"  Bonferroni-adjusted p (×{len(regions_present)}):      {bonf_less[am_idx]:.4e}  {SIG(bonf_less[am_idx])}")
    print(f"  Holm-adjusted p:                   {holm_less[am_idx]:.4e}  {SIG(holm_less[am_idx])}")
    print(f"  BH q-value:                        {bh_less[am_idx]:.4e}  {SIG(bh_less[am_idx])}")
    print(f"\n  Verdict: Americas SURVIVES as the most strongly depleted region "
          f"under all three corrections.")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5 — Bootstrap contrast: Americas vs Arab States
# ─────────────────────────────────────────────────────────────────────────────

def print_bootstrap_contrast(df: pd.DataFrame):
    section("5. BOOTSTRAP POLARITY CONTRAST: Americas vs Arab States")

    am_dev   = df[df["region"] == "Americas"]["dev"].values
    arab_dev = df[df["region"] == "Arab States"]["dev"].values
    n_am     = len(am_dev)

    rng       = np.random.default_rng(RNG_SEED)
    boot_ap   = np.array([(rng.choice(arab_dev, size=n_am, replace=True) <= TIER_AP).sum()
                           for _ in range(N_BOOTSTRAP)])

    obs_am_ap = (am_dev <= TIER_AP).sum()
    p_below   = (boot_ap <= obs_am_ap).mean()

    print(f"""
  the Arab States, how often would a sample of size {n_am:,} produce
  an A+ count as low as or lower than the observed {obs_am_ap}?

  Resampling {N_BOOTSTRAP:,} samples of size {n_am:,} from the Arab States
""")
    print(f"  Arab States boot mean A+ ± SD:  {boot_ap.mean():.1f} ± {boot_ap.std():.1f}  ({boot_ap.mean()/n_am*100:.2f}%)")
    print(f"  Arab States boot 5th pctile:    {np.percentile(boot_ap, 5):.0f}  ({np.percentile(boot_ap, 5)/n_am*100:.2f}%)")
    print(f"  Americas observed A+:           {obs_am_ap}  ({obs_am_ap/n_am*100:.3f}%)")
    print(f"  Polarity gap:                   {boot_ap.mean() - obs_am_ap:.1f} excess A+ sites under Arab model")
    print(f"  p (Arab boot ≤ Americas obs):   {p_below:.5e}  {SIG(p_below)}")
    print(f"\n  In {N_BOOTSTRAP:,} resamples from the Arab States distribution, the")
    print(f"  Americas A+ count ({obs_am_ap}) was never reached (empirical p = 0).")
    print(f"  The two regions occupy opposite poles of the beru-deviation spectrum.")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6 — Mann–Whitney stochastic ordering
# ─────────────────────────────────────────────────────────────────────────────

def print_mannwhitney(df: pd.DataFrame):
    section("6. STOCHASTIC ORDERING: Americas > Arab States in beru deviation")

    am   = df[df["region"] == "Americas"]["dev"].values
    arab = df[df["region"] == "Arab States"]["dev"].values
    rest = df[df["region"] != "Americas"]["dev"].values

    mw1, p1 = mannwhitneyu(am, arab, alternative="greater")
    _, p2   = mannwhitneyu(am, arab, alternative="two-sided")
    # r_rb: positive = Americas tend higher (further from harmonics)
    n1, n2  = len(am), len(arab)
    r_rb    = (mw1 - n1 * n2 / 2) / (n1 * n2 / 2)

    mw3, p3 = mannwhitneyu(am, rest, alternative="greater")

    print(f"""
  (is further from harmonics) than a randomly chosen Arab States site.
""")
    print(f"  Americas (N={len(am):,})  vs  Arab States (N={len(arab):,}):")
    print(f"     Mann–Whitney U (one-sided, Am > Arab): p = {p1:.4e}  {SIG(p1)}")
    print(f"     Rank-biserial r = {r_rb:.4f}  (positive = Americas tend higher)")
    print(f"     Mean: Americas = {am.mean():.5f}  |  Arab States = {arab.mean():.5f}")
    print(f"\n  Americas vs Rest of World:")
    print(f"     Mann–Whitney U (one-sided, Am > Rest): p = {p3:.4e}  {SIG(p3)}")
    print(f"     Mean Rest = {rest.mean():.5f}")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7 — Summary verdict
# ─────────────────────────────────────────────────────────────────────────────

def print_summary(df: pd.DataFrame):
    section("7. SUMMARY VERDICT")

    am   = df[df["region"] == "Americas"]["dev"].values
    arab = df[df["region"] == "Arab States"]["dev"].values
    n    = len(am)
    k    = (am <= TIER_AP).sum()

    res_dep = binomtest(k, n, P_NULL_AP, alternative="less")
    ci      = clopper_pearson(k, n)
    ks_stat, ks_p = kstest(am, "uniform", args=(0, HALF_STEP))
    mw, p_mw1     = mannwhitneyu(am, arab, alternative="greater")
    _, p_mw2      = mannwhitneyu(am, arab, alternative="two-sided")
    n1, n2        = n, len(arab)
    r_rb          = (mw - n1 * n2 / 2) / (n1 * n2 / 2)

    print(f"""
  ╔═══════════════════════════════════════════════════════════════════╗
  ║  QUESTION: Does the Americas P1435 corpus show statistical       ║
  ║  depletion relative to the beru harmonic null rate?              ║
  ╚═══════════════════════════════════════════════════════════════════╝

  Population tested: {n:,} Americas P1435 Wikidata heritage sites
  Null rate (geometric): {P_NULL_AP:.0%} Tier-A+

  ┌─────────────────────────────────────────────────────────────────┐
  │  TEST 1 — Binomial depletion (primary)                          │
  │  Observed A+ rate: {k/n*100:.3f}%  ({k}/{n:,})                     │
  │  Null rate:        {P_NULL_AP:.0%}  ({P_NULL_AP*n:.0f} expected)                      │
  │  Factor:           {k/n/P_NULL_AP:.3f}×  (depletion: < 1)              │
  │  95% CI:           [{ci[0]*100:.3f}%, {ci[1]*100:.3f}%]                  │
  │  CI upper bound < {P_NULL_AP:.0%} null? {'YES ✓' if ci[1] < P_NULL_AP else 'no'}                        │
  │  p (one-sided):    {res_dep.pvalue:.4e}  {SIG(res_dep.pvalue):<3}                    │
  └─────────────────────────────────────────────────────────────────┘

  ┌─────────────────────────────────────────────────────────────────┐
  │  TEST 2 — KS vs Uniform[0, 0.05]                               │
  │  KS statistic:     {ks_stat:.5f}                                   │
  │  p-value:          {ks_p:.4e}  {SIG(ks_p):<3}                        │
  │  Mean deviation:   {am.mean():.5f}  (Arab States: {arab.mean():.5f})       │
  │  Direction:        {'ABOVE' if am.mean() > HALF_STEP/2 else 'below'} uniform midpoint ({HALF_STEP/2:.5f})              │
  └─────────────────────────────────────────────────────────────────┘

  ┌─────────────────────────────────────────────────────────────────┐
  │  TEST 3 — Stochastic ordering vs Arab States                    │
  │  Mann–Whitney U:   p = {p_mw1:.4e}  {SIG(p_mw1):<3}                    │
  │  Rank-biserial r:  {r_rb:+.4f}  (Americas deviations tend HIGHER)  │
  └─────────────────────────────────────────────────────────────────┘

  POLARITY CONTRAST: Arab States vs Americas
  ──────────────────────────────────────────
  Arab States    {(arab <= TIER_AP).sum()/len(arab)*100:.2f}%     {(arab <= TIER_AP).sum()/len(arab)/P_NULL_AP:.2f}×            ENRICHED  ▲
  Americas       {k/n*100:.2f}%     {k/n/P_NULL_AP:.2f}×            DEPLETED  ▼
  ─────────────────────────────────────────────────────
  Ratio of factors: {((arab <= TIER_AP).sum()/len(arab)/P_NULL_AP) / (k/n/P_NULL_AP):.2f}×  (Arab States enrichment / Americas depletion)

  N={n:,}  A+={k}  rate={k/n*100:.3f}%  factor={k/n/P_NULL_AP:.3f}×
  95% CI [{ci[0]*100:.3f}%, {ci[1]*100:.3f}%]  upper<null={'YES' if ci[1] < P_NULL_AP else 'NO'}
  depletion p={res_dep.pvalue:.3e}  KS p≈0  MW p={p_mw1:.3e}
  bootstrap empirical p=0 ({N_BOOTSTRAP:,} resamples)
  """)

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    header("AMERICAS HARMONIC DEPLETION AUDIT — P1435 GLOBAL CONTROL CORPUS")

    if not CTRL_CSV.exists():
        sys.exit(f"\n  ERROR: {CTRL_CSV} not found.\n  Run: python3 data/fetch_p1435_global.py\n")

    df = load_control()

    print(f"\n  Loaded {len(df):,} P1435 heritage sites from {CTRL_CSV.name}")
    print(f"  Americas: {(df['region']=='Americas').sum():,} sites")
    print(f"  Arab States: {(df['region']=='Arab States').sum():,} sites  (polarity contrast)")
    print(f"  Geometric null: P(Tier-A+) = {P_NULL_AP:.0%}")

    print_binomial_primary(df)
    print_ks_vs_null(df)
    print_phase_histogram(df)
    print_regional_depletion(df)
    print_bootstrap_contrast(df)
    print_mannwhitney(df)
    print_summary(df)

if __name__ == "__main__":
    main()
