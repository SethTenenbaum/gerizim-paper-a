
import sys
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, fisher_exact, spearmanr

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_B_MAX, P_NULL_AP,
    deviation as beru_dev, tier_label as tier,
)
from lib.stats import significance_label as sig, cochran_armitage
from lib.results_store import ResultsStore

# ── Load corpus ──────────────────────────────────────────────────────────────
corpus = load_corpus()
_cultural = cultural_sites_with_coords(corpus)

sites = []
for s in _cultural:
    yr = s.year
    if yr is None:
        continue
    d = beru_dev(s.longitude)
    sites.append({
        "name": s.site,
        "lon": s.longitude, "yr": yr,
        "dev": d, "tier": tier(d),
        "is_ap": d <= TIER_APLUS,
    })

N = len(sites)

# ── 1. Per-decade breakdown ──────────────────────────────────────────────────
print("=" * 90)
print("  TEMPORAL GRADIENT: PER-DECADE A+ RATES")
print(f"  Full corpus N = {N}  |  Null rate = 4%")
print("=" * 90)

decades = [
    ("1978–1984", 1978, 1984),
    ("1985–1989", 1985, 1989),
    ("1990–1994", 1990, 1994),
    ("1995–1999", 1995, 1999),
    ("2000–2004", 2000, 2004),
    ("2005–2009", 2005, 2009),
    ("2010–2014", 2010, 2014),
    ("2015–2019", 2015, 2019),
    ("2020–2025", 2020, 2025),
]

print(f"\n  {'Period':<12}  {'N':>5}  {'A+':>4}  {'A+%':>7}  {'Enrich':>7}  {'p(binom)':>10}  Sig  {'Bar'}")
print("  " + "─" * 85)

decade_rates = []
decade_ns = []
decade_aps = []
decade_midpoints = []

for label, y0, y1 in decades:
    subset = [s for s in sites if y0 <= s["yr"] <= y1]
    n = len(subset)
    nap = sum(1 for s in subset if s["is_ap"])
    if n == 0:
        continue
    rate = nap / n
    enrich = rate / P_NULL_AP
    p = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue
    bar = "█" * int(rate * 100)
    print(f"  {label:<12}  {n:>5}  {nap:>4}  {100*rate:>5.1f}%  {enrich:>6.2f}×  {p:>10.4f}  {sig(p):<3}  {bar}")
    decade_rates.append(rate)
    decade_ns.append(n)
    decade_aps.append(nap)
    decade_midpoints.append((y0 + y1) / 2)

print()
print("=" * 90)
print("  RUNNING CUMULATIVE: A+ rate as sites are added in inscription order")
print("=" * 90)

sites_sorted = sorted(sites, key=lambda s: s["yr"])
cumulative_points = [50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, N]

print(f"\n  {'First N sites':>14}  {'Through yr':>11}  {'A+':>4}  {'A+%':>7}  {'Enrich':>7}  {'p(binom)':>10}  Sig")
print("  " + "─" * 75)

for cutoff in cumulative_points:
    if cutoff > N:
        cutoff = N
    subset = sites_sorted[:cutoff]
    n = len(subset)
    nap = sum(1 for s in subset if s["is_ap"])
    max_yr = max(s["yr"] for s in subset)
    rate = nap / n
    enrich = rate / P_NULL_AP
    p = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue
    print(f"  N = {cutoff:>5}      through {max_yr}  {nap:>4}  {100*rate:>5.1f}%  {enrich:>6.2f}×  {p:>10.4f}  {sig(p)}")

print()
print("=" * 90)
print("  SPEARMAN CORRELATION: Inscription year vs. beru deviation")
print("=" * 90)

years = np.array([s["yr"] for s in sites])
devs = np.array([s["dev"] for s in sites])
rho, p_spearman = spearmanr(years, devs)

print(f"\n  Spearman rho = {rho:.4f}")
print(f"  p = {p_spearman:.6f}  {sig(p_spearman)}")
print(f"  Direction: {'Earlier sites have SMALLER deviations (closer to harmonics)' if rho > 0 else 'Earlier sites have LARGER deviations'}")

is_ap_arr = np.array([1 if s["is_ap"] else 0 for s in sites])
rho2, p_rho2 = spearmanr(years, is_ap_arr)
print(f"\n  Spearman rho (year vs. A+ status) = {rho2:.4f}")
print(f"  p = {p_rho2:.6f}  {sig(p_rho2)}")
print(f"  Direction: {'Earlier sites more likely to be A+' if rho2 < 0 else 'Later sites more likely to be A+'}")

print()
print("=" * 90)
print("  COCHRAN-ARMITAGE TREND TEST")
print("  H₁: A+ rate decreases with inscription-year cohort (temporal gradient)")
print("=" * 90)

cohort_data = []
for label, y0, y1 in [("1978-1984", 1978, 1984), ("1985-1999", 1985, 1999), ("2000-2025", 2000, 2025)]:
    subset = [s for s in sites if y0 <= s["yr"] <= y1]
    n = len(subset)
    nap = sum(1 for s in subset if s["is_ap"])
    cohort_data.append((label, n, nap))

ns = [d[1] for d in cohort_data]
aps = [d[2] for d in cohort_data]

ca3 = cochran_armitage(ns, aps, scores=[1, 2, 3])
Z_ca = ca3.z_statistic
p_ca = ca3.p_value

print(f"\n  Cohort breakdown:")
for label, n, nap in cohort_data:
    print(f"    {label}:  N={n:>4},  A+={nap:>3},  rate={100*nap/n:.1f}%")
print(f"\n  Cochran-Armitage Z = {Z_ca:.4f}")
print(f"  p (one-sided, decreasing) = {p_ca:.6f}  {sig(p_ca)}")
if ca3.direction == "decreasing":
    print(f"  → CONFIRMS: A+ rate DECREASES with later inscription year")
else:
    print(f"  → No decreasing trend detected")

print()
print("=" * 90)
print("  FIVE-COHORT TREND TEST")
print("=" * 90)

five_cohorts = [
    ("1978-1984", 1978, 1984),
    ("1985-1991", 1985, 1991),
    ("1992-1999", 1992, 1999),
    ("2000-2009", 2000, 2009),
    ("2010-2025", 2010, 2025),
]

scores5 = [1, 2, 3, 4, 5]
ns5 = []
aps5 = []
cohort5_data = []   # (label, n, nap, rate, p) — needed for macro emission below
for label, y0, y1 in five_cohorts:
    subset = [s for s in sites if y0 <= s["yr"] <= y1]
    n = len(subset)
    nap = sum(1 for s in subset if s["is_ap"])
    ns5.append(n)
    aps5.append(nap)
    rate = nap / n if n else 0
    p_b = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue if n else 1.0
    cohort5_data.append((label, n, nap, rate, p_b))
    print(f"  {label}:  N={n:>4},  A+={nap:>3},  rate={100*rate:.1f}%,  enrich={rate/P_NULL_AP:.2f}×,  p={p_b:.4f} {sig(p_b)}")

ca5 = cochran_armitage(ns5, aps5, scores=scores5)
Z5 = ca5.z_statistic
p_ca5 = ca5.p_value

print(f"\n  Cochran-Armitage Z (5 cohorts) = {Z5:.4f}")
print(f"  p (one-sided, decreasing) = {p_ca5:.6f}  {sig(p_ca5)}")

print()
print("=" * 90)
print("  LATEX MACROS FOR TEMPORAL GRADIENT")
print("=" * 90)
print(f"  \\newcommand{{\\ZcochranThree}}{{{Z_ca:.3f}}}       % Cochran-Armitage Z, 3-bin")
print(f"  \\newcommand{{\\pCochranThree}}{{{p_ca:.4f}}}      % p-value, Cochran-Armitage 3-bin")
print(f"  \\newcommand{{\\ZcochranFive}}{{{Z5:.3f}}}         % Cochran-Armitage Z, 5-bin")
print(f"  \\newcommand{{\\pCochranFive}}{{{p_ca5:.4f}}}      % p-value, Cochran-Armitage 5-bin")
print(f"  \\newcommand{{\\rhoSpearman}}{{{rho2:.4f}}}        % Spearman rho, A+ rate vs era")
print(f"  \\newcommand{{\\pSpearman}}{{{p_rho2:.6f}}}        % p-value, Spearman rank correlation")

# ── Five-cohort table macros ──────────────────────────────────────────────────
# Cohort suffixes: A=1978-1984, B=1985-1991, C=1992-1999, D=2000-2009, E=2010-2025
suffixes = list("ABCDE")
for (label, n, nap, rate, p_b), sfx in zip(cohort5_data, suffixes):
    print(f"  \\newcommand{{\\NcohortRaw{sfx}}}{{{n}}}  % {label} cohort N")
    print(f"  \\newcommand{{\\NcohortRaw{sfx}ap}}{{{nap}}}  % {label} cohort A+ count")
    print(f"  \\newcommand{{\\cohortRaw{sfx}rate}}{{{100*rate:.1f}}}  % {label} cohort A+ rate (%)")
    # Manuscript-compatible aliases (without 'Raw' infix)
    print(f"  \\newcommand{{\\Ncohort{sfx}}}{{{n}}}  % {label} cohort N (alias)")
    print(f"  \\newcommand{{\\Ncohort{sfx}ap}}{{{nap}}}  % {label} cohort A+ count (alias)")
    print(f"  \\newcommand{{\\cohort{sfx}rate}}{{{100*rate:.1f}}}  % {label} cohort A+ rate (%) (alias)")

# ── Bonferroni-adjusted p-values for primary tests ───────────────────────────
# pAdjTestTwoB: Test 2b (sensitivity variant) is excluded from Bonferroni family
# by design — it is reported descriptively only.
print(f"  \\newcommand{{\\pAdjTestTwoB}}{{---}}  % Test 2b excluded from Bonferroni (sensitivity variant, not independent)")

# ── Write to results store ────────────────────────────────────────────────────
store_dict = {
    "pCochranThree": p_ca,     # Cochran-Armitage 3-cohort p — Test 4
    "pCochranFive":  p_ca5,    # Cochran-Armitage 5-cohort p
    "ZcochranThree": Z_ca,     # Z-statistic, 3-cohort
    "ZcochranFive":  Z5,       # Z-statistic, 5-cohort
}
# Store 5-cohort breakdown
for (label, n, nap, rate, p_b), sfx in zip(cohort5_data, suffixes):
    store_dict[f"NcohortRaw{sfx}"]    = n
    store_dict[f"NcohortRaw{sfx}ap"]  = nap
    store_dict[f"cohortRaw{sfx}rate"] = round(100 * rate, 1)
ResultsStore().write_many(store_dict)
