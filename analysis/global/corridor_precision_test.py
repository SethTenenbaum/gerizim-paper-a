import sys
import math
import numpy as np
from pathlib import Path
from collections import defaultdict
from scipy.stats import binomtest, fisher_exact, chi2

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APP, TIER_APLUS, TIER_B_MAX, P_NULL_AP,
    deviation as beru_deviation, tier_label, is_aplus, is_a_or_better,
)
from lib.stats import significance_label as sig

# ── Corridor definition ───────────────────────────────────────────────────────
LUMBINI_LON      = 83.2761        # UNESCO XML coordinate
LUMBINI_BERU     = abs(LUMBINI_LON - GERIZIM) / BERU   # = 1.6001
CORRIDOR_MAX     = 1.6            # nearest integer-tenth beru to Lumbini
CORRIDOR_HARMONICS = [round(n * 0.1, 1) for n in range(0, int(CORRIDOR_MAX * 10) + 1)]
N_HARMONICS      = len(CORRIDOR_HARMONICS)   # 17

SURVEY_MAX        = 5.5           # beru; covers all 1248 UNESCO sites
SURVEY_HARMONICS  = [round(n * 0.1, 1) for n in range(0, int(SURVEY_MAX * 10) + 1)]
N_SURVEY          = len(SURVEY_HARMONICS)

RUN_PERMTEST = "--permtest" in sys.argv
N_PERM       = 100_000

# ── Site loading ──────────────────────────────────────────────────────────────
corpus = load_corpus()
all_sites = cultural_sites_with_coords(corpus)   # Cultural + Mixed only

def best_site_at_harmonic(harmonic: float, sites):
    best = None
    for s in sites:
        bv = abs(s.longitude - GERIZIM) / BERU
        nearest = round(bv * 10) / 10
        if nearest != harmonic:
            continue
        dev = beru_deviation(s.longitude)
        if dev <= TIER_APLUS:
            if best is None or dev < best["dev"]:
                direction = "E" if s.longitude >= GERIZIM else "W"
                best = {
                    "site":      s.site,
                    "lon":       s.longitude,
                    "dev":       dev,
                    "km":        dev * BERU * 111,
                    "tier":      tier_label(dev),
                    "cat":       s.category,
                    "direction": direction,
                }
    return best

def all_sites_at_harmonic(harmonic: float, sites, max_tier: str = "A"):
    out = []
    for s in sites:
        bv = abs(s.longitude - GERIZIM) / BERU
        nearest = round(bv * 10) / 10
        if nearest != harmonic:
            continue
        dev = beru_deviation(s.longitude)
        tier = tier_label(dev)
        if is_aplus(tier):
            direction = "E" if s.longitude >= GERIZIM else "W"
            out.append({
                "site":      s.site,
                "lon":       s.longitude,
                "dev":       dev,
                "km":        dev * BERU * 111,
                "tier":      tier,
                "direction": direction,
            })
    return sorted(out, key=lambda x: x["dev"])

# ── Header ────────────────────────────────────────────────────────────────────
print("=" * 110)
print("  GERIZIM → LUMBINI CORRIDOR PRECISION TEST")
print(f"  Corridor: {GERIZIM}°E (0.0 beru) → {LUMBINI_LON}°E (1.6001 beru = 1.6 harmonic)")
print(f"  {N_HARMONICS} harmonics: {CORRIDOR_HARMONICS[0]:.1f} to {CORRIDOR_HARMONICS[-1]:.1f} beru")
print(f"  A+ threshold: ≤ {TIER_APLUS} beru (≤ {TIER_APLUS*BERU*111:.1f} km from harmonic)")
print(f"  A++ threshold: ≤ {TIER_APP} beru (≤ {TIER_APP*BERU*111:.0f} m from harmonic)")
print(f"  Null rate per harmonic: P₀ = {P_NULL_AP:.2f} (uniform)")
print(f"  Corpus: {len(all_sites)} UNESCO Cultural/Mixed sites with coordinates")
print("=" * 110)

# ── Corridor table ────────────────────────────────────────────────────────────
print(f"\n  {'Harm':>5}  {'Harmonic°E':>10}  {'Dir':>3}  {'Best A+ site':<50}  {'Dev':>8}  {'km':>5}  T  Additional A+ at same harmonic")
print("  " + "─" * 120)

corridor_hits    = []   # harmonics with >= 1 A+ site
corridor_misses  = []   # harmonics with no A+ site

for n in CORRIDOR_HARMONICS:
    harmonic_e_lon = GERIZIM + n * BERU   # eastward harmonic
    harmonic_w_lon = GERIZIM - n * BERU   # westward harmonic (same beru distance)
    best = best_site_at_harmonic(n, all_sites)
    all_ap_here = all_sites_at_harmonic(n, all_sites)

    if best:
        corridor_hits.append({**best, "harmonic": n, "harmonic_lon": harmonic_e_lon, "all_ap": all_ap_here})
        marker = " ◀◀" if best["tier"] == "A++" else " ◀"
        dir_str = best["direction"]
        best_harmonic_lon = harmonic_e_lon if dir_str == "E" else harmonic_w_lon
        extra = "  [+" + ", ".join(f"{s['site'][:25]}({s['direction']})" for s in all_ap_here[1:3]) + "]" if len(all_ap_here) > 1 else ""
        print(f"  {n:>5.1f}  {best_harmonic_lon:>10.3f}  {dir_str:>3}  {best['site']:<50}  "
              f"{best['dev']:>8.6f}  {best['km']:>5.1f}  {best['tier']}{marker}{extra}")
    else:
        corridor_misses.append(n)
        print(f"  {n:>5.1f}  {harmonic_e_lon:>10.3f}   E  {'—  (no A+ site)':<50}")

n_hits  = len(corridor_hits)
n_app   = sum(1 for h in corridor_hits if h["tier"] == "A++")
n_ap    = sum(1 for h in corridor_hits if h["tier"] == "A+")

print(f"\n  Corridor occupancy: {n_hits}/{N_HARMONICS} harmonics have ≥ 1 A+ site")
print(f"  Of which A++: {n_app}  |  A+: {n_ap}")

print(f"\n{'='*110}")
print("  FULL A+ SITE CATALOGUE PER CORRIDOR HARMONIC")
print(f"{'='*110}")
for h in corridor_hits:
    harmonic_lon_e = GERIZIM + h["harmonic"] * BERU
    harmonic_lon_w = GERIZIM - h["harmonic"] * BERU
    print(f"\n  ── {h['harmonic']:.1f} beru  (E-harmonic: {harmonic_lon_e:.3f}°E  |  W-harmonic: {harmonic_lon_w:.3f}°E) ──")
    for s in h["all_ap"]:
        marker = " ◀◀ A++" if s["tier"] == "A++" else " ◀ A+"
        dir_str = s["direction"]
        print(f"      {dir_str}  {s['site']:<55}  lon={s['lon']:>9.4f}  dev={s['dev']:.6f}  {s['km']:>5.1f} km{marker}")

# ── Primary statistical test ──────────────────────────────────────────────────
print(f"\n{'='*110}")
print("  STATISTICAL TESTS")
print(f"{'='*110}")

bt = binomtest(n_hits, N_HARMONICS, P_NULL_AP, alternative="greater")
p_all = P_NULL_AP ** N_HARMONICS   # exact probability of k=n=17 under H₀
print(f"""
  TEST 1: Binomial (k={n_hits} hits in n={N_HARMONICS} harmonics, p₀={P_NULL_AP})
  ──────────────────────────────────────────────────────────────────────────────────
  H₀: each harmonic independently occupied with probability P₀ = {P_NULL_AP:.3f}
  H₁: corridor occupancy rate > P₀
  
  Observed: {n_hits}/{N_HARMONICS} = {100*n_hits/N_HARMONICS:.0f}% occupied   (expected: {N_HARMONICS*P_NULL_AP:.2f})
  Binomial p (one-tailed):  {bt.pvalue:.2e}  {sig(bt.pvalue)}
  
  Exact P(all {N_HARMONICS} occupied | p={P_NULL_AP}):  {p_all:.2e}
  (This is a lower bound; the binomial p above is the correct test.)""")

log_p_sum = 0.0
chi_obs   = 0.0
for h in corridor_hits:
    p_i    = h["dev"] / TIER_B_MAX   # CDF at observed deviation
    p_i    = max(p_i, 1e-15)         # numerical floor
    log_p_sum += math.log(p_i)
    chi_obs   += -2 * math.log(p_i)

df_fisher = 2 * n_hits
from scipy.stats import chi2 as chi2_dist
p_fisher  = chi2_dist.sf(chi_obs, df_fisher)
joint_p   = math.exp(log_p_sum)

print(f"""
  TEST 2: Fisher's combined probability method
  ──────────────────────────────────────────────────────────────────────────────────
  For each occupied harmonic, P_i = dev_i / {TIER_B_MAX} (CDF under uniform null).
  Fisher statistic: X² = -2 × Σ ln(P_i) ~ χ²(2k), k = {n_hits}

  X² = {chi_obs:.2f}  (df = {df_fisher})
  Fisher combined p: {p_fisher:.2e}  {sig(p_fisher)}
  Joint probability (∏ P_i):  {joint_p:.2e}""")

corridor_harmonic_set = set(CORRIDOR_HARMONICS)
survey_non_corridor   = [n for n in SURVEY_HARMONICS if n not in corridor_harmonic_set]

n_corr_app    = n_app
n_corr_nohit  = N_HARMONICS - n_hits
n_ncorr_app   = 0
n_ncorr_total = len(survey_non_corridor)

for n in survey_non_corridor:
    b = best_site_at_harmonic(n, all_sites)
    if b and b["tier"] == "A++":
        n_ncorr_app += 1

ft_table = [
    [n_corr_app,  N_HARMONICS - n_corr_app],
    [n_ncorr_app, n_ncorr_total - n_ncorr_app],
]
ft_OR, ft_p = fisher_exact(ft_table, alternative="greater")

print(f"""
  ──────────────────────────────────────────────────────────────────────────────────
  Corridor   ({N_HARMONICS:>2} harmonics): {n_corr_app} A++  of {N_HARMONICS} harmonics
  Non-corridor ({n_ncorr_total:>2} harmonics): {n_ncorr_app} A++  of {n_ncorr_total} harmonics
  Fisher exact (one-tailed):  OR = {ft_OR:.2f},  p = {ft_p:.4f}  {sig(ft_p)}""")

# Test 4: Mean deviation comparison
corr_devs    = [h["dev"] for h in corridor_hits]
all_ap_devs  = [beru_deviation(s.longitude) for s in all_sites
                if is_aplus(tier_label(beru_deviation(s.longitude)))]
corr_mean    = sum(corr_devs) / len(corr_devs)
global_mean  = sum(all_ap_devs) / len(all_ap_devs)

from scipy.stats import mannwhitneyu
from scipy.stats import fisher_exact, norm
mw_stat, mw_p = mannwhitneyu(corr_devs, all_ap_devs, alternative="less")

print(f"""
  ──────────────────────────────────────────────────────────────────────────────────
  Corridor A+ mean deviation:    {corr_mean:.6f} beru  ({corr_mean*BERU*111:.2f} km)
  All A+ sites mean deviation:   {global_mean:.6f} beru  ({global_mean*BERU*111:.2f} km)
  Mann-Whitney U (corridor < global, one-tailed): p = {mw_p:.4f}  {sig(mw_p)}
  Corridor A+ sites are {global_mean/corr_mean:.2f}× more precise than A+ average""")

print(f"\n{'='*110}")
print("  ANCHOR COMPARISON: CORRIDOR OCCUPANCY")
print(f"{'='*110}")

anchors = [
    ("Gerizim",   GERIZIM,  "(proposed anchor)"),
    ("Jerusalem", 35.2317,  "(35.232°E)"),
    ("Megiddo",   35.1833,  "(35.183°E, Tel Megiddo)"),
    ("Bethel",    35.2250,  "(35.225°E)"),
]

def corridor_occupancy(anchor_lon, sites, n_harmonics=17, max_beru=1.6):
    harmonics = [round(n * 0.1, 1) for n in range(0, int(max_beru * 10) + 1)]
    hits = 0
    app  = 0
    for n in harmonics:
        best_dev = None
        for s in sites:
            arc = abs(s.longitude - anchor_lon)
            bv  = arc / BERU
            near = round(bv * 10) / 10
            if near != n:
                continue
            dev = abs(bv - near)
            if dev <= TIER_APLUS:
                if best_dev is None or dev < best_dev:
                    best_dev = dev
        if best_dev is not None:
            hits += 1
            if best_dev <= TIER_APP:
                app += 1
    return hits, app

print(f"\n  {'Anchor':<20}  {'Lon°E':>8}  {'Hits/17':>8}  {'A++':>5}  {'Occupancy%':>11}  {'Notes'}")
print("  " + "─" * 80)
for name, lon, notes in anchors:
    hits, app = corridor_occupancy(lon, all_sites, n_harmonics=17)
    pct = 100 * hits / 17
    bar = "█" * hits + "░" * (17 - hits)
    print(f"  {name:<20}  {lon:>8.4f}  {hits:>6}/17  {app:>5}  {pct:>10.0f}%  {bar}  {notes}")

# ── Corridor extension check ──────────────────────────────────────────────────
print(f"\n{'='*110}")
print("  CORRIDOR ROBUSTNESS: EXTENSION TO 2.0 BERU")
print(f"{'='*110}")
print(f"\n  (Additional harmonics 1.7 → 2.0 beru, not part of the primary test)")
print(f"  {'Harm':>5}  {'Lon°E':>8}  {'Best A+ site':<52}  {'Dev':>8}  {'km':>5}  T")
print("  " + "─" * 90)
for n_tenth in range(17, 21):
    n = round(n_tenth * 0.1, 1)
    harmonic_lon = GERIZIM + n * BERU
    best = best_site_at_harmonic(n, all_sites)
    if best:
        marker = " ◀◀" if best["tier"] == "A++" else " ◀"
        print(f"  {n:>5.1f}  {harmonic_lon:>8.3f}  {best['site']:<52}  "
              f"{best['dev']:>8.6f}  {best['km']:>5.1f}  {best['tier']}{marker}")
    else:
        print(f"  {n:>5.1f}  {harmonic_lon:>8.3f}  {'—':<52}")

if RUN_PERMTEST:
    print(f"\n{'='*110}")
    print(f"  PERMUTATION TEST  ({N_PERM:,} shuffles)")
    print(f"{'='*110}")
    rng = np.random.default_rng(42)
    lons = np.array([s.longitude for s in all_sites])
    perm_counts = np.zeros(N_PERM, dtype=int)
    for i in range(N_PERM):
        shuffled = rng.permutation(lons)
        count = 0
        for n in CORRIDOR_HARMONICS:
            harm_lon = GERIZIM + n * BERU
            for lon in shuffled:
                arc = abs(lon - GERIZIM)
                bv  = arc / BERU
                near = round(bv * 10) / 10
                dev  = abs(bv - near)
                if near == n and dev <= TIER_APLUS:
                    count += 1
                    break
        perm_counts[i] = count
    perm_p   = float(np.mean(perm_counts >= n_hits))
    perm_mean = float(np.mean(perm_counts))
    perm_std  = float(np.std(perm_counts))
    perm_z    = (n_hits - perm_mean) / perm_std if perm_std > 0 else float("inf")
    print(f"\n  Observed corridor occupancy: {n_hits}/17")
    print(f"  Permutation mean:            {perm_mean:.2f} ± {perm_std:.2f}")
    print(f"  Z-score:                     {perm_z:.2f}")
    print(f"  Permutation p:               {perm_p:.6f}  {sig(perm_p) if perm_p > 0 else '***'}")
else:
    print(f"\n  (Permutation test skipped — run with --permtest flag)")

# ── Verdict ───────────────────────────────────────────────────────────────────
print(f"\n{'='*110}")
print("  VERDICT")
print(f"{'='*110}")
corridor_set = set(CORRIDOR_HARMONICS)
non_corridor = [n for n in SURVEY_HARMONICS if n not in corridor_set]

n_corr_ap = n_hits                       # already computed: corridor harmonics with ≥1 A+
n_corr_noap = N_HARMONICS - n_corr_ap

n_ncorr_ap = 0
for n in non_corridor:
    b = best_site_at_harmonic(n, all_sites)
    if b:
        n_ncorr_ap += 1
n_ncorr_total = len(non_corridor)

ft_table_ap = [
    [n_corr_ap,      n_corr_noap],
    [n_ncorr_ap,     n_ncorr_total - n_ncorr_ap],
]
ft_OR_ap, ft_p_ap = fisher_exact(ft_table_ap, alternative="greater")

p1 = n_corr_ap / N_HARMONICS
p2 = n_ncorr_ap / n_ncorr_total
p_pool = (n_corr_ap + n_ncorr_ap) / (N_HARMONICS + n_ncorr_total)
se_pool = math.sqrt(p_pool * (1 - p_pool) * (1 / N_HARMONICS + 1 / n_ncorr_total)) if (N_HARMONICS > 0 and n_ncorr_total > 0) else 0.0
z_stat = (p1 - p2) / se_pool if se_pool > 0 else float("inf")
p_z_one_tailed = norm.sf(z_stat)

print(f"""
  ──────────────────────────────────────────────────────────────────────────────
  Corridor   ({N_HARMONICS:>2} harmonics): {n_corr_ap} harmonics with ≥1 A+
  Non-corridor ({n_ncorr_total:>2} harmonics): {n_ncorr_ap} harmonics with ≥1 A+
  Fisher exact (one-tailed):  OR = {ft_OR_ap:.2f},  p = {ft_p_ap:.4f}  {sig(ft_p_ap)}
  Proportion z-test (corridor > non-corridor): z = {z_stat:.2f},  p = {p_z_one_tailed:.4f}  {sig(p_z_one_tailed)}
""")

# ── LaTeX macros (GROUP 15) ───────────────────────────────────────────────────
# Jerusalem corridor occupancy
jer_hits, _ = corridor_occupancy(35.2317, all_sites, n_harmonics=17)

print("  % LaTeX macros (GROUP 15):")
print(f"  \\newcommand{{\\corridorPrecN}}{{{N_HARMONICS}}}            % harmonics tested in corridor precision")
print(f"  \\newcommand{{\\corridorPrecHits}}{{{n_hits}}}           % harmonics with corridor hit")
print(f"  \\newcommand{{\\corridorPrecApp}}{{{n_app}}}            % A++ hits in corridor precision test")
print(f"  \\newcommand{{\\corridorPrecApOnly}}{{{n_hits - n_app}}}         % A+-only hits (not A++) in corridor")
print(f"  \\newcommand{{\\pCorridorBinom}}{{{bt.pvalue:.1e}}}       % p-value, corridor binomial test")
print(f"  \\newcommand{{\\pCorridorFisher}}{{{p_fisher:.1e}}}       % p-value, corridor Fisher test")
print(f"  \\newcommand{{\\corridorAppOR}}{{{ft_OR:.2f}}}           % Fisher OR, A++ vs A+ in corridor")
print(f"  \\newcommand{{\\pCorridorAppFisher}}{{{ft_p:.4f}}}       % p-value, A++ Fisher test in corridor")
print(f"  \\newcommand{{\\corridorMeanDevKm}}{{{corr_mean * BERU * 111:.2f}}}        % mean deviation within corridor (km)")
print(f"  \\newcommand{{\\corridorGlobalMeanKm}}{{{global_mean * BERU * 111:.2f}}}      % mean deviation global (km)")
print(f"  \\newcommand{{\\corridorPrecRatio}}{{{global_mean / max(corr_mean, 1e-10):.2f}}}        % precision ratio, global/corridor deviation")
print(f"  \\newcommand{{\\pCorridorMW}}{{{mw_p:.3f}}}           % p-value, Mann-Whitney corridor precision")
print(f"  \\newcommand{{\\corridorJeruHits}}{{{jer_hits}}}          % corridor hits with Jerusalem anchor")
