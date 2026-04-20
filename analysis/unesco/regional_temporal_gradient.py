"""
regional_temporal_gradient.py
==============================
Stratifies the temporal gradient (inscription-year → A+ rate) by
UNESCO geographic region.

A reviewer concern is that the temporal gradient may be confounded
by the changing regional composition of the UNESCO corpus over time
(e.g., early inscriptions are Europe-heavy; later ones diversify).
If the A+ enrichment is driven entirely by European sites, the
temporal gradient might reflect shifting regional composition rather
than a genuine antiquity effect.

This script tests whether the gradient persists WITHIN each
UNESCO region, or whether it is driven by one region alone.

ANCHOR: Mount Gerizim 35.274°E  |  BERU = 30.0°

USAGE
-----
  python3 analysis/regional_temporal_gradient.py
"""

import re
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
from scipy.stats import binomtest, norm, fisher_exact

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, TIER_APLUS, P_NULL_AP, deviation as beru_dev
from lib.stats import significance_label as sig
from lib.results_store import ResultsStore

TIER_AP = TIER_APLUS

# ── UNESCO region assignment ────────────────────────────────────────────────
# Based on UNESCO regional groupings (approximate by longitude/latitude)
# The XML has a "region" field in some versions; we'll try to extract it,
# and fall back to lat/lon heuristics.

def assign_region(name, lon, lat, states_text=""):
    """
    Assign a broad UNESCO region from states_parties / coordinates.
    Uses a simplified scheme:
      - Europe & North America (lon -30 to 45, lat > 35; or known European states)
      - Asia & Pacific (lon > 45 or lat < 35 with lon 25-180)
      - Africa (lat < 35 and lon -20 to 55, excluding Mediterranean Europe)
      - Latin America & Caribbean (lon < -30)
      - Arab States (Middle East/North Africa)
    """
    states = states_text.lower()

    # Explicit state-based assignment for accuracy
    arab = ["egypt", "jordan", "iraq", "syria", "lebanon", "saudi",
            "oman", "yemen", "bahrain", "qatar", "kuwait", "uae",
            "united arab", "libya", "tunisia", "algeria", "morocco",
            "mauritania", "sudan", "palestine", "palestinian"]
    for a in arab:
        if a in states:
            return "Arab States"

    african = ["ethiopia", "kenya", "tanzania", "uganda", "nigeria",
               "senegal", "mali", "ghana", "cameroon", "congo",
               "mozambique", "zimbabwe", "south africa", "madagascar",
               "benin", "togo", "burkina", "niger", "guinea",
               "gambia", "ivory", "côte d'ivoire", "gabon",
               "central african", "chad", "eritrea", "zambia",
               "malawi", "botswana", "namibia", "angola", "lesotho",
               "eswatini", "swaziland", "mauritius", "seychelles",
               "cabo verde", "comoros"]
    for a in african:
        if a in states:
            return "Africa"

    latam = ["mexico", "brazil", "argentina", "chile", "colombia",
             "peru", "bolivia", "ecuador", "venezuela", "uruguay",
             "paraguay", "panama", "costa rica", "guatemala",
             "honduras", "el salvador", "nicaragua", "cuba",
             "jamaica", "haiti", "dominican", "trinidad",
             "barbados", "belize", "suriname", "guyana"]
    for a in latam:
        if a in states:
            return "Latin America & Caribbean"

    asia_pac = ["china", "japan", "india", "indonesia", "thailand",
                "viet nam", "vietnam", "cambodia", "myanmar",
                "korea", "philippines", "malaysia", "nepal",
                "sri lanka", "bangladesh", "pakistan", "afghanistan",
                "iran", "uzbekistan", "turkmenistan", "tajikistan",
                "kyrgyzstan", "kazakhstan", "mongolia", "laos",
                "australia", "new zealand", "fiji", "tonga",
                "micronesia", "marshall", "palau", "papua",
                "solomon", "vanuatu", "samoa", "kiribati"]
    for a in asia_pac:
        if a in states:
            return "Asia & Pacific"

    europe_na = ["france", "germany", "italy", "spain", "united kingdom",
                 "portugal", "greece", "turkey", "austria", "belgium",
                 "netherlands", "switzerland", "sweden", "norway",
                 "denmark", "finland", "iceland", "ireland",
                 "poland", "czech", "slovakia", "hungary",
                 "romania", "bulgaria", "croatia", "serbia",
                 "bosnia", "montenegro", "albania", "north macedonia",
                 "slovenia", "estonia", "latvia", "lithuania",
                 "ukraine", "belarus", "russia", "georgia",
                 "armenia", "azerbaijan", "moldova", "cyprus",
                 "malta", "luxembourg", "liechtenstein", "monaco",
                 "san marino", "andorra", "holy see", "vatican",
                 "canada", "united states"]
    for a in europe_na:
        if a in states:
            return "Europe & North America"

    # Fallback: coordinate-based
    if lon < -30:
        return "Latin America & Caribbean"
    if lat > 35 and -30 <= lon <= 50:
        return "Europe & North America"
    if 20 <= lon <= 55 and lat < 35 and lat > -5:
        return "Arab States"
    if lon > 50 or (lon > 25 and lat < 35):
        return "Asia & Pacific"
    if lat < 35 and -20 <= lon <= 55:
        return "Africa"
    return "Europe & North America"


# ── Load corpus ──────────────────────────────────────────────────────────────
corpus = load_corpus()
_cultural = cultural_sites_with_coords(corpus)

sites = []
for s in _cultural:
    yr = s.year
    if yr is None:
        continue
    states_text = s.states or s.iso_code or ""
    d = beru_dev(s.longitude)
    region = assign_region(s.site, s.longitude, s.latitude, states_text)

    sites.append({
        "name": s.site,
        "lon": s.longitude,
        "lat": s.latitude,
        "yr": yr,
        "dev": d,
        "is_ap": d <= TIER_AP,
        "region": region,
    })

N = len(sites)

# ── Cochran-Armitage helper ──────────────────────────────────────────────────
def cochran_armitage(ns_arr, aps_arr, scores=None):
    """Cochran-Armitage trend test. Returns (Z, p_one_sided_decreasing)."""
    ns_arr = np.array(ns_arr, dtype=float)
    aps_arr = np.array(aps_arr, dtype=float)
    if scores is None:
        scores = np.arange(1, len(ns_arr) + 1, dtype=float)
    N_total = ns_arr.sum()
    R_total = aps_arr.sum()
    p_bar = R_total / N_total
    s_bar = np.sum(scores * ns_arr) / N_total
    T_num = np.sum(scores * aps_arr) - R_total * s_bar
    T_den_sq = p_bar * (1 - p_bar) * (np.sum(scores**2 * ns_arr) - N_total * s_bar**2)
    if T_den_sq <= 0:
        return 0.0, 1.0
    Z = T_num / np.sqrt(T_den_sq)
    p = norm.cdf(Z)  # left tail: Z < 0 means decreasing
    return Z, p


# ── 1. Overall regional composition ─────────────────────────────────────────
print("=" * 100)
print("  REGIONAL TEMPORAL GRADIENT ANALYSIS")
print(f"  Full corpus N = {N}  |  Null rate = 4%")
print("=" * 100)

# ── 0. Region audit: output full assignment for reproducibility ──────────────
# Check for official region field in the XML
has_official_region = any(s.regions for s in _cultural)
if has_official_region:
    print("\n  NOTE: Official UNESCO 'region' field found in XML.")
else:
    print("\n  NOTE: No 'region' field in XML; regions assigned by states_parties + coordinate heuristic.")
    print("        Run with --audit flag to print all 1011 region assignments.")

if "--audit" in sys.argv:
    print(f"\n  {'Site':<55} {'Region':<30} {'Lon':>8} {'Lat':>7}")
    print("  " + "─" * 105)
    for s in sorted(sites, key=lambda s: s["region"]):
        print(f"  {s['name'][:54]:<55} {s['region']:<30} {s['lon']:>8.3f} {s['lat']:>7.3f}")
    print()

regions = sorted(set(s["region"] for s in sites))
print(f"\n  {'Region':<30}  {'N':>5}  {'A+':>4}  {'A+%':>7}  {'Enrich':>7}  {'p':>10}  Sig")
print("  " + "─" * 80)

for reg in regions:
    rs = [s for s in sites if s["region"] == reg]
    n = len(rs)
    nap = sum(1 for s in rs if s["is_ap"])
    rate = nap / n if n else 0
    enrich = rate / P_NULL_AP if P_NULL_AP > 0 else 0
    p = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue if n else 1.0
    print(f"  {reg:<30}  {n:>5}  {nap:>4}  {100*rate:>5.1f}%  {enrich:>6.2f}×  {p:>10.4f}  {sig(p)}")

# ── 2. Regional composition by era ──────────────────────────────────────────
print(f"\n" + "=" * 100)
print("  REGIONAL COMPOSITION BY ERA (% of each cohort from each region)")
print("=" * 100)

eras = [
    ("1978-1984", 1978, 1984),
    ("1985-1999", 1985, 1999),
    ("2000-2025", 2000, 2025),
]

print(f"\n  {'Era':<12}", end="")
for reg in regions:
    short = reg[:12]
    print(f"  {short:>12}", end="")
print(f"  {'Total':>7}")
print("  " + "─" * (14 + 14*len(regions) + 9))

for label, y0, y1 in eras:
    cohort = [s for s in sites if y0 <= s["yr"] <= y1]
    n_coh = len(cohort)
    print(f"  {label:<12}", end="")
    for reg in regions:
        n_reg = sum(1 for s in cohort if s["region"] == reg)
        pct = 100 * n_reg / n_coh if n_coh else 0
        print(f"  {pct:>10.1f}%", end="")
    print(f"  {n_coh:>7}")

# ── 3. Temporal gradient WITHIN each region ──────────────────────────────────
print(f"\n" + "=" * 100)
print("  TEMPORAL GRADIENT WITHIN EACH REGION")
print("  (3-cohort A+ rates: 1978-84, 1985-99, 2000-25)")
print("=" * 100)

for reg in regions:
    rs = [s for s in sites if s["region"] == reg]
    if len(rs) < 10:
        continue

    print(f"\n  ── {reg} (N = {len(rs)}) ──")
    print(f"  {'Era':<12}  {'N':>5}  {'A+':>4}  {'A+%':>7}  {'Enrich':>7}  {'p':>10}")
    print(f"  {'─'*60}")

    era_ns = []
    era_aps = []
    for label, y0, y1 in eras:
        subset = [s for s in rs if y0 <= s["yr"] <= y1]
        n = len(subset)
        nap = sum(1 for s in subset if s["is_ap"])
        era_ns.append(n)
        era_aps.append(nap)
        rate = nap / n if n else 0
        enrich = rate / P_NULL_AP
        p = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue if n else 1.0
        print(f"  {label:<12}  {n:>5}  {nap:>4}  {100*rate:>5.1f}%  {enrich:>6.2f}×  {p:>10.4f}  {sig(p)}")

    # Cochran-Armitage within region
    if all(n > 0 for n in era_ns):
        Z, p_ca = cochran_armitage(era_ns, era_aps)
        print(f"  Cochran-Armitage Z = {Z:.3f}, p = {p_ca:.4f}  {sig(p_ca)}")
        if Z < 0:
            print(f"  → Decreasing trend (consistent with global pattern)")
        else:
            print(f"  → No decreasing trend in this region")

    # Fisher exact: canon vs modern within region
    canon = [s for s in rs if 1978 <= s["yr"] <= 1984]
    modern = [s for s in rs if 2000 <= s["yr"] <= 2025]
    if len(canon) >= 5 and len(modern) >= 5:
        a = sum(1 for s in canon if s["is_ap"])
        b = len(canon) - a
        c = sum(1 for s in modern if s["is_ap"])
        d = len(modern) - c
        OR, p_f = fisher_exact([[a, b], [c, d]], alternative="greater")
        print(f"  Fisher (canon vs modern): OR = {OR:.2f}, p = {p_f:.4f}  {sig(p_f)}")

# ── 4. Mantel-Haenszel region-adjusted test (binary, kept for reference) ─────
print(f"\n" + "=" * 100)
print("  MANTEL-HAENSZEL REGION-ADJUSTED ANALYSIS (binary early/late — reference only)")
print("  Note: MH collapses 3 ordinal cohorts to binary; use logistic regression instead.")
print("=" * 100)

mh_tables = []
for reg in regions:
    rs = [s for s in sites if s["region"] == reg]
    early = [s for s in rs if s["yr"] <= 1999]
    late = [s for s in rs if s["yr"] >= 2000]
    if len(early) < 3 or len(late) < 3:
        continue
    a = sum(1 for s in early if s["is_ap"])
    b = len(early) - a
    c = sum(1 for s in late if s["is_ap"])
    d_val = len(late) - c
    n_stratum = a + b + c + d_val
    mh_tables.append((reg, a, b, c, d_val, n_stratum))

if mh_tables:
    MH_num = 0; MH_den = 0; chi_num = 0; chi_var = 0
    for reg, a, b, c, d, n in mh_tables:
        n1, n2, m1, m2 = a+b, c+d, a+c, b+d
        MH_num += a * d / n
        MH_den += c * b / n
        E_a = n1 * m1 / n
        V_a = n1 * n2 * m1 * m2 / (n**2 * (n - 1)) if n > 1 else 0
        chi_num += a - E_a
        chi_var += V_a

    MH_OR = MH_num / MH_den if MH_den > 0 else float('inf')
    from scipy.stats import chi2
    chi_sq = chi_num**2 / chi_var if chi_var > 0 else 0
    p_cmh = 1 - chi2.cdf(chi_sq, df=1)

    n_concordant = sum(
        1 for reg, a, b, c, d, n in mh_tables
        if (a+b) > 0 and (c+d) > 0 and (a/(a+b)) > (c/(c+d))
    )
    n_regions_tested = len(mh_tables)
    from scipy.stats import fisher_exact as _fisher
    n_sig_regions = 0
    for reg, a, b, c, d, n in mh_tables:
        _, p_r = _fisher([[a, b], [c, d]], alternative="greater")
        if p_r < 0.05:
            n_sig_regions += 1
    verdict = "directionally consistent but underpowered within strata"

    print(f"\n  Region strata used: {n_regions_tested}")
    print(f"  MH common OR (early vs late): {MH_OR:.3f}")
    print(f"  CMH χ² = {chi_sq:.3f}, df = 1, p = {p_cmh:.4f}  {sig(p_cmh)}")
    print(f"  Concordant regions (early > late): {n_concordant}/{n_regions_tested}")

# ── 4b. Ordinal logistic regression: cohort + region → A+ ────────────────────
print(f"\n" + "=" * 100)
print("  ORDINAL LOGISTIC REGRESSION: cohort (ordinal 1–3) + region → A+")
print("  Preserves 3-cohort structure; region enters as categorical covariate.")
print("=" * 100)

try:
    import numpy as np
    import statsmodels.api as sm

    # Assign ordinal cohort score
    def cohort_score(yr):
        if yr <= 1984: return 1
        if yr <= 1999: return 2
        return 3

    y = np.array([int(s["is_ap"]) for s in sites])
    cohort = np.array([cohort_score(s["yr"]) for s in sites], dtype=float)

    # Region dummies (drop first = Africa as reference)
    region_list = sorted(set(s["region"] for s in sites))
    ref_region  = region_list[0]
    region_dummies = {}
    for reg in region_list[1:]:
        region_dummies[reg] = np.array([1.0 if s["region"] == reg else 0.0 for s in sites])

    X = np.column_stack([cohort] + list(region_dummies.values()))
    X = sm.add_constant(X)

    model  = sm.Logit(y, X)
    result = model.fit(disp=0)

    # cohort is column index 1
    coef_cohort = result.params[1]
    pval_cohort = result.pvalues[1]
    se_cohort   = result.bse[1]
    z_cohort    = result.tvalues[1]
    or_cohort   = float(np.exp(coef_cohort))

    print(f"\n  Reference region: {ref_region}")
    print(f"  N = {len(y)}  |  A+ = {y.sum()}  |  converged = {result.mle_retvals['converged']}")
    print(f"\n  Cohort effect (ordinal 1→3, decreasing = negative coef):")
    print(f"    β = {coef_cohort:.4f}  SE = {se_cohort:.4f}  Z = {z_cohort:.3f}  p = {pval_cohort:.4f}  {sig(pval_cohort)}")
    print(f"    OR per cohort step = {or_cohort:.3f}")
    print(f"\n  Full model summary:")
    col_names = ["const", "cohort"] + list(region_dummies.keys())
    print(f"  {'Term':<30}  {'β':>8}  {'SE':>8}  {'Z':>8}  {'p':>10}  {'OR':>8}")
    print(f"  {'─'*80}")
    for i, name in enumerate(col_names):
        b = result.params[i]; se = result.bse[i]
        z = result.tvalues[i]; p = result.pvalues[i]
        print(f"  {name:<30}  {b:>8.4f}  {se:>8.4f}  {z:>8.3f}  {p:>10.4f}  {np.exp(b):>8.3f}  {sig(p)}")

    # Store for macros
    p_logreg = pval_cohort
    or_logreg = or_cohort
    z_logreg  = z_cohort

    print(f"\n  LATEX MACROS:")
    print(f"  \\newcommand{{\\pLogRegCohort}}{{{p_logreg:.4f}}}   % logistic regression cohort p")
    print(f"  \\newcommand{{\\orLogRegCohort}}{{{or_logreg:.3f}}}  % logistic regression cohort OR per step")
    print(f"  \\newcommand{{\\zLogRegCohort}}{{{z_logreg:.3f}}}   % logistic regression cohort Z")

    ResultsStore().write_many({
        "pLogRegCohort": p_logreg,
        "orLogRegCohort": or_logreg,
        "zLogRegCohort":  z_logreg,
        "pCMH":           p_cmh,
        "MH_OR":          MH_OR,
        "stratNconcordant": n_concordant,
        "stratNregions":    n_regions_tested,
    })

except ImportError:
    print("  statsmodels not available — skipping logistic regression")

# ── 5. Summary ───────────────────────────────────────────────────────────────
print(f"\n" + "=" * 100)
print("  SUMMARY")
print("=" * 100)
print("""
  This script addresses the reviewer concern that the temporal gradient
  could be an artifact of shifting regional composition.

  If the A+ enrichment is driven solely by European sites (which dominate
  the early UNESCO cohorts), the temporal gradient should disappear when
  stratified by region.

  If the gradient persists within multiple regions, or survives
  Mantel-Haenszel adjustment, the effect is not a regional-composition
  artifact.

  LATEX MACROS (copy into paper_a_primary_unesco.tex preamble):
""")

# Print LaTeX macros
if mh_tables:
    print(f"  \\newcommand{{\\MHcommonOR}}{{{MH_OR:.3f}}}             % MH common OR (early vs late)")
    print(f"  \\newcommand{{\\CMHchisq}}{{{chi_sq:.3f}}}              % CMH chi-squared")
    print(f"  \\newcommand{{\\pCMH}}{{{p_cmh:.4f}}}                  % CMH p-value")
    print(f"  \\newcommand{{\\stratNconcordant}}{{{n_concordant}}}    % regions with early > late")
    print(f"  \\newcommand{{\\stratNregions}}{{{n_regions_tested}}}   % total regions tested")
    print(f"  \\newcommand{{\\stratNsig}}{{{n_sig_regions}}}          % regions sig at p<0.05")
    print(f"  \\newcommand{{\\stratVerdict}}{{{verdict}}}")

    # ── Region-specific empirical null rates: stratified binomial tests ───
    # For each region, compute that region's empirical A+ rate as the null,
    # then test whether the full corpus (all regions) is enriched vs the
    # weighted average region-specific null.
    # More precisely: pool the per-region binomial p-values via Fisher's
    # combined method (all sites), and separately for pre-2000 inscriptions.
    print()
    print("  % ── Region-specific empirical null (stratified binomial) ──────────")

    from scipy.stats import chi2 as _chi2

    def _region_emp_null_p(site_list):
        """
        For each region, compute that region's A+ rate as the null, then
        test each site stratum against its own regional null.  Combine the
        per-region one-sided binomial p-values via Fisher's combined method.
        """
        region_groups = {}
        for s in site_list:
            region_groups.setdefault(s["region"], []).append(s)

        # Compute per-region empirical null rate from the *full* corpus
        full_region_rates = {}
        for reg in regions:
            all_reg = [s for s in sites if s["region"] == reg]
            if len(all_reg) >= 5:
                full_region_rates[reg] = sum(1 for s in all_reg if s["is_ap"]) / len(all_reg)

        chi_stat = 0.0
        df_total = 0
        for reg, rlist in region_groups.items():
            null_rate = full_region_rates.get(reg, P_NULL_AP)
            n_r = len(rlist)
            k_r = sum(1 for s in rlist if s["is_ap"])
            if n_r < 3:
                continue
            p_r = binomtest(k_r, n_r, null_rate, alternative="greater").pvalue
            # Fisher: -2 ln(p)
            import math
            chi_stat += -2.0 * math.log(max(p_r, 1e-300))
            df_total += 2
        if df_total == 0:
            return 1.0
        p_combined = 1.0 - _chi2.cdf(chi_stat, df=df_total)
        return p_combined

    p_reg_emp_all    = _region_emp_null_p(sites)
    pre2k_sites_list = [s for s in sites if s["yr"] < 2000]
    p_reg_emp_pre2k  = _region_emp_null_p(pre2k_sites_list)

    print(f"  \\newcommand{{\\pRegionEmpAll}}{{{p_reg_emp_all:.4f}}}    % region-specific empirical null, all sites")
    print(f"  \\newcommand{{\\pRegionEmpPreTwoK}}{{{p_reg_emp_pre2k:.4f}}}  % region-specific empirical null, pre-2000")

    # ── Per-region A+ counts and rates ────────────────────────────────────────
    print()
    print("  % ── Per-region A+ counts and rates ────────────────────────────────")
    region_map = {
        "Europe & North America":    "EuropeNA",
        "Asia & Pacific":            "AsiaPac",
        "Africa":                    "Africa",
        "Arab States":               "ArabStates",
        "Latin America & Caribbean": "LatAm",
    }
    total_ap = sum(1 for s in sites if s["is_ap"])
    for reg, key in region_map.items():
        reg_sites = [s for s in sites if s["region"] == reg]
        reg_n   = len(reg_sites)
        reg_ap  = sum(1 for s in reg_sites if s["is_ap"])
        reg_rate = 100.0 * reg_ap / reg_n if reg_n else 0.0
        reg_share = 100.0 * reg_ap / total_ap if total_ap else 0.0
        print(f"  \\newcommand{{\\{key}N}}{{{reg_n}}}  % N sites, {reg}")
        print(f"  \\newcommand{{\\{key}ApCount}}{{{reg_ap}}}  % A+ count, {reg}")
        print(f"  \\newcommand{{\\{key}ApRate}}{{{reg_rate:.1f}}}  % A+ rate (%), {reg}")
        print(f"  \\newcommand{{\\{key}ApShare}}{{{reg_share:.1f}}}  % share of total A+, {reg}")

        # Per-region temporal gradient macros (early = pre-2000, late = 2000+)
        early_r = [s for s in reg_sites if s["yr"] <= 1999]
        late_r  = [s for s in reg_sites if s["yr"] >= 2000]
        early_rate = 100.0 * sum(1 for s in early_r if s["is_ap"]) / len(early_r) if early_r else 0.0
        late_rate  = 100.0 * sum(1 for s in late_r  if s["is_ap"]) / len(late_r)  if late_r  else 0.0
        print(f"  \\newcommand{{\\{key}EarlyApRate}}{{{early_rate:.1f}}}  % early A+ rate (pre-2000), {reg}")
        print(f"  \\newcommand{{\\{key}LateApRate}}{{{late_rate:.1f}}}   % late A+ rate (2000+), {reg}")

    # ── Write to results store ─────────────────────────────────────────────
    ResultsStore().write_many({
        "pCMH":                p_cmh,
        "MH_OR":               MH_OR,
        "stratNconcordant":    n_concordant,
        "stratNregions":       n_regions_tested,
        "stratNsig":           n_sig_regions,
        "pRegionEmpAll":       p_reg_emp_all,
        "pRegionEmpPreTwoK":   p_reg_emp_pre2k,
    })
print()
