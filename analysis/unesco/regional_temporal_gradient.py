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

# ── 4. Mantel-Haenszel region-adjusted test ──────────────────────────────────
print(f"\n" + "=" * 100)
print("  MANTEL-HAENSZEL REGION-ADJUSTED ANALYSIS")
print("  Tests whether the temporal gradient persists after stratifying by region")
print("=" * 100)

# For each region, compute 2×2 table: early (1978-1999) vs late (2000-2025) × A+ vs not
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

# Mantel-Haenszel common odds ratio
if mh_tables:
    MH_num = 0
    MH_den = 0
    chi_num = 0
    chi_var = 0

    for reg, a, b, c, d, n in mh_tables:
        n1 = a + b  # early total
        n2 = c + d  # late total
        m1 = a + c  # A+ total
        m2 = b + d  # non-A+ total

        MH_num += a * d / n
        MH_den += c * b / n

        # CMH test statistic components
        E_a = n1 * m1 / n
        V_a = n1 * n2 * m1 * m2 / (n**2 * (n - 1)) if n > 1 else 0
        chi_num += a - E_a
        chi_var += V_a

    MH_OR = MH_num / MH_den if MH_den > 0 else float('inf')

    # CMH chi-squared
    chi_sq = chi_num**2 / chi_var if chi_var > 0 else 0
    from scipy.stats import chi2
    p_cmh = 1 - chi2.cdf(chi_sq, df=1)

    print(f"\n  Region strata used: {len(mh_tables)}")
    print(f"  Mantel-Haenszel common OR (early vs late): {MH_OR:.3f}")
    print(f"  CMH χ² = {chi_sq:.3f}, df = 1, p = {p_cmh:.4f}  {sig(p_cmh)}")
    if MH_OR > 1:
        print(f"  → Early-inscribed sites have HIGHER A+ rate than late-inscribed,")
        print(f"    even after adjusting for regional composition changes.")
    else:
        print(f"  → No consistent early > late pattern after regional adjustment.")

    print(f"\n  Per-region summary:")
    print(f"  {'Region':<30}  {'Early A+%':>10}  {'Late A+%':>10}  {'N_early':>8}  {'N_late':>8}")
    print(f"  {'─'*75}")
    for reg, a, b, c, d, n in mh_tables:
        early_pct = 100 * a / (a + b) if (a + b) > 0 else 0
        late_pct = 100 * c / (c + d) if (c + d) > 0 else 0
        direction = "↓" if early_pct > late_pct else "↑" if early_pct < late_pct else "="
        print(f"  {reg:<30}  {early_pct:>8.1f}%  {late_pct:>8.1f}%  {a+b:>8}  {c+d:>8}  {direction}")

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

  LATEX MACROS (update tenenbaum_2026_gerizim_beru.tex):
""")

# Print LaTeX macros
if mh_tables:
    print(f"  \\newcommand{{\\MHcommonOR}}{{{MH_OR:.3f}}}          % Mantel-Haenszel common OR")
    print(f"  \\newcommand{{\\CMHchisq}}{{{chi_sq:.3f}}}           % Cochran-Mantel-Haenszel chi-sq")
    print(f"  \\newcommand{{\\pCMH}}{{{p_cmh:.4f}}}               % p-value, CMH test")
print()
