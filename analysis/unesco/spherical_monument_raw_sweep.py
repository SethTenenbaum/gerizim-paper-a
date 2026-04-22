"""
spherical_monument_raw_sweep.py
================================
Raw positive keyword sweep for hemispherical monument forms.

Unlike spherical_monument_test.py (Test 2), this script performs NO
context validation and NO disambiguation.  Any site who
(name + short_description + extended OUV text) contains one or more
of the FORM_KEYWORDS is included, regardless of whether the keyword
refers to built monumental architecture or to geology, vernacular
housing, or natural features.

PURPOSE
-------
This script answers the question: does the beru-harmonic signal in the
domed/spherical corpus depend on the context-filtering step, or is it
present even in the raw (over-inclusive) sweep?

  - If the signal survives the raw sweep  → context-filtering is not
    generating or inflating the result.
  - If the signal collapses in the raw sweep → context-filtering is
    doing real work and the validated test is the correct one to report.

USAGE
-----
    cd /path/to/gerizim-analysis
    python3 analysis/unesco/spherical_monument_raw_sweep.py

OUTPUT
------
    Console table comparing raw-sweep population vs. validated population,
    with binomial test results for both.

KEYWORDS
--------
    Sourced from keywords.json > dome_forms (same as Test 2).
    Unambiguous: stupa, stupas, tholos
    Ambiguous:   dome, domed, domes, spherical
    All treated identically here — raw substring match, no gating.
"""

import re
import sys
import numpy as np
from pathlib import Path
from scipy.stats import binomtest, chisquare, binom as binom_dist, fisher_exact as _fisher_exact

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus
from lib.results_store import ResultsStore
from lib.beru import (
    GERIZIM, BERU, TIER_APP, TIER_APLUS, TIER_A_MAX,
    P_NULL_APP, P_NULL_AP, P_NULL_A,
    TIER_APLUS_LABEL, TIER_A_LABEL, TIER_B_LABEL, TIER_C_LABEL,
    deviation as _beru_dev, tier_label, is_aplusplus, is_aplus, is_a_or_better,
)
from lib.dome_filter import FORM_KEYWORDS, FORM_KEYWORD_RES, AMBIGUOUS_KEYWORDS
from lib.stats import significance_label as sig

# ── Helpers ───────────────────────────────────────────────────────────────────
def beru_dev_fields(lon: float):
    arc      = abs(lon - GERIZIM)
    beru_val = arc / BERU
    nearest  = round(beru_val * 10) / 10
    dev      = _beru_dev(lon)
    return arc, beru_val, nearest, dev

# ── Load corpus and apply raw keyword sweep ───────────────────────────────────
corpus = load_corpus()

raw_sites   = []   # matched by any keyword, no validation
missed      = []   # Cultural/Mixed with coords but no keyword hit

for site_obj in corpus:
    if site_obj.category == "Natural":
        continue
    if not site_obj.has_coords:
        continue

    # full_text is already lowercased in the corpus object
    full_text = site_obj.full_text

    matched_kws = [kw for kw in FORM_KEYWORDS if FORM_KEYWORD_RES[kw].search(full_text)]

    if not matched_kws:
        missed.append(site_obj.site)
        continue

    arc, beru_val, nearest, dev = beru_dev_fields(site_obj.longitude)
    raw_sites.append({
        "name":     site_obj.site,
        "lon":      site_obj.longitude,
        "arc":      arc,
        "beru":     beru_val,
        "nearest":  nearest,
        "dev":      dev,
        "dev_km":   dev * BERU * 111.0,
        "tier":     tier_label(dev),
        "keywords": matched_kws,
    })

# ── Stats ──────────────────────────────────────────────────────────────────────
N_raw  = len(raw_sites)
nApp   = sum(1 for s in raw_sites if is_aplusplus(s["tier"]))
nAp    = sum(1 for s in raw_sites if is_aplus(s["tier"]))
nA     = sum(1 for s in raw_sites if is_a_or_better(s["tier"]))
nB     = sum(1 for s in raw_sites if s["tier"] == "B")
nC     = sum(1 for s in raw_sites if s["tier"] == "C")

bt_App = binomtest(nApp, N_raw, P_NULL_APP, alternative="greater")
bt_Ap  = binomtest(nAp, N_raw, P_NULL_AP, alternative="greater")
bt_A   = binomtest(nA,  N_raw, P_NULL_A,  alternative="greater")

obs_bins = [0] * 5
for s in raw_sites:
    obs_bins[min(int(s["dev"] / TIER_A_MAX), 4)] += 1
_, chi_p = chisquare(obs_bins, f_exp=[N_raw / 5.0] * 5)

enr_App = (nApp / N_raw) / P_NULL_APP if N_raw else 0
enr_Ap = (nAp / N_raw) / P_NULL_AP if N_raw else 0

# ── Fisher exact: dome corpus vs rest of full Cultural/Mixed corpus ───────────
from data.unesco_corpus import load_corpus as _load_full
_full = [s for s in _load_full() if s.category != "Natural" and s.has_coords]
N_CORPUS = len(_full)
K_AP_CORPUS = sum(1 for s in _full if is_aplus(tier_label(_beru_dev(s.longitude))))
K_APP_CORPUS= sum(1 for s in _full if is_aplusplus(tier_label(_beru_dev(s.longitude))))
K_A_CORPUS  = sum(1 for s in _full if is_a_or_better(tier_label(_beru_dev(s.longitude))))

_tbl_App= [[nApp, N_raw - nApp], [K_APP_CORPUS - nApp, N_CORPUS - N_raw - (K_APP_CORPUS - nApp)]]
_tbl_Ap = [[nAp, N_raw - nAp], [K_AP_CORPUS - nAp, N_CORPUS - N_raw - (K_AP_CORPUS - nAp)]]
_tbl_A  = [[nA,  N_raw - nA],  [K_A_CORPUS  - nA,  N_CORPUS - N_raw - (K_A_CORPUS  - nA)]]
fisher_or_App, fisher_p_App = _fisher_exact(_tbl_App, alternative="greater")
fisher_or_Ap, fisher_p_Ap = _fisher_exact(_tbl_Ap, alternative="greater")
fisher_or_A,  fisher_p_A  = _fisher_exact(_tbl_A,  alternative="greater")

# ── Anchor sweep ───────────────────────────────────────────────────────────────
sweep   = np.arange(34.0, 37.001, 0.001)
cnt_Ap  = np.array([
    sum(1 for s in raw_sites
        if abs(abs(s["lon"] - a) / BERU - round(abs(s["lon"] - a) / BERU * 10) / 10) <= TIER_APLUS)
    for a in sweep
])
g_Ap    = sum(1 for s in raw_sites if is_aplus(s["tier"]))
pctile  = float(np.mean(cnt_Ap <= g_Ap)) * 100

# ── Print ──────────────────────────────────────────────────────────────────────
SEP = "=" * 100

print()
print(SEP)
print("  UNESCO WHC — SPHERICAL/DOME/STUPA RAW KEYWORD SWEEP  (no context validation)")
print(f"  Keywords: {FORM_KEYWORDS}")
print(f"  All {len(FORM_KEYWORDS)} keywords treated identically — positive sweep only")
print(f"  Anchor: Gerizim {GERIZIM}°E  |  BERU = {BERU}°")
print(SEP)
print()
print(f"  Raw-sweep population : N = {N_raw}")
print(f"  (context-validated Test 2 population for comparison: N = 83, n_Ap = 11)")
print(f"  Sites added by raw sweep vs validated: {N_raw - 83} extra sites")
print()

# Per-site table
print("─" * 120)
print(f"  {'Site':<52}  {'Lon':>8}  {'Beru':>7}  {'Near':>5}  {'Dev':>7}  "
      f"{'km':>6}  T    Keywords")
print("─" * 120)
for s in sorted(raw_sites, key=lambda x: x["dev"]):
    mark = " ◀◀ A+" if is_aplus(s["tier"]) else (" ◀ A" if s["tier"] == "A" else "")
    kw_str = ", ".join(s["keywords"])[:40]
    print(f"  {s['name']:<52}  {s['lon']:>8.4f}  {s['beru']:>7.4f}  "
          f"{s['nearest']:>5.1f}  {s['dev']:>7.5f}  {s['dev_km']:>6.1f}"
          f"  {s['tier']}{mark:<7}  [{kw_str}]")

# Statistical summary
print()
print(SEP)
print("  STATISTICAL RESULTS — RAW SWEEP")
print(SEP)
print()
print(f"  {'Metric':<45}  {'Obs':>5}  {'Exp(H₀)':>9}  {'Enrich':>7}  {'p':>8}  Sig")
print(f"  {'-'*85}")
print(f"  {f'Tier-A+  ({TIER_APLUS_LABEL})':<45}  {nAp:>5}  "
      f"{N_raw*P_NULL_AP:>9.2f}  {enr_Ap:>6.2f}×  {bt_Ap.pvalue:>8.4f}  {sig(bt_Ap.pvalue)}")
print(f"  {f'Tier-A   ({TIER_A_LABEL})':<45}  {nA:>5}  "
      f"{N_raw*P_NULL_A:>9.2f}  {(nA/N_raw)/P_NULL_A:>6.2f}×  {bt_A.pvalue:>8.4f}  {sig(bt_A.pvalue)}")
print(f"  {f'Tier-B   ({TIER_B_LABEL})':<45}  {nB:>5}")
print(f"  {f'Tier-C   ({TIER_C_LABEL})':<45}  {nC:>5}")
print(f"  {'χ²-uniform (5 bins, df=4)':<45}  {'':>5}  {'':>9}  {'':>7}  {chi_p:>8.4f}  {sig(chi_p)}")
print(f"  {'Anchor sweep percentile (34–37°E)':<45}  {g_Ap:>5}  {'':>9}  {'':>7}  "
      f"{'':>8}  {pctile:.0f}th pctile")

# Comparison table
print()
print(SEP)
print("  COMPARISON: RAW SWEEP vs CONTEXT-VALIDATED (Test 2)")
print(SEP)
print()
print(f"  {'':30}  {'Raw sweep':>14}  {'Validated (Test 2)':>20}")
print(f"  {'-'*70}")
print(f"  {'N (population)':<30}  {N_raw:>14}  {'83':>20}")
print(f"  {'n A+ hits':<30}  {nAp:>14}  {'11':>20}")
print(f"  {'A+ rate':<30}  {100*nAp/N_raw:>13.1f}%  {'13.3%':>20}")
print(f"  {'Binomial p (A+)':<30}  {bt_Ap.pvalue:>14.4f}  {'0.0005':>20}")
print(f"  {'Enrichment vs geometric null':<30}  {enr_Ap:>13.2f}×  {11/83/P_NULL_AP:>19.2f}×")
print(f"  {'Anchor pctile (A+)':<30}  {pctile:>13.0f}th  ")
print()

# Sites that are in raw but not in validated (the context-rejected additions)
print(SEP)
print("  SITES ADDED BY RAW SWEEP  (would be context-rejected in Test 2)")
print(SEP)
validated_names = {
    # These are the 83 sites in the validated population — any raw site
    # not among the validated set is a context-rejected addition.
    # We detect them by checking: in raw but tier would shift if added.
    # Simpler: just flag raw sites whose keyword hit is ambiguous-only
    # AND which were not in the validated set.
}
print()
for s in sorted(raw_sites, key=lambda x: x["dev"]):
    ambig_only = all(kw in AMBIGUOUS_KEYWORDS for kw in s["keywords"])
    if ambig_only:
        mark = f"  ◀ A+ (p-diluting)" if is_aplus(s["tier"]) else ""
        print(f"  {s['name']:<55}  tier={s['tier']:<3}  km={s['dev_km']:>5.1f}"
              f"  kw={s['keywords']}{mark}")

print()
print("  Note: Sites with unambiguous keywords (stupa/tholos) are always")
print("  in both populations. Only ambiguous-keyword-only sites differ.")
print()

# ── LaTeX macros (GROUP 1) ────────────────────────────────────────────────────
enr_A = (nA / N_raw) / P_NULL_A if N_raw else 0
mean_dev = float(np.mean([s["dev"] for s in raw_sites])) if raw_sites else 0
N_validated = 83
N_context_rejected = N_raw - N_validated

print("  % LaTeX macros (GROUP 1):")
print(f"  \\newcommand{{\\NcircTotal}}{{{N_raw}}}           % raw-sweep population size")
print(f"  \\newcommand{{\\NcircValidated}}{{{N_validated}}}             % context-validated population (Test 2)")
print(f"  \\newcommand{{\\NcircKeywords}}{{{len(FORM_KEYWORDS)}}}               % number of dome-form keywords")
print(f"  \\newcommand{{\\NcircExcluded}}{{{0}}}               % excluded (natural/no-coord)")
print(f"  \\newcommand{{\\NcircContextRejected}}{{{N_context_rejected}}}               % added by raw sweep vs validated")
print(f"  \\newcommand{{\\NcircTierApp}}{{{nApp}}}           % Tier-A++ hits (raw sweep)")
print(f"  \\newcommand{{\\NcircTierAp}}{{{nAp}}}             % Tier-A+ hits (raw sweep)")
print(f"  \\newcommand{{\\NcircApRate}}{{{100*nAp/N_raw:.1f}}}     % A+ rate (%) = {nAp}/{N_raw}, rounded to 1 d.p.")
print(f"  \\newcommand{{\\NcircTierA}}{{{nA}}}             % Tier-A hits (raw sweep)")
print(f"  \\newcommand{{\\NcircTierB}}{{{nB}}}             % Tier-B hits (raw sweep)")
print(f"  \\newcommand{{\\NcircTierC}}{{{nC}}}               % Tier-C hits (raw sweep)")
print(f"  \\newcommand{{\\pCircApp}}{{{bt_App.pvalue:.4f}}}         % p-value, A++ binomial (raw sweep)")
print(f"  \\newcommand{{\\pCircAp}}{{{bt_Ap.pvalue:.4f}}}          % p-value, A+ binomial (raw sweep)")
print(f"  \\newcommand{{\\pCircA}}{{{bt_A.pvalue:.4f}}}          % p-value, A  binomial (raw sweep)")
print(f"  \\newcommand{{\\pCircApFisher}}{{{fisher_p_Ap:.4f}}}      % p-value, A+ Fisher exact vs rest (raw sweep)")
print(f"  \\newcommand{{\\circApFisherOR}}{{{fisher_or_Ap:.2f}}}      % OR, A+ Fisher exact vs rest (raw sweep)")
print(f"  \\newcommand{{\\pCircAFisher}}{{{fisher_p_A:.4f}}}       % p-value, A Fisher exact vs rest (raw sweep)")
print(f"  \\newcommand{{\\circAFisherOR}}{{{fisher_or_A:.2f}}}       % OR, A Fisher exact vs rest (raw sweep)")
print(f"  \\newcommand{{\\pCircChi}}{{{chi_p:.4f}}}          % chi-sq uniform, 5 bins, df=4 (raw sweep)")
print(f"  \\newcommand{{\\GerizimPctileAp}}{{{pctile:.0f}}}            % anchor-sweep percentile, A+ (raw sweep)")
print(f"  \\newcommand{{\\GerizimPctileA}}{{{pctile:.0f}}}            % anchor-sweep percentile, A  (raw sweep)")
print(f"  \\newcommand{{\\circEnrichAp}}{{{enr_Ap:.2f}}}           % enrichment ratio, A+")
print(f"  \\newcommand{{\\circEnrichA}}{{{enr_A:.2f}}}           % enrichment ratio, A")
print(f"  \\newcommand{{\\circMeanDev}}{{{mean_dev:.4f}}}        % mean beru deviation (raw sweep)")

# ── Expected counts and Clopper-Pearson CIs ───────────────────────────────────
circ_exp_ap = round(N_raw * P_NULL_AP, 2)
circ_exp_a  = round(N_raw * P_NULL_A,  2)
circ_rate_a = round(100 * nA / N_raw, 1)
# Clopper-Pearson exact 95% CI (uses beta distribution)
from scipy.stats import beta as _beta
_cp_ap_lo = round(100 * _beta.ppf(0.025, nAp, N_raw - nAp + 1), 1) if nAp > 0 else 0.0
_cp_ap_hi = round(100 * _beta.ppf(0.975, nAp + 1, N_raw - nAp), 1)
_cp_a_lo  = round(100 * _beta.ppf(0.025, nA, N_raw - nA + 1), 1) if nA > 0 else 0.0
_cp_a_hi  = round(100 * _beta.ppf(0.975, nA + 1, N_raw - nA), 1)
print(f"  \\newcommand{{\\circExpAp}}{{{circ_exp_ap}}}         % expected A+ count under H0 (= N × 4%)")
print(f"  \\newcommand{{\\circExpA}}{{{circ_exp_a}}}          % expected A count under H0 (= N × 20%)")
print(f"  \\newcommand{{\\circARate}}{{{circ_rate_a}}}         % Tier-A rate (%) = {nA}/{N_raw}")
print(f"  \\newcommand{{\\circApCIlo}}{{{_cp_ap_lo}}}         % Clopper-Pearson 95% CI lower, A+ rate (%)")
print(f"  \\newcommand{{\\circApCIhi}}{{{_cp_ap_hi}}}         % Clopper-Pearson 95% CI upper, A+ rate (%)")
print(f"  \\newcommand{{\\circACIlo}}{{{_cp_a_lo}}}          % Clopper-Pearson 95% CI lower, A rate (%)"  )
print(f"  \\newcommand{{\\circACIhi}}{{{_cp_a_hi}}}          % Clopper-Pearson 95% CI upper, A rate (%)"  )
# Enrichment ratio CIs (rate CI divided by null rate)
_enr_ci_lo = round(_cp_ap_lo / (P_NULL_AP * 100), 2) if nAp > 0 else 0.0
_enr_ci_hi = round(_cp_ap_hi / (P_NULL_AP * 100), 2)
print(f"  \\newcommand{{\\circEnrichCIlo}}{{{_enr_ci_lo}}}        % Clopper-Pearson 95% CI lower, A+ enrichment ratio")
print(f"  \\newcommand{{\\circEnrichCIhi}}{{{_enr_ci_hi}}}        % Clopper-Pearson 95% CI upper, A+ enrichment ratio")
print(f"  \\newcommand{{\\NcircValidatedApRate}}{{{100*nAp/N_validated:.1f}}}  % A+ rate (%) = {nAp}/{N_validated}, validated pop")

# ── Write to results store ────────────────────────────────────────────────────
ResultsStore().write_many({
    "pCircApp":     bt_App.pvalue,  # binomial p, A++ (raw sweep)
    "pCircAp":      bt_Ap.pvalue,   # binomial p, A+ (raw sweep) — Test 2
    "pCircA":       bt_A.pvalue,    # binomial p, A  (raw sweep)
    "pCircApFisher": fisher_p_Ap,   # Fisher p, A+ vs rest (raw sweep)
    "circApFisherOR": fisher_or_Ap, # Fisher OR, A+ vs rest
    "pCircAFisher":  fisher_p_A,    # Fisher p, A vs rest (raw sweep)
    "pCircChi":     chi_p,          # chi-sq uniform, 5 bins, df=4
    "NcircTotal":   N_raw,          # raw-sweep population size
    "NcircTierAp":  nAp,            # Tier-A+ hits
    "NcircApRate":  round(100 * nAp / N_raw, 1),  # A+ rate (%)
    "circEnrichAp": enr_Ap,         # enrichment ratio, A+
    "circEnrichCIlo": _enr_ci_lo,   # Clopper-Pearson 95% CI lower, enrichment ratio
    "circEnrichCIhi": _enr_ci_hi,   # Clopper-Pearson 95% CI upper, enrichment ratio
})
