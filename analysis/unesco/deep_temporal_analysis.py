
import re
import sys
from pathlib import Path
from scipy.stats import binomtest, fisher_exact, spearmanr, norm
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords, strip_html
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_B_MAX, P_NULL_AP,
    deviation as beru_dev, tier_label,
)
from lib.stats import significance_label as sig
from lib.results_store import ResultsStore

TIER_AP = TIER_APLUS

def tier(d):
    if d <= TIER_AP: return "A+"
    if d <= TIER_B_MAX:  return "B"
    return "C"

# ── Load corpus ──────────────────────────────────────────────────────────────
corpus = load_corpus()
_cultural = cultural_sites_with_coords(corpus)

sites = []
for s in _cultural:
    yr = s.year
    if yr is None:
        continue
    d = beru_dev(s.longitude)
    text = (s.short_description + " " + s.site).lower()
    sites.append({
        "name": s.site, "lon": s.longitude, "yr": yr,
        "dev": d, "tier": tier(d),
        "is_ap": d <= TIER_AP,
        "text": text,
        "desc": s.short_description.lower(),
        "km": d * BERU * 111.0,
    })

N = len(sites)

# ════════════════════════════════════════════════════════════════════════════════
# PART A: PRE-YEAR-0 ANTIQUITY ANALYSIS
# ════════════════════════════════════════════════════════════════════════════════
print("=" * 95)
print("  PART A: ANTIQUITY ANALYSIS — Do pre-year-0 sites snap to the grid?")
print(f"  Full corpus N = {N}  |  Null A+ rate = 4%")
print("=" * 95)

def extract_earliest_date(text):
    bce_matches = re.findall(r'(\d{1,5})\s*(?:b\.?c\.?e?\.?|bc\b)', text)
    # Pattern 2: "Xth/Xst/Xnd/Xrd century BCE/BC"
    century_bce = re.findall(r'(\d{1,2})(?:st|nd|rd|th)\s*(?:century|cent\.?)\s*(?:b\.?c\.?e?\.?|bc\b)', text)
    # Pattern 3: "Xth millennium BCE/BC"  
    millennium_bce = re.findall(r'(\d{1,2})(?:st|nd|rd|th)\s*millennium\s*(?:b\.?c\.?e?\.?|bc\b)', text)
    century_ce = re.findall(r'(\d{1,2})(?:st|nd|rd|th)\s*(?:century|cent\.?)\b(?!\s*(?:b\.?c|bce))', text)
    millennium_ce = re.findall(r'(\d{1,2})(?:st|nd|rd|th)\s*millennium\b(?!\s*(?:b\.?c|bce))', text)
    ce_dates = re.findall(r'(?:founded\s+in|built\s+in|dating\s+(?:from|to)|erected\s+in|constructed\s+in|established\s+in)\s+(\d{3,4})\b', text)
    
    dates = []
    
    for m in bce_matches:
        try:
            dates.append(-int(m))
        except: pass
    
    for m in century_bce:
        try:
            c = int(m)
            dates.append(-(c * 100 - 50))  # midpoint of century
        except: pass
    
    for m in millennium_bce:
        try:
            c = int(m)
            dates.append(-(c * 1000 - 500))  # midpoint
        except: pass
    
    for m in century_ce:
        try:
            c = int(m)
            if c <= 21:  # reasonable century
                dates.append(c * 100 - 50)  # midpoint
        except: pass
    
    for m in millennium_ce:
        try:
            c = int(m)
            if c <= 5:
                dates.append(c * 1000 - 500)
        except: pass
    
    for m in ce_dates:
        try:
            y = int(m)
            if 100 <= y <= 2025:
                dates.append(y)
        except: pass
    
    if 'neolithic' in text or 'neolithique' in text:
        dates.append(-5000)
    if 'palaeolithic' in text or 'paleolithic' in text:
        dates.append(-50000)
    if 'bronze age' in text:
        dates.append(-2000)
    if 'iron age' in text:
        dates.append(-800)
    if 'prehistoric' in text and not dates:
        dates.append(-3000)
    
    return min(dates) if dates else None

for s in sites:
    s["est_date"] = extract_earliest_date(s["text"])

n_with_date = sum(1 for s in sites if s["est_date"] is not None)
print(f"\n  Sites with extractable founding-era dates: {n_with_date}/{N} ({100*n_with_date/N:.1f}%)")

ERA_BINS = [
    ("Deep antiquity (pre-1000 BCE)",  None, -1000),
    ("Classical antiquity (-1000 to 0 CE)", -1000, 0),
    ("Pre-CE combined (before 1 CE)", None, 0),
    ("Late antiquity (1–500 CE)",      0, 500),
    ("Medieval (500–1500 CE)",         500, 1500),
    ("Early modern (1500–1800 CE)",    1500, 1800),
    ("Modern (post-1800 CE)",          1800, None),
]

print(f"\n  {'Era':<42}  {'N':>5}  {'A+':>4}  {'A+%':>7}  {'Enrich':>7}  {'p(binom)':>10}  Sig")
print("  " + "─" * 90)

era_results = {}
for label, lo, hi in ERA_BINS:
    subset = []
    for s in sites:
        if s["est_date"] is None:
            continue
        d = s["est_date"]
        if lo is not None and d < lo:
            continue
        if hi is not None and d >= hi:
            continue
        subset.append(s)
    n = len(subset)
    if n == 0:
        continue
    nap = sum(1 for s in subset if s["is_ap"])
    rate = nap / n
    enrich = rate / P_NULL_AP
    p = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue
    print(f"  {label:<42}  {n:>5}  {nap:>4}  {100*rate:>5.1f}%  {enrich:>6.2f}×  {p:>10.4f}  {sig(p)}")
    era_results[label] = {"n": n, "nap": nap, "rate": rate, "p": p, "enrich": enrich}

print()
print("  " + "─" * 90)
print("  FISHER EXACT: Pre-CE sites vs. Post-CE sites")
print("  " + "─" * 90)

dated = [s for s in sites if s["est_date"] is not None]
pre_ce = [s for s in dated if s["est_date"] < 0]
post_ce = [s for s in dated if s["est_date"] >= 0]
no_date = [s for s in sites if s["est_date"] is None]

n_pre = len(pre_ce)
nap_pre = sum(1 for s in pre_ce if s["is_ap"])
n_post = len(post_ce)
nap_post = sum(1 for s in post_ce if s["is_ap"])
n_nodate = len(no_date)
nap_nodate = sum(1 for s in no_date if s["is_ap"])

print(f"\n  Pre-CE sites:    N = {n_pre:>4},  A+ = {nap_pre:>3}  ({100*nap_pre/n_pre:.1f}%)")
print(f"  Post-CE sites:   N = {n_post:>4},  A+ = {nap_post:>3}  ({100*nap_post/n_post:.1f}%)")
print(f"  No date found:   N = {n_nodate:>4},  A+ = {nap_nodate:>3}  ({100*nap_nodate/n_nodate:.1f}%)")

table = [[nap_pre, n_pre - nap_pre],
         [nap_post, n_post - nap_post]]
or_val, p_fisher = fisher_exact(table, alternative="greater")
print(f"\n  Fisher exact (pre-CE A+ rate > post-CE A+ rate):")
print(f"    Odds ratio = {or_val:.2f}")
print(f"    p = {p_fisher:.6f}  {sig(p_fisher)}")

p_pre_binom = binomtest(nap_pre, n_pre, P_NULL_AP, alternative="greater").pvalue
print(f"\n  Binomial (pre-CE A+ rate > 4% null):")
print(f"    {nap_pre}/{n_pre} = {100*nap_pre/n_pre:.1f}%,  p = {p_pre_binom:.6f}  {sig(p_pre_binom)}")

print()
print("  " + "─" * 90)
print("  A+ SITE ANTIQUITY PROFILE")
print("  " + "─" * 90)

ap_sites = sorted([s for s in sites if s["is_ap"]], key=lambda x: x["dev"])
ap_with_date = [s for s in ap_sites if s["est_date"] is not None]
ap_pre_ce = [s for s in ap_sites if s["est_date"] is not None and s["est_date"] < 0]
ap_deep = [s for s in ap_sites if s["est_date"] is not None and s["est_date"] < -1000]

print(f"\n  Total A+ sites:                {len(ap_sites)}")
print(f"  A+ with extractable date:      {len(ap_with_date)}")
print(f"  A+ sites with pre-CE founding: {len(ap_pre_ce)} ({100*len(ap_pre_ce)/len(ap_sites):.1f}% of all A+)")
print(f"  A+ sites pre-1000 BCE:         {len(ap_deep)} ({100*len(ap_deep)/len(ap_sites):.1f}% of all A+)")

nap_total = len(ap_sites)
non_ap = [s for s in sites if not s["is_ap"]]
non_ap_pre_ce = [s for s in non_ap if s["est_date"] is not None and s["est_date"] < 0]
n_non_ap = len(non_ap)

print(f"\n  Pre-CE fraction of A+ sites:     {len(ap_pre_ce)}/{len(ap_with_date)} = {100*len(ap_pre_ce)/len(ap_with_date):.1f}% (of dated A+)")
non_ap_with_date = [s for s in non_ap if s["est_date"] is not None]
print(f"  Pre-CE fraction of non-A+ sites: {len(non_ap_pre_ce)}/{len(non_ap_with_date)} = {100*len(non_ap_pre_ce)/len(non_ap_with_date):.1f}% (of dated non-A+)")

table2 = [[len(ap_pre_ce), len(ap_with_date) - len(ap_pre_ce)],
          [len(non_ap_pre_ce), len(non_ap_with_date) - len(non_ap_pre_ce)]]
or2, p2 = fisher_exact(table2, alternative="greater")
print(f"\n  Fisher exact (A+ sites are more likely to be pre-CE):")
print(f"    Odds ratio = {or2:.2f}")
print(f"    p = {p2:.6f}  {sig(p2)}")

print()
print("  " + "─" * 90)
print("  CATALOGUE: All pre-CE Tier-A+ sites")
print("  " + "─" * 90)
for s in sorted(ap_pre_ce, key=lambda x: x["est_date"]):
    print(f"    {s['est_date']:>7} CE  |  {s['km']:>5.1f} km  |  {s['name'][:60]}")

print()
print("  " + "─" * 90)
print("  SPEARMAN: Founding date vs. beru deviation (among dated sites)")
print("  " + "─" * 90)

dated_years = np.array([s["est_date"] for s in dated])
dated_devs = np.array([s["dev"] for s in dated])
rho, p_spear = spearmanr(dated_years, dated_devs)
print(f"  Spearman rho = {rho:.4f}")
print(f"  p = {p_spear:.6f}  {sig(p_spear)}")
if rho > 0:
    print(f"  → Older sites have SMALLER beru deviations (closer to harmonics)")
else:
    print(f"  → No clear relationship")

dated_ap = np.array([1 if s["is_ap"] else 0 for s in dated])
rho2, p_sp2 = spearmanr(dated_years, dated_ap)
print(f"\n  Spearman rho (founding date vs. A+ status) = {rho2:.4f}")
print(f"  p = {p_sp2:.6f}  {sig(p_sp2)}")
if rho2 < 0:
    print(f"  → Older sites are MORE likely to be A+ ✓")

# ════════════════════════════════════════════════════════════════════════════════
# PART B: SEQUENTIAL INSCRIPTION ANALYSIS
# ════════════════════════════════════════════════════════════════════════════════
print()
print()
print("=" * 95)
print("  PART B: SEQUENTIAL INSCRIPTION ANALYSIS")
print("  How does the A+ enrichment evolve as UNESCO adds sites year by year?")
print(f"  Full corpus N = {N}  |  Null A+ rate = 4%")
print("=" * 95)

sites_sorted = sorted(sites, key=lambda s: (s["yr"], s["name"]))

# ── B1. Year-by-year cumulative ──────────────────────────────────────────────
print(f"\n  {'Through Year':>13}  {'Cum N':>6}  {'Cum A+':>7}  {'A+ Rate':>8}  {'Enrich':>7}  {'p(binom)':>10}  {'Sig':>4}  {'─log10(p)':>10}  Bar")
print("  " + "─" * 100)

yearly_data = []
cum_n = 0
cum_ap = 0
prev_year = None

unique_years = sorted(set(s["yr"] for s in sites_sorted))

peak_neglog = 0
peak_year = None
peak_n = 0
peak_ap = 0
peak_rate = 0
peak_p = 1.0

for year in unique_years:
    year_sites = [s for s in sites_sorted if s["yr"] == year]
    cum_n += len(year_sites)
    cum_ap += sum(1 for s in year_sites if s["is_ap"])
    
    if cum_n < 10:
        continue
    
    rate = cum_ap / cum_n
    enrich = rate / P_NULL_AP
    p = binomtest(cum_ap, cum_n, P_NULL_AP, alternative="greater").pvalue
    neglog = -np.log10(p) if p > 0 else 30
    
    yearly_data.append({
        "year": year, "cum_n": cum_n, "cum_ap": cum_ap,
        "rate": rate, "enrich": enrich, "p": p, "neglog": neglog,
    })
    
    if neglog > peak_neglog:
        peak_neglog = neglog
        peak_year = year
        peak_n = cum_n
        peak_ap = cum_ap
        peak_rate = rate
        peak_p = p
    
    bar = "█" * min(int(neglog * 3), 50)
    print(f"  Through {year:>4}  {cum_n:>6}  {cum_ap:>7}  {100*rate:>6.1f}%  {enrich:>6.2f}×  {p:>10.6f}  {sig(p):<4}  {neglog:>8.2f}  {bar}")

print(f"\n  ★ PEAK SIGNIFICANCE: Through {peak_year}")
print(f"    N = {peak_n}, A+ = {peak_ap}, rate = {100*peak_rate:.1f}%, p = {peak_p:.8f}  {sig(peak_p)}")
print(f"    -log10(p) = {peak_neglog:.3f}")

print()
print("=" * 95)
print("  KEY CHECKPOINTS: Significance at milestone iterations")
print("=" * 95)

checkpoints = [
    ("First 50 sites", 50),
    ("First 100 sites", 100),
    ("First 138 sites (through 1984)", 138),
    ("First 200 sites", 200),
    ("First 300 sites", 300),
    ("First 400 sites", 400),
    ("First 500 sites (through ~1999)", 500),
    ("First 600 sites", 600),
    ("First 700 sites", 700),
    ("First 800 sites", 800),
    ("First 900 sites", 900),
    ("Full corpus", N),
]

print(f"\n  {'Checkpoint':<40}  {'N':>5}  {'A+':>4}  {'Rate':>6}  {'Enrich':>7}  {'p':>10}  {'Sig':>4}  {'Year':>5}")
print("  " + "─" * 90)

for label, cutoff in checkpoints:
    if cutoff > N:
        cutoff = N
    subset = sites_sorted[:cutoff]
    n = len(subset)
    nap = sum(1 for s in subset if s["is_ap"])
    max_yr = max(s["yr"] for s in subset)
    rate = nap / n
    enrich = rate / P_NULL_AP
    p = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue
    print(f"  {label:<40}  {n:>5}  {nap:>4}  {100*rate:>4.1f}%  {enrich:>6.2f}×  {p:>10.6f}  {sig(p):<4}  {max_yr:>5}")

print()
print("=" * 95)
print("  SIGNIFICANCE TRAJECTORY: When does p < 0.05 first appear, and does it persist?")
print("=" * 95)

first_sig = None
sig_lost = None
for yd in yearly_data:
    if yd["p"] < 0.05 and first_sig is None:
        first_sig = yd
    if first_sig is not None and yd["p"] >= 0.05 and sig_lost is None:
        sig_lost = yd

if first_sig:
    print(f"\n  First significant (p < 0.05): Through {first_sig['year']}, N = {first_sig['cum_n']}, A+ = {first_sig['cum_ap']}, rate = {100*first_sig['rate']:.1f}%, p = {first_sig['p']:.6f}")
if sig_lost:
    print(f"  Significance lost:           Through {sig_lost['year']}, N = {sig_lost['cum_n']}, A+ = {sig_lost['cum_ap']}, rate = {100*sig_lost['rate']:.1f}%, p = {sig_lost['p']:.6f}")
else:
    print(f"  Significance NEVER LOST after first appearing — signal persists through all {N} sites")

# Final significance
final = yearly_data[-1]
print(f"\n  Final (full corpus): N = {final['cum_n']}, A+ = {final['cum_ap']}, rate = {100*final['rate']:.1f}%, p = {final['p']:.6f}")

print()
print("=" * 95)
print("  A+ SITES PER INSCRIPTION YEAR")
print("=" * 95)

print(f"\n  {'Year':>5}  {'New sites':>10}  {'New A+':>7}  {'A+ rate (year)':>15}  A+ site names")
print("  " + "─" * 100)

for year in unique_years:
    year_sites = [s for s in sites_sorted if s["yr"] == year]
    n_yr = len(year_sites)
    ap_yr = [s for s in year_sites if s["is_ap"]]
    nap_yr = len(ap_yr)
    if nap_yr > 0:
        names = "; ".join(f"{s['name'][:35]} ({s['km']:.1f}km)" for s in sorted(ap_yr, key=lambda x: x["dev"]))
        print(f"  {year:>5}  {n_yr:>10}  {nap_yr:>7}  {100*nap_yr/n_yr:>13.1f}%  {names[:100]}")

print()
print("=" * 95)
print("  SLIDING WINDOW: A+ rate in successive blocks of 100 sites")
print("=" * 95)

window = 100
print(f"\n  {'Block':>15}  {'Sites':>6}  {'A+':>4}  {'Rate':>6}  {'Enrich':>7}  {'p':>10}  Sig")
print("  " + "─" * 65)

for i in range(0, N, window):
    block = sites_sorted[i:i+window]
    n = len(block)
    nap = sum(1 for s in block if s["is_ap"])
    yr_lo = min(s["yr"] for s in block)
    yr_hi = max(s["yr"] for s in block)
    rate = nap / n if n else 0
    enrich = rate / P_NULL_AP if rate else 0
    p = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue if n >= 5 else 1.0
    label = f"#{i+1}–{i+n}"
    print(f"  {label:>15}  ({yr_lo}-{yr_hi})  {nap:>4}  {100*rate:>4.1f}%  {enrich:>6.2f}×  {p:>10.4f}  {sig(p)}")

print()
print("=" * 95)
print("  HALF-CORPUS COMPARISON: First half (most important) vs. Second half (less important)")
print("=" * 95)

half = N // 2
first_half = sites_sorted[:half]
second_half = sites_sorted[half:]

n1 = len(first_half)
nap1 = sum(1 for s in first_half if s["is_ap"])
n2 = len(second_half)
nap2 = sum(1 for s in second_half if s["is_ap"])
yr1_max = max(s["yr"] for s in first_half)

p1 = binomtest(nap1, n1, P_NULL_AP, alternative="greater").pvalue
p2 = binomtest(nap2, n2, P_NULL_AP, alternative="greater").pvalue

print(f"\n  First half (sites 1–{half}, through {yr1_max}):")
print(f"    N = {n1}, A+ = {nap1}, rate = {100*nap1/n1:.1f}%, enrich = {nap1/n1/P_NULL_AP:.2f}×, p = {p1:.6f}  {sig(p1)}")
print(f"\n  Second half (sites {half+1}–{N}):")
print(f"    N = {n2}, A+ = {nap2}, rate = {100*nap2/n2:.1f}%, enrich = {nap2/n2/P_NULL_AP:.2f}×, p = {p2:.6f}  {sig(p2)}")

table3 = [[nap1, n1 - nap1], [nap2, n2 - nap2]]
or3, p3 = fisher_exact(table3, alternative="greater")
print(f"\n  Fisher exact (first half A+ > second half A+):")
print(f"    Odds ratio = {or3:.2f}, p = {p3:.6f}  {sig(p3)}")

print()
print("=" * 95)
print("  MONOTONIC DILUTION: -log10(p) as corpus grows")
print("  (Higher = more significant; the PEAK tells us when the purest signal existed)")
print("=" * 95)

sorted_by_sig = sorted(yearly_data, key=lambda x: -x["neglog"])
print(f"\n  TOP 10 most significant corpus snapshots:")
print(f"  {'Rank':>5}  {'Through Year':>13}  {'N':>5}  {'A+':>4}  {'Rate':>6}  {'p':>12}  {'-log10(p)':>10}")
print("  " + "─" * 65)
for i, yd in enumerate(sorted_by_sig[:10]):
    print(f"  {i+1:>5}  {yd['year']:>13}  {yd['cum_n']:>5}  {yd['cum_ap']:>4}  {100*yd['rate']:>4.1f}%  {yd['p']:>12.8f}  {yd['neglog']:>10.3f}")

# ════════════════════════════════════════════════════════════════════════════════
# ════════════════════════════════════════════════════════════════════════════════
print()
print()
print("=" * 95)
print("  PART C: INTERACTION — Pre-CE founding date × Early inscription year")
print("  The DOUBLE FILTER: sites that are both ancient AND early-inscribed")
print("=" * 95)

canon_pre_ce = [s for s in sites if s["yr"] <= 1984 and s["est_date"] is not None and s["est_date"] < 0]
canon_post_ce = [s for s in sites if s["yr"] <= 1984 and s["est_date"] is not None and s["est_date"] >= 0]
modern_pre_ce = [s for s in sites if s["yr"] >= 2000 and s["est_date"] is not None and s["est_date"] < 0]
modern_post_ce = [s for s in sites if s["yr"] >= 2000 and s["est_date"] is not None and s["est_date"] >= 0]

groups = [
    ("Canon (≤1984) + Pre-CE", canon_pre_ce),
    ("Canon (≤1984) + Post-CE", canon_post_ce),
    ("Modern (≥2000) + Pre-CE", modern_pre_ce),
    ("Modern (≥2000) + Post-CE", modern_post_ce),
]

print(f"\n  {'Group':<35}  {'N':>5}  {'A+':>4}  {'Rate':>6}  {'Enrich':>7}  {'p':>10}  Sig")
print("  " + "─" * 80)

for label, subset in groups:
    n = len(subset)
    if n == 0:
        print(f"  {label:<35}  {0:>5}  {0:>4}  {'N/A':>6}")
        continue
    nap = sum(1 for s in subset if s["is_ap"])
    rate = nap / n
    enrich = rate / P_NULL_AP
    p = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue
    print(f"  {label:<35}  {n:>5}  {nap:>4}  {100*rate:>4.1f}%  {enrich:>6.2f}×  {p:>10.6f}  {sig(p)}")

# ════════════════════════════════════════════════════════════════════════════════
# LATEX MACROS
# ════════════════════════════════════════════════════════════════════════════════
print()
print("=" * 95)
print("  LATEX MACROS FOR THE MANUSCRIPT")
print("=" * 95)

# Pre-CE stats
pre_ce_rate = nap_pre / n_pre if n_pre else 0
pre_ce_enrich = pre_ce_rate / P_NULL_AP
print(f"\n  % Pre-year-0 (antiquity) analysis")
print(f"  \\newcommand{{\\NdatedSites}}{{{n_with_date}}}          % sites with known founding date")
print(f"  \\newcommand{{\\NpreCE}}{{{n_pre}}}             % pre-CE (antiquity) sites")
print(f"  \\newcommand{{\\NpreCEap}}{{{nap_pre}}}            % A+ sites pre-CE")
print(f"  \\newcommand{{\\preCErate}}{{{100*pre_ce_rate:.1f}}}           % A+ rate in pre-CE sites (%)")
print(f"  \\newcommand{{\\preCEenrich}}{{{pre_ce_enrich:.2f}}}          % enrichment, pre-CE A+")
print(f"  \\newcommand{{\\pPreCEbinom}}{{{p_pre_binom:.4f}}}        % p-value, pre-CE A+ binomial")
print(f"  \\newcommand{{\\NpostCE}}{{{n_post}}}            % post-CE sites")
print(f"  \\newcommand{{\\NpostCEap}}{{{nap_post}}}           % A+ sites post-CE")
print(f"  \\newcommand{{\\postCErate}}{{{100*nap_post/n_post:.1f}}}           % A+ rate in post-CE sites (%)")
print(f"  \\newcommand{{\\preCEfisherOR}}{{{or_val:.2f}}}          % Fisher OR, pre- vs post-CE A+")
print(f"  \\newcommand{{\\pPreCEfisher}}{{{p_fisher:.4f}}}         % p-value, Fisher pre- vs post-CE")

print(f"\n  % A+ antiquity profile")
print(f"  \\newcommand{{\\NapPreCE}}{{{len(ap_pre_ce)}}}           % A+ sites with pre-CE founding date")
print(f"  \\newcommand{{\\apPreCEpct}}{{{100*len(ap_pre_ce)/len(ap_sites):.0f}}}            % % of all A+ that are pre-CE")
print(f"  \\newcommand{{\\NapDeep}}{{{len(ap_deep)}}}             % A+ sites deeply pre-CE (antiquity)")
print(f"  \\newcommand{{\\apPreCEofDated}}{{{100*len(ap_pre_ce)/len(ap_with_date):.0f}}}           % % of dated A+ that are pre-CE")
frac_non_ap_pre = len(non_ap_pre_ce)/len(non_ap_with_date) if len(non_ap_with_date) > 0 else 0
print(f"  \\newcommand{{\\nonApPreCEpct}}{{{100*frac_non_ap_pre:.0f}}}           % % of dated non-A+ that are pre-CE")
print(f"  \\newcommand{{\\preCEantiqFisherOR}}{{{or2:.2f}}}       % Fisher OR, A+ vs non-A+ pre-CE antiquity")
print(f"  \\newcommand{{\\pPreCEantiqFisher}}{{{p2:.4f}}}         % p-value, A+ vs non-A+ antiquity Fisher")

print(f"\n  % Spearman: founding date vs A+ status")
print(f"  \\newcommand{{\\rhoFoundDate}}{{{rho2:.4f}}}          % Spearman rho, founding date vs A+ status")
print(f"  \\newcommand{{\\pFoundDateSpearman}}{{{p_sp2:.4f}}}      % p-value, founding date Spearman")

# Sequential/peak analysis
print(f"\n  % Sequential inscription analysis")
print(f"  \\newcommand{{\\peakSigYear}}{{{peak_year}}}             % UNESCO inscription year with peak A+ signal")
print(f"  \\newcommand{{\\peakSigN}}{{{peak_n}}}              % cumulative sites at peak signal year")
print(f"  \\newcommand{{\\peakSigAp}}{{{peak_ap}}}             % cumulative A+ at peak signal year")
print(f"  \\newcommand{{\\peakSigRate}}{{{100*peak_rate:.1f}}}           % A+ rate at peak signal year (%)")
print(f"  \\newcommand{{\\peakSigP}}{{{peak_p:.8f}}}     % p-value at peak signal year")
print(f"  \\newcommand{{\\peakNegLogP}}{{{peak_neglog:.2f}}}          % -log10(p) at peak signal year")

if first_sig:
    print(f"  \\newcommand{{\\firstSigYear}}{{{first_sig['year']}}}           % first year A+ signal reached p<0.05")
    print(f"  \\newcommand{{\\firstSigN}}{{{first_sig['cum_n']}}}             % cumulative sites at first significant year")

print(f"\n  % Half-corpus comparison")
print(f"  \\newcommand{{\\NfirstHalf}}{{{n1}}}             % sites inscribed in first half of UNESCO record")
print(f"  \\newcommand{{\\NapFirstHalf}}{{{nap1}}}            % A+ sites in first half")
print(f"  \\newcommand{{\\firstHalfRate}}{{{100*nap1/n1:.1f}}}           % A+ rate in first half (%)")
print(f"  \\newcommand{{\\pFirstHalf}}{{{p1:.6f}}}       % p-value, first-half A+ binomial")
print(f"  \\newcommand{{\\NsecondHalf}}{{{n2}}}            % sites inscribed in second half of UNESCO record")
print(f"  \\newcommand{{\\NapSecondHalf}}{{{nap2}}}           % A+ sites in second half")
print(f"  \\newcommand{{\\secondHalfRate}}{{{100*nap2/n2:.1f}}}           % A+ rate in second half (%)")
print(f"  \\newcommand{{\\pSecondHalf}}{{{p2:.6f}}}       % p-value, second-half A+ binomial")
print(f"  \\newcommand{{\\halfFisherOR}}{{{or3:.2f}}}          % Fisher OR, first vs second half A+")
print(f"  \\newcommand{{\\pHalfFisher}}{{{p3:.6f}}}       % p-value, Fisher half-corpus comparison")

# ── Write to results store ────────────────────────────────────────────────────
ResultsStore().write_many({
    "pFoundDateSpearman": p_sp2,       # Spearman p, founding date vs A+
    "pPreCEfisher":       p_fisher,    # Fisher pre- vs post-CE A+
    "pPreCEbinom":        p_pre_binom, # binomial p, pre-CE A+
    "peakSigP":           peak_p,      # p at peak sequential signal
    "pFirstHalf":         p1,          # first-half A+ binomial p
    "pSecondHalf":        p2,          # second-half A+ binomial p
    "pHalfFisher":        p3,          # Fisher p, half-corpus comparison
})

print("\n  DONE.")
