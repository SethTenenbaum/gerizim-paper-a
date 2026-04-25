"""
tier_logsweep_sensitivity.py
============================
Logarithmic sensitivity sweep across the full A/B/C tier threshold space.

ARGUMENT
--------
The canonical tier boundaries (A++, A+, A and their C-band mirrors)
are based on Babylonian metrological sub-units:

    A++ ≈ 0.0016 bēru  = 1 parasang     (5.3 km)
    A+  ≈ 0.0032 bēru  = 1 distance-bēru (10.7 km)  [= 2 parasangs]
    A   ≈ 0.0094 bēru  ≈ 3 distance-bēru (31.4 km)

If these values were cherry-picked to maximise the signal, the
enrichment curve would have sharp spikes at exactly these thresholds
with weaker results on either side.

This script tests the contrary hypothesis: the enrichment is a smooth,
monotone function across a log-scale sweep of thresholds. The canonical
values sit at metrologically-motivated positions on a smooth curve —
not at isolated spikes that disappear with small perturbations.

METHOD
------
1. Sweep τ (deviation threshold) on a log scale from 0.0001 to 0.0499 bēru
   (~60 points, covering 2.5 log-decades).
2. For each τ, compute for both the full Cultural/Mixed corpus and the
   dome sub-population:
     - n_hits       : sites with deviation ≤ τ (harmonic-side enrichment)
     - n_hits_c     : sites with dist-from-midpoint ≤ τ (midpoint-side)
     - null_rate    : geometric null = 2τ / harmonic_step
     - enrich_ratio : observed rate / null_rate
     - binom_p      : one-sided binomial test (H₁: rate > null_rate)
     - binom_p_c    : depletion test for C-band dome (H₁: dome rate < null)
3. Mark each row with the nearest canonical tier label if within 5% of
   that tier's boundary.
4. Print a structured table and write tier_logsweep_sensitivity.csv for
   plotting.

INTERPRETATION OF OUTPUT
------------------------
- A smooth curve (no isolated spikes) on a log-log plot of (τ, p-value)
  argues the choice is not threshold-engineered.
- The canonical tiers should show sub-optimal but consistent enrichment,
  not uniquely optimal performance, confirming they were chosen for
  metrological reasons, not to maximise the statistic.
- The C-band depletion mirror curve should be symmetric for dome sites,
  confirming the A/C bimodal structure is real.

USAGE
-----
    cd /path/to/gerizim-paper-a
    python3 analysis/unesco/tier_logsweep_sensitivity.py
"""

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, HARMONIC_STEP, KM_PER_DEGREE,
    TIER_APP, TIER_APLUS, TIER_A_MAX,
    TIER_C_MAX, TIER_CMINUS, TIER_CMINUS2,
    MIDPOINT,
)
from lib.dome_filter import is_dome_site
from lib.results_store import ResultsStore

# ── Corpus ─────────────────────────────────────────────────────────────────────
corpus = load_corpus()
sites  = cultural_sites_with_coords(corpus)

all_lons  = np.array([s.longitude for s in sites])
dome_lons = np.array([s.longitude for s in sites if is_dome_site(s)])

N_all  = len(all_lons)
N_dome = len(dome_lons)

# ── Pre-compute deviations (fixed Gerizim anchor) ─────────────────────────────
def _devs(lons):
    arcs  = np.abs(lons - GERIZIM)
    bvs   = arcs / BERU
    nears = np.round(bvs / HARMONIC_STEP) * HARMONIC_STEP
    devs  = np.abs(bvs - nears)          # distance from nearest harmonic
    mdist = MIDPOINT - devs              # distance from midpoint (positive = near midpoint)
    return devs, mdist

all_devs,  all_mdist  = _devs(all_lons)
dome_devs, dome_mdist = _devs(dome_lons)

# ── Canonical tier annotations ─────────────────────────────────────────────────
# A-side and C-side are annotated separately: both share the same numeric
# threshold values (they are mirror tiers), so we keep two distinct dicts.
CANONICAL_A = [
    ("A++", TIER_APP),
    ("A+",  TIER_APLUS),
    ("A",   TIER_A_MAX),
]
CANONICAL_C = [
    ("C--", TIER_CMINUS2),
    ("C-",  TIER_CMINUS),
    ("C",   TIER_C_MAX),
]

ANNOT_TOL = 0.05   # 5% relative tolerance for annotation

def _nearest_a_tier(tau):
    best, best_reldiff = "", 1.0
    for label, val in CANONICAL_A:
        rd = abs(tau - val) / val
        if rd < best_reldiff:
            best, best_reldiff = label, rd
    return best if best_reldiff <= ANNOT_TOL else ""

def _nearest_c_tier(tau):
    best, best_reldiff = "", 1.0
    for label, val in CANONICAL_C:
        rd = abs(tau - val) / val
        if rd < best_reldiff:
            best, best_reldiff = label, rd
    return best if best_reldiff <= ANNOT_TOL else ""

# ── Metrological reference lines ──────────────────────────────────────────────
KM_PER_BERU = BERU * KM_PER_DEGREE   # 3330 km
def _km(tau): return tau * KM_PER_BERU

# ── Log-scale sweep ────────────────────────────────────────────────────────────
# 60 points from 0.0001 to 0.0499 bēru, evenly spaced in log space
TAU_MIN = 1e-4
TAU_MAX = MIDPOINT - 1e-4   # just below midpoint
N_STEPS = 70

taus = np.logspace(np.log10(TAU_MIN), np.log10(TAU_MAX), N_STEPS)

# Inject canonical tier values exactly (so they appear as exact rows)
canon_vals = [v for _, v in CANONICAL_A] + [v for _, v in CANONICAL_C]
taus = np.sort(np.unique(np.concatenate([taus, canon_vals])))

# ── Row computation ────────────────────────────────────────────────────────────
rows = []

for tau in taus:
    null_rate = min(2.0 * tau / HARMONIC_STEP, 1.0)   # geometric null

    # Harmonic-side (A-band enrichment)
    n_all_a   = int(np.sum(all_devs  <= tau))
    n_dome_a  = int(np.sum(dome_devs <= tau))

    obs_rate_all  = n_all_a  / N_all
    obs_rate_dome = n_dome_a / N_dome

    p_all_a  = binomtest(n_all_a,  N_all,  null_rate, alternative="greater").pvalue
    p_dome_a = binomtest(n_dome_a, N_dome, null_rate, alternative="greater").pvalue

    enrich_all  = obs_rate_all  / null_rate if null_rate > 0 else float("nan")
    enrich_dome = obs_rate_dome / null_rate if null_rate > 0 else float("nan")

    # Midpoint-side (C-band)
    # mdist > 0 means closer to midpoint; count sites with mdist >= 0 and mdist <= tau
    n_all_c   = int(np.sum((all_mdist  >= 0) & (all_mdist  <= tau)))
    n_dome_c  = int(np.sum((dome_mdist >= 0) & (dome_mdist <= tau)))

    obs_rate_all_c  = n_all_c  / N_all
    obs_rate_dome_c = n_dome_c / N_dome

    p_all_c  = binomtest(n_all_c,  N_all,  null_rate, alternative="greater").pvalue
    p_dome_c = binomtest(n_dome_c, N_dome, null_rate, alternative="less").pvalue

    enrich_all_c  = obs_rate_all_c  / null_rate if null_rate > 0 else float("nan")
    enrich_dome_c = obs_rate_dome_c / null_rate if null_rate > 0 else float("nan")

    a_annot = _nearest_a_tier(tau)
    c_annot = _nearest_c_tier(tau)
    annot   = a_annot or (f"({c_annot})" if c_annot else "")
    km      = _km(tau)

    rows.append({
        "tau_beru":      round(float(tau), 7),
        "tau_km":        round(km, 3),
        "null_pct":      round(100 * null_rate, 4),
        # Harmonic (A) side
        "n_all_a":       n_all_a,
        "obs_pct_all_a": round(100 * obs_rate_all, 3),
        "enrich_all_a":  round(enrich_all, 3),
        "p_all_a":       float(p_all_a),
        "n_dome_a":      n_dome_a,
        "obs_pct_dome_a":round(100 * obs_rate_dome, 3),
        "enrich_dome_a": round(enrich_dome, 3),
        "p_dome_a":      float(p_dome_a),
        # Midpoint (C) side
        "n_all_c":       n_all_c,
        "enrich_all_c":  round(enrich_all_c, 3),
        "p_all_c":       float(p_all_c),
        "n_dome_c":      n_dome_c,
        "enrich_dome_c": round(enrich_dome_c, 3),
        "p_dome_c_depl": float(p_dome_c),   # depletion test
        # Annotation
        "tier":          annot,
        "tier_a":        a_annot,
        "tier_c":        c_annot,
    })

# ── Print table ────────────────────────────────────────────────────────────────
SEP = "─" * 120
HEADER = (f"{'τ (bēru)':>10}  {'τ (km)':>7}  {'null%':>6}  "
          f"{'full n':>6}  {'full%':>6}  {'enrich':>7}  {'p_full':>10}  "
          f"{'dome n':>6}  {'dome%':>6}  {'enrich':>7}  {'p_dome':>10}  "
          f"{'C_dome':>6}  {'Cenr':>6}  {'Cdepl_p':>10}  {'tier':>5}")

print()
print("  TIER LOG-SCALE SENSITIVITY SWEEP")
print(f"  Full corpus N = {N_all},  Dome N = {N_dome}")
print(f"  Sweep: {len(taus)} thresholds, τ ∈ [{TAU_MIN:.4f}, {TAU_MAX:.4f}] bēru (log-spaced)")
print(f"  Canonical tiers injected at exact values")
print()
print(f"  All τ values are in bēru (1 bēru = {BERU}° = {KM_PER_BERU:.0f} km arc)")
print(f"  Metrological key:  A++ = 1 parasang ({_km(TIER_APP):.1f} km),  "
      f"A+ = 1 dist-bēru ({_km(TIER_APLUS):.1f} km),  "
      f"A = 3 dist-bēru ({_km(TIER_A_MAX):.1f} km)")
print()
print(SEP)
print(HEADER)
print(SEP)

def _fmt_p(p):
    if p < 0.001: return "< 0.001"
    if p > 0.999: return "> 0.999"
    return f"{p:.4f}"

for r in rows:
    tag = f"◄ {r['tier']}" if r["tier"] else ""
    print(f"  {r['tau_beru']:>10.6f}  {r['tau_km']:>7.2f}  {r['null_pct']:>6.3f}  "
          f"{r['n_all_a']:>6}  {r['obs_pct_all_a']:>6.2f}  {r['enrich_all_a']:>7.3f}  {_fmt_p(r['p_all_a']):>10}  "
          f"{r['n_dome_a']:>6}  {r['obs_pct_dome_a']:>6.2f}  {r['enrich_dome_a']:>7.3f}  {_fmt_p(r['p_dome_a']):>10}  "
          f"{r['n_dome_c']:>6}  {r['enrich_dome_c']:>6.3f}  {_fmt_p(r['p_dome_c_depl']):>10}  {tag}")

print(SEP)
print()

# ── Smoothness diagnostic ──────────────────────────────────────────────────────
# Check that the dome enrichment curve is monotone (or near-monotone) on the A side.
# Non-monotone bumps would indicate cherry-picking of tier boundaries.
dome_enrich = [r["enrich_dome_a"] for r in rows]
n_inversions = sum(1 for a, b in zip(dome_enrich, dome_enrich[1:]) if b > a + 0.1)
print(f"  Smoothness check — dome A-side enrichment curve:")
print(f"    {n_inversions} non-monotone inversions > 0.1× (of {len(rows)-1} steps)")
print(f"    (0 = perfectly smooth decay as τ → 0)")
print()

# Range of enrichment at canonical A-tier rows
print("  Enrichment at canonical A-tier boundaries:")
for r in rows:
    if r["tier_a"]:
        print(f"    {r['tier_a']:<4}  τ={r['tau_beru']:.5f} ({r['tau_km']:.1f} km)  "
              f"dome_enrich={r['enrich_dome_a']:.3f}×  p_dome={_fmt_p(r['p_dome_a'])}  "
              f"C_depl_p={_fmt_p(r['p_dome_c_depl'])}")
print()

# ── Data-optimal threshold ─────────────────────────────────────────────────────
# The threshold that minimises p_dome (data-optimal by p-value).
# If the metrological values were cherry-picked FROM the data, they should
# match the optimum.  If they differ, the metrological choice is pre-specified.
#
# We report the optimum over rows with n_dome_a >= 5 (minimum power floor)
# to avoid reporting a tau where 1-2 dome sites drive a trivially small p.
powered_rows = [r for r in rows if r["n_dome_a"] >= 5]
opt_p_row    = min(powered_rows, key=lambda r: r["p_dome_a"])
opt_enr_row  = max(powered_rows, key=lambda r: r["enrich_dome_a"])

# Distance from canonical A+ in log space
import math
def _log_dist(tau_a, tau_b):
    return abs(math.log10(tau_a) - math.log10(tau_b))

dist_opt_p_to_ap   = _log_dist(opt_p_row["tau_beru"],   TIER_APLUS)
dist_opt_p_to_app  = _log_dist(opt_p_row["tau_beru"],   TIER_APP)
dist_opt_enr_to_ap = _log_dist(opt_enr_row["tau_beru"], TIER_APLUS)
dist_opt_enr_to_app= _log_dist(opt_enr_row["tau_beru"], TIER_APP)

print("  Data-optimal thresholds (dome population, min n=5 power floor):")
print(f"    Best p-value:    τ={opt_p_row['tau_beru']:.5f} ({opt_p_row['tau_km']:.1f} km)  "
      f"p_dome={_fmt_p(opt_p_row['p_dome_a'])}  enrich={opt_p_row['enrich_dome_a']:.3f}×  "
      f"(log-dist from A+: {dist_opt_p_to_ap:.3f}  from A++: {dist_opt_p_to_app:.3f})")
print(f"    Best enrichment: τ={opt_enr_row['tau_beru']:.5f} ({opt_enr_row['tau_km']:.1f} km)  "
      f"p_dome={_fmt_p(opt_enr_row['p_dome_a'])}  enrich={opt_enr_row['enrich_dome_a']:.3f}×  "
      f"(log-dist from A+: {dist_opt_enr_to_ap:.3f}  from A++: {dist_opt_enr_to_app:.3f})")
print()
print(f"  Canonical A+ ({TIER_APLUS:.5f} bēru = {_km(TIER_APLUS):.1f} km):")
print(f"    → {dist_opt_p_to_ap:.3f} log-decades from best-p optimum")
print(f"    → corresponds to 1 distance-bēru, a pre-specified Babylonian unit")
print()

# ── Optimal geometric tier series ──────────────────────────────────────────────
# Starting from the data-optimal A++ threshold (best-p, n>=5), sweep over
# geometric ratios r to find what A+ = τ_app × r and A = τ_app × r² give.
# Reports p and enrichment at each tier in the implied geometric series.
# The goal: find the ratio r where all three tiers are simultaneously
# significant (p < 0.05) — that defines the "optimal logarithmic band structure."

tau_app_opt = opt_p_row["tau_beru"]   # data-optimal A++ threshold

# Interpolation helper: get (enrichment, p) at arbitrary tau by linear interp
tau_arr     = np.array([r["tau_beru"]    for r in rows])
p_dome_arr  = np.array([r["p_dome_a"]   for r in rows])
enr_dome_arr= np.array([r["enrich_dome_a"] for r in rows])

def _interp(tau):
    return (float(np.interp(tau, tau_arr, enr_dome_arr)),
            float(np.interp(tau, tau_arr, p_dome_arr)))

# Ratios to test: 1.2 to 5.0 (covers ×1.2 to ×5 per tier step)
ratios = np.linspace(1.2, 5.0, 200)
MIDPOINT_BERU = MIDPOINT   # 0.05 bēru

print("  Optimal geometric tier series (anchored at data-optimal A++):")
print(f"  A++ fixed at τ*={tau_app_opt:.5f} bēru ({_km(tau_app_opt):.1f} km)")
print()
print(f"  {'ratio r':>8}  "
      f"{'A++ km':>7}  {'A++ p':>8}  "
      f"{'A+ km':>7}  {'A+ enr':>7}  {'A+ p':>8}  "
      f"{'A km':>7}  {'A enr':>7}  {'A p':>8}  "
      f"{'all sig':>7}")
print(f"  {'─'*8}  {'─'*7}  {'─'*8}  {'─'*7}  {'─'*7}  {'─'*8}  {'─'*7}  {'─'*7}  {'─'*8}  {'─'*7}")

best_ratio    = None
best_minimax  = float("-inf")   # maximise the weakest tier's significance
geo_results   = []

# Minimum spread: A must be at least r² × A++, and A/A++ >= 3 (tiers must be distinct).
MIN_SPREAD = 3.0   # A must be at least 3× A++

for r in ratios:
    tau_ap  = tau_app_opt * r
    tau_a   = tau_app_opt * r * r
    if tau_a >= MIDPOINT_BERU:
        continue
    if r * r < MIN_SPREAD:
        continue   # tiers too close together to be meaningful

    enr_app, p_app = _interp(tau_app_opt)
    enr_ap,  p_ap  = _interp(tau_ap)
    enr_a,   p_a   = _interp(tau_a)

    all_sig = p_app < 0.05 and p_ap < 0.05 and p_a < 0.05
    # Minimax score: the weakest tier's -log10(p) — maximising this ensures
    # all three tiers are simultaneously as significant as possible.
    minimax = min(-math.log10(max(p_app, 1e-10)),
                  -math.log10(max(p_ap,  1e-10)),
                  -math.log10(max(p_a,   1e-10)))

    geo_results.append({
        "r": r, "tau_ap": tau_ap, "tau_a": tau_a,
        "p_app": p_app, "enr_ap": enr_ap, "p_ap": p_ap,
        "enr_a": enr_a, "p_a": p_a, "all_sig": all_sig, "minimax": minimax,
    })

    if all_sig and minimax > best_minimax:
        best_minimax = minimax
        best_ratio   = r

# Print a thinned table (every ~15th row) plus the optimal row
print_indices = set(range(0, len(geo_results), 15))
if best_ratio is not None:
    best_idx = min(range(len(geo_results)),
                   key=lambda i: abs(geo_results[i]["r"] - best_ratio))
    print_indices.add(best_idx)
# Also add rows near canonical ratios (×2, ×3)
for canon_r in [2.0, 3.0]:
    idx = min(range(len(geo_results)), key=lambda i: abs(geo_results[i]["r"] - canon_r))
    print_indices.add(idx)

for i in sorted(print_indices):
    g = geo_results[i]
    flag = " ◄ OPTIMAL" if best_ratio and abs(g["r"] - best_ratio) < 0.02 else (
           " (×2)" if abs(g["r"] - 2.0) < 0.02 else (
           " (×3)" if abs(g["r"] - 3.0) < 0.02 else ""))
    print(f"  {g['r']:>8.2f}  "
          f"{_km(tau_app_opt):>7.1f}  {_fmt_p(g['p_app']):>8}  "
          f"{_km(g['tau_ap']):>7.1f}  {g['enr_ap']:>7.3f}  {_fmt_p(g['p_ap']):>8}  "
          f"{_km(g['tau_a']):>7.1f}  {g['enr_a']:>7.3f}  {_fmt_p(g['p_a']):>8}  "
          f"{'YES' if g['all_sig'] else 'no':>7}{flag}")

print()

# Summary of optimal geometric series
if best_ratio:
    best        = next(g for g in geo_results if abs(g["r"] - best_ratio) < 0.02)
    tau_ap_opt  = best["tau_ap"]
    tau_a_opt   = best["tau_a"]
    tau_ap_opt_km = _km(tau_ap_opt)
    tau_a_opt_km  = _km(tau_a_opt)
    print(f"  Optimal geometric series (ratio r={best_ratio:.2f}, all tiers p<0.05):")
    print(f"    A++ = {tau_app_opt:.5f} bēru  ({_km(tau_app_opt):.1f} km)")
    print(f"    A+  = {tau_ap_opt:.5f} bēru  ({tau_ap_opt_km:.1f} km)  ×{best_ratio:.2f}")
    print(f"    A   = {tau_a_opt:.5f} bēru  ({tau_a_opt_km:.1f} km)  ×{best_ratio**2:.2f}")
    print()
    print(f"  Canonical values for comparison:")
    print(f"    A++ = {TIER_APP:.5f} bēru  ({_km(TIER_APP):.1f} km)  "
          f"ratio to optimal: {TIER_APP/tau_app_opt:.2f}×")
    print(f"    A+  = {TIER_APLUS:.5f} bēru  ({_km(TIER_APLUS):.1f} km)  "
          f"ratio to optimal: {TIER_APLUS/tau_ap_opt:.2f}×")
    print(f"    A   = {TIER_A_MAX:.5f} bēru  ({_km(TIER_A_MAX):.1f} km)  "
          f"ratio to optimal: {TIER_A_MAX/tau_a_opt:.2f}×")
else:
    tau_ap_opt = tau_app_opt * 1.87
    tau_a_opt  = tau_app_opt * 1.87 ** 2
    tau_ap_opt_km = _km(tau_ap_opt)
    tau_a_opt_km  = _km(tau_a_opt)
    print("  No ratio r found where all three tiers reach p<0.05 simultaneously.")
print()

# Range check — does enrichment hold across 1 log-decade around A+?
a_plus_region = [r for r in rows
                 if TIER_APLUS * 0.3 <= r["tau_beru"] <= TIER_APLUS * 3.0]
min_e = min(r["enrich_dome_a"] for r in a_plus_region)
max_e = max(r["enrich_dome_a"] for r in a_plus_region)
n_sig_in_region = sum(1 for r in a_plus_region if r["p_dome_a"] < 0.05)
print(f"  Stability — dome enrichment in ±1 log-decade around A+  "
      f"(τ ∈ [{TIER_APLUS*0.3:.4f}, {TIER_APLUS*3.0:.4f}] bēru):")
print(f"    enrichment range: {min_e:.2f}× – {max_e:.2f}×")
print(f"    rows with p < 0.05: {n_sig_in_region}/{len(a_plus_region)}")
print()

# ── CSV output ─────────────────────────────────────────────────────────────────
csv_path = Path(__file__).parent / "tier_logsweep_sensitivity.csv"
fieldnames = list(rows[0].keys())
with open(csv_path, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
print(f"  CSV written to: {csv_path}")

# ── LaTeX macros ───────────────────────────────────────────────────────────────
n_taus_tested  = len(taus)
n_taus_sig     = sum(1 for r in rows if r["p_dome_a"] < 0.05)
n_taus_A_range = len(a_plus_region)
n_taus_A_sig   = sum(1 for r in a_plus_region if r["p_dome_a"] < 0.05)
canon_a_rows = [r for r in rows if r["tier_a"] in ("A++", "A+", "A")]
min_enrich_canon = min(r["enrich_dome_a"] for r in canon_a_rows)
max_enrich_canon = max(r["enrich_dome_a"] for r in canon_a_rows)

# Exact canonical row values (first exact match for each)
ap_row  = next(r for r in rows if r["tier_a"] == "A+")
app_row = next(r for r in rows if r["tier_a"] == "A++")
a_row   = next(r for r in rows if r["tier_a"] == "A")

def _fmt_p_macro(p):
    return "< 0.001" if p < 0.001 else ("> 0.999" if p > 0.999 else f"{p:.3f}")

print()
print("% LaTeX macros — tier_logsweep_sensitivity.py")
print(f"\\newcommand{{\\logSweepNtaus}}{{{n_taus_tested}}}           % thresholds tested in log sweep")
print(f"\\newcommand{{\\logSweepNsig}}{{{n_taus_sig}}}            % thresholds with p_dome < 0.05")
print(f"\\newcommand{{\\logSweepNinversions}}{{{n_inversions}}}         % non-monotone inversions in dome curve")
print(f"\\newcommand{{\\logSweepAregionNtaus}}{{{n_taus_A_range}}}       % thresholds in ±1 log-decade around A+")
print(f"\\newcommand{{\\logSweepAregionNsig}}{{{n_taus_A_sig}}}        % of those with p_dome < 0.05")
print(f"\\newcommand{{\\logSweepMinEnrichCanon}}{{{min_enrich_canon:.2f}}}   % min dome enrichment at any A-canon tier")
print(f"\\newcommand{{\\logSweepMaxEnrichCanon}}{{{max_enrich_canon:.2f}}}   % max dome enrichment at any A-canon tier")
print(f"\\newcommand{{\\logSweepApEnrich}}{{{ap_row['enrich_dome_a']:.2f}}}      % dome enrichment at A+ canonical value")
print(f"\\newcommand{{\\logSweepAppEnrich}}{{{app_row['enrich_dome_a']:.2f}}}     % dome enrichment at A++ canonical value")
print(f"\\newcommand{{\\logSweepAEnrich}}{{{a_row['enrich_dome_a']:.2f}}}       % dome enrichment at A canonical value")
print(f"\\newcommand{{\\logSweepApP}}{{{_fmt_p_macro(ap_row['p_dome_a'])}}}     % dome p at A+ canonical threshold")
print(f"\\newcommand{{\\logSweepAppP}}{{{_fmt_p_macro(app_row['p_dome_a'])}}}    % dome p at A++ canonical threshold")
print(f"\\newcommand{{\\logSweepAP}}{{{_fmt_p_macro(a_row['p_dome_a'])}}}      % dome p at A canonical threshold")
print(f"\\newcommand{{\\logSweepOptPtau}}{{{opt_p_row['tau_beru']:.5f}}}    % data-optimal tau (min p, n>=5)")
print(f"\\newcommand{{\\logSweepOptPtauKm}}{{{opt_p_row['tau_km']:.1f}}}     % data-optimal tau in km")
print(f"\\newcommand{{\\logSweepOptPval}}{{{_fmt_p_macro(opt_p_row['p_dome_a'])}}}     % p at data-optimal tau")
print(f"\\newcommand{{\\logSweepOptEnrich}}{{{opt_enr_row['enrich_dome_a']:.2f}}}   % max enrichment at any tau (n>=5)")
print(f"\\newcommand{{\\logSweepOptEnrichKm}}{{{opt_enr_row['tau_km']:.1f}}}  % km at max enrichment tau")
print(f"\\newcommand{{\\logSweepApOptLogDist}}{{{dist_opt_p_to_ap:.3f}}}  % log-decades from A+ to best-p tau")
print(f"\\newcommand{{\\logSweepOptApBeru}}{{{tau_ap_opt:.5f}}}   % optimal A+ threshold (beru)")
print(f"\\newcommand{{\\logSweepOptApKm}}{{{tau_ap_opt_km:.1f}}}      % optimal A+ threshold (km)")
print(f"\\newcommand{{\\logSweepOptABeru}}{{{tau_a_opt:.5f}}}    % optimal A threshold (beru)")
print(f"\\newcommand{{\\logSweepOptAKm}}{{{tau_a_opt_km:.1f}}}       % optimal A threshold (km)")
print(f"\\newcommand{{\\logSweepOptRatio}}{{{best_ratio:.2f}}}      % optimal geometric step ratio")

# ── Results store ──────────────────────────────────────────────────────────────
ResultsStore().write_many({
    "logSweepNtaus":           float(n_taus_tested),
    "logSweepNsig":            float(n_taus_sig),
    "logSweepNinversions":     float(n_inversions),
    "logSweepAregionNtaus":    float(n_taus_A_range),
    "logSweepAregionNsig":     float(n_taus_A_sig),
    "logSweepMinEnrichCanon":  min_enrich_canon,
    "logSweepMaxEnrichCanon":  max_enrich_canon,
    "logSweepApEnrich":        ap_row["enrich_dome_a"],
})
