"""
unit_sweep_montecarlo.py
========================
Monte Carlo test: how many spacing peaks should we expect under
geographic concentration alone?

Candidate spacing values are expressed in degrees (not beru) to avoid
circular reasoning. 3.0 degrees corresponds to the beru, but it is
tested alongside arbitrary non-beru alternatives on an even grid.

The observed number of significant spacings in the dome corpus is compared
to the null distribution obtained by drawing random dome-sized subsets
from the full UNESCO longitude distribution.

Key question: is finding significant spacing peaks surprising, or trivially
expected from the geographic distribution alone?

Run from repo root:
    python3 analysis/unesco/unit_sweep_montecarlo.py
"""

import sys
import math
from fractions import Fraction
import numpy as np
from pathlib import Path
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, TIER_APLUS, deviation_at_spacing
from lib.dome_filter import is_dome_site
from lib.stats import significance_label as sig, null_rate_at_spacing
from lib.results_store import ResultsStore

# Candidate spacing values in degrees — evenly spaced, NOT denominated in beru.
# 3.0 deg (= 1 beru sub-unit) is one candidate among many arbitrary alternatives.
CANDIDATE_DEG = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
                 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0]
BERU_PER_DEG  = 1.0 / 30.0
COARSE_SPACINGS = [d * BERU_PER_DEG for d in CANDIDATE_DEG]  # convert to beru for math

N_PERM = 10_000
ALPHA  = 0.05  # threshold for "significant spacing"
RNG    = np.random.default_rng(42)

# ── Load corpus ───────────────────────────────────────────────────────────────
corpus   = load_corpus()
cultural = cultural_sites_with_coords(corpus)

all_lons  = np.array([s.longitude for s in cultural])
dome_mask = np.array([is_dome_site(s) for s in cultural])
dome_lons = all_lons[dome_mask]

N_all  = len(all_lons)
N_dome = len(dome_lons)

print("=" * 80)
print("  UNIT SWEEP MONTE CARLO")
print(f"  N_all = {N_all}  |  N_dome = {N_dome}  |  N_perm = {N_PERM}")
print(f"  Candidate spacings (deg): {CANDIDATE_DEG}")
print("=" * 80)

# ── Helper: count significant spacings for a given longitude array ─────────
def spacing_pvals(lons):
    n = len(lons)
    out = []
    for sp_beru, sp_deg in zip(COARSE_SPACINGS, CANDIDATE_DEG):
        nr   = null_rate_at_spacing(TIER_APLUS, sp_beru)
        hits = sum(1 for lon in lons
                   if deviation_at_spacing(lon, sp_beru) <= TIER_APLUS)
        p = binomtest(hits, n, nr, alternative="greater").pvalue
        out.append((sp_deg, p, hits))
    return out


def lcm_deg(a: float, b: float) -> float:
    a_frac = Fraction(str(a)).limit_denominator()
    b_frac = Fraction(str(b)).limit_denominator()
    num = math.lcm(a_frac.numerator, b_frac.numerator)
    den = math.gcd(a_frac.denominator, b_frac.denominator)
    return float(Fraction(num, den))


def is_multiple_of_three(sp_deg: float, tol: float = 1e-9) -> bool:
    return abs(round(sp_deg / 3.0) * 3.0 - sp_deg) < tol


def overlap_with_three(sp_deg: float) -> float:
    return lcm_deg(3.0, sp_deg)


def hits_mask(lons, sp_deg: float):
    sp_beru = sp_deg * BERU_PER_DEG
    return np.array([deviation_at_spacing(lon, sp_beru) <= TIER_APLUS for lon in lons], dtype=bool)


def node_hits_mask(lons):
    node_beru = 0.1
    return np.array([deviation_at_spacing(lon, node_beru) <= TIER_APLUS for lon in lons], dtype=bool)


def spacing_pval(lons, sp_deg: float):
    sp_beru = sp_deg * BERU_PER_DEG
    nr = null_rate_at_spacing(TIER_APLUS, sp_beru)
    hits = sum(1 for lon in lons
               if deviation_at_spacing(lon, sp_beru) <= TIER_APLUS)
    p = binomtest(hits, len(lons), nr, alternative="greater").pvalue
    return sp_deg, p, hits


def n_sig(pvals, alpha=ALPHA):
    return sum(1 for _, p, _ in pvals if p < alpha)


def significant_nonthree_overlap(lons, pvals):
    node_mask = node_hits_mask(lons)
    out = []
    n = len(lons)
    for sp_deg, p, hits in pvals:
        if sp_deg == 3.0 or p >= ALPHA:
            continue
        sp_beru = sp_deg * BERU_PER_DEG
        nr = null_rate_at_spacing(TIER_APLUS, sp_beru)
        spacing_mask = hits_mask(lons, sp_deg)
        overlap_hits = int(np.count_nonzero(spacing_mask & node_mask))
        off_hits = int(np.count_nonzero(spacing_mask & ~node_mask))
        off_p = (binomtest(off_hits, hits, nr, alternative="greater").pvalue
                 if hits > 0 else 1.0)
        overlap_period = overlap_with_three(sp_deg)
        out.append({
            "spacing": sp_deg,
            "overlap_period": overlap_period,
            "total_hits": hits,
            "overlap_hits": overlap_hits,
            "off_hits": off_hits,
            "p": p,
            "off_p": off_p,
        })
    return out


# ── Observed ──────────────────────────────────────────────────────────────────
obs_dome_pvals = spacing_pvals(dome_lons)
obs_full_pvals = spacing_pvals(all_lons)

obs_dome_nsig = n_sig(obs_dome_pvals)
obs_full_nsig = n_sig(obs_full_pvals)

# Track whether non-3° spacings are significant only because they fall on a 3° multiple.
nonthree_overlap_results = significant_nonthree_overlap(dome_lons, obs_dome_pvals)
nonthree_overlap_results_full = significant_nonthree_overlap(all_lons, obs_full_pvals)

for label, pvals, nsig in [("DOME", obs_dome_pvals, obs_dome_nsig),
                             ("FULL CORPUS", obs_full_pvals, obs_full_nsig)]:
    print(f"\n  OBSERVED ({label}):")
    print(f"  {'Spacing':>8}  {'Null%':>6}  {'Hits':>5}  {'p':>10}  Sig")
    print("  " + "-" * 45)
    for sp_deg, p, hits in pvals:
        sp_beru = sp_deg * BERU_PER_DEG
        nr = null_rate_at_spacing(TIER_APLUS, sp_beru)
        print(f"  {sp_deg:>7.1f}°  {100*nr:>5.1f}%  {hits:>5}  {p:>10.4f}  {sig(p)}")
    print(f"  Significant spacings (p<{ALPHA}): {nsig} / {len(COARSE_SPACINGS)}")

if nonthree_overlap_results:
    print("\n  NON-3° SIGNIFICANT SPACINGS NOT DRIVEN BY 3° NODE OVERLAP (DOME):")
    print(f"  {'Spacing':>8}  {'Overlap':>8}  {'TotalHits':>9}  {'OvlHits':>8}  {'OffHits':>8}  {'p(total)':>10}  {'p(off)':>12}  Sig")
    print("  " + "-" * 92)
    for res in nonthree_overlap_results:
        sp_deg = res["spacing"]
        overlap_deg = res["overlap_period"]
        hits = res["total_hits"]
        overlap_hits = res["overlap_hits"]
        off_hits = res["off_hits"]
        p = res["p"]
        off_p = res["off_p"]
        print(f"  {sp_deg:>7.1f}°  {overlap_deg:>7.1f}°  {hits:>9}  {overlap_hits:>8}  {off_hits:>8}  {p:>10.4f}  {off_p:>12.4f}  {sig(off_p)}")

if nonthree_overlap_results_full:
    print("\n  NON-3° SIGNIFICANT SPACINGS NOT DRIVEN BY 3° NODE OVERLAP (FULL CORPUS):")
    print(f"  {'Spacing':>8}  {'Overlap':>8}  {'TotalHits':>9}  {'OvlHits':>8}  {'OffHits':>8}  {'p(total)':>10}  {'p(off)':>12}  Sig")
    print("  " + "-" * 92)
    for res in nonthree_overlap_results_full:
        sp_deg = res["spacing"]
        overlap_deg = res["overlap_period"]
        hits = res["total_hits"]
        overlap_hits = res["overlap_hits"]
        off_hits = res["off_hits"]
        p = res["p"]
        off_p = res["off_p"]
        print(f"  {sp_deg:>7.1f}°  {overlap_deg:>7.1f}°  {hits:>9}  {overlap_hits:>8}  {off_hits:>8}  {p:>10.4f}  {off_p:>12.4f}  {sig(off_p)}")

# ── Permutation null ──────────────────────────────────────────────────────────
# Dome null: draw N_dome without replacement from all_lons (label shuffle)
print(f"\n  Running {N_PERM} permutations (dome label shuffle)...")

perm_dome_nsig = np.zeros(N_PERM, dtype=int)

for i in range(N_PERM):
    perm_dome = all_lons[RNG.choice(N_all, size=N_dome, replace=False)]
    perm_dome_nsig[i] = n_sig(spacing_pvals(perm_dome))

def summarise(label, obs_nsig, perm_nsig):
    mean_p = perm_nsig.mean()
    sd_p   = perm_nsig.std()
    p_zero = (perm_nsig == 0).mean()
    p_ge   = (perm_nsig >= obs_nsig).mean()
    print(f"\n  NULL DISTRIBUTION ({label}, N_perm={N_PERM}):")
    print(f"  E[N_sig] = {mean_p:.2f}  SD = {sd_p:.2f}")
    print(f"  P(N_sig=0) = {p_zero:.4f}   P(N_sig>={obs_nsig}) = {p_ge:.4f}  {sig(p_ge)}")
    print(f"  Observed N_sig = {obs_nsig}")
    for k in range(max(perm_nsig) + 2):
        frac = (perm_nsig == k).mean()
        if frac > 0.001:
            print(f"    N_sig={k}: {100*frac:5.1f}%  {'█'*int(frac*40)}")
    return mean_p, sd_p, p_zero, p_ge

dome_mean, dome_sd, dome_pzero, dome_pge = summarise("DOME", obs_dome_nsig, perm_dome_nsig)

# ── LaTeX macros ──────────────────────────────────────────────────────────────
n_candidates = len(CANDIDATE_DEG)

# Extract p-values for specific spacings of interest
def get_p(pvals, deg):
    for sp_deg, p, _ in pvals:
        if abs(sp_deg - deg) < 0.01:
            return p
    return None

p_four_dome  = get_p(obs_dome_pvals, 4.0)
p_eight_dome = get_p(obs_dome_pvals, 8.0)

print(f"\n  LATEX MACROS:")
print(f"  \\newcommand{{\\mcDomeNsigObs}}{{{obs_dome_nsig}}}        % dome: observed sig spacings")
print(f"  \\newcommand{{\\mcDomeNsigMean}}{{{dome_mean:.2f}}}      % dome: E[N_sig] under null")
print(f"  \\newcommand{{\\mcDomePZero}}{{{dome_pzero:.4f}}}        % dome: P(N_sig=0)")
print(f"  \\newcommand{{\\mcDomePObsGe}}{{{dome_pge:.4f}}}         % dome: P(N_sig>=obs)")
print(f"  \\newcommand{{\\mcNcandidates}}{{{n_candidates}}}         % number of candidate spacings")
print(f"  \\newcommand{{\\mcNperms}}{{{N_PERM}}}               % number of permutations")
print(f"  \\newcommand{{\\mcFullNsigObs}}{{{obs_full_nsig}}}        % full corpus: observed sig spacings")
print(f"  \\newcommand{{\\mcFourDegDomeP}}{{{p_four_dome:.4f}}}    % dome p-value at 4 deg")
print(f"  \\newcommand{{\\mcEightDegDomeP}}{{{p_eight_dome:.4f}}}  % dome p-value at 8 deg")
print(f"  \\newcommand{{\\mcNonThreeOverlapCount}}{{{len(nonthree_overlap_results)}}}  % non-3° spacings overlapping a 3° multiple")

ResultsStore().write_many({
    "mcDomeNsigObs":   obs_dome_nsig,
    "mcDomeNsigMean":  round(dome_mean, 2),
    "mcDomePZero":     round(dome_pzero, 4),
    "mcDomePObsGe":    round(dome_pge, 4),
    "mcNcandidates":   n_candidates,
    "mcFullNsigObs":   obs_full_nsig,
    "mcFourDegDomeP":  round(p_four_dome, 4),
    "mcEightDegDomeP": round(p_eight_dome, 4),
    "mcNonThreeOverlapCount": len(nonthree_overlap_results),
})
