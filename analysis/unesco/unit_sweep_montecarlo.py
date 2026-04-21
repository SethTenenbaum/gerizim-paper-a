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
from lib.beru import GERIZIM, BERU, TIER_APLUS
from lib.dome_filter import is_dome_site
from lib.stats import significance_label as sig
from lib.results_store import ResultsStore

# Candidate spacing values in degrees — evenly spaced, NOT denominated in beru.
CANDIDATE_DEG = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
                 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0]
BERU_PER_DEG  = 1.0 / 30.0
COARSE_SPACINGS = np.array([d * BERU_PER_DEG for d in CANDIDATE_DEG])  # beru

N_PERM = 10_000
ALPHA  = 0.05
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

# ── Vectorized helpers ────────────────────────────────────────────────────────

def null_rates(spacings: np.ndarray, threshold: float = TIER_APLUS) -> np.ndarray:
    """Geometric null rate for each spacing: 2*threshold / spacing (clipped to 1)."""
    return np.minimum(2.0 * threshold / spacings, 1.0)

# Pre-compute null rates once
NULL_RATES = null_rates(COARSE_SPACINGS)  # shape (S,)

def deviation_matrix(lons: np.ndarray, spacings: np.ndarray = COARSE_SPACINGS) -> np.ndarray:
    """
    Compute deviation from nearest harmonic for every (lon, spacing) pair.

    Returns array of shape (len(lons), len(spacings)) in beru.
    """
    arc = np.abs(lons - GERIZIM) / BERU          # shape (N,)  — in beru
    # For each spacing: nearest = round(arc / sp) * sp
    # Vectorized over both axes simultaneously
    arc_col  = arc[:, np.newaxis]                  # (N, 1)
    sp_row   = spacings[np.newaxis, :]             # (1, S)
    nearest  = np.round(arc_col / sp_row) * sp_row
    return np.abs(arc_col - nearest)               # (N, S)

def hits_and_pvals(lons: np.ndarray):
    """
    Return (hits array shape (S,), p-values array shape (S,)).
    Fully vectorized over sites and spacings simultaneously.
    """
    dev   = deviation_matrix(lons)                 # (N, S)
    mask  = dev <= TIER_APLUS                      # (N, S) bool
    hits  = mask.sum(axis=0)                       # (S,)
    n     = len(lons)
    pvals = np.array([
        binomtest(int(h), n, float(nr), alternative="greater").pvalue
        for h, nr in zip(hits, NULL_RATES)
    ])
    return hits, pvals

def n_sig_vec(pvals: np.ndarray) -> int:
    return int((pvals < ALPHA).sum())


# ── Vectorized permutation null (all perms at once) ─────────────────────────

def perm_nsig_vectorized(lons_pool: np.ndarray, sample_size: int,
                         n_perm: int, replace: bool = False) -> np.ndarray:
    """
    For each permutation, draw a sample from lons_pool and count
    how many spacings are significant.

    Core computation is fully vectorized:
      - Build (n_perm × sample_size) index matrix at once
      - Compute deviation matrix for ALL permuted samples simultaneously:
        shape (n_perm * sample_size, S)  → reshape to (n_perm, sample_size, S)
      - Sum hits over sites axis, apply binomtest over (n_perm, S) matrix

    Returns nsig array of shape (n_perm,).
    """
    N_pool = len(lons_pool)
    # Sample indices: (n_perm, sample_size)
    if replace:
        idx = RNG.integers(0, N_pool, size=(n_perm, sample_size))
    else:
        idx = np.array([RNG.choice(N_pool, size=sample_size, replace=False)
                        for _ in range(n_perm)])

    # Flatten, compute deviations, reshape
    flat_lons = lons_pool[idx.ravel()]            # (n_perm * sample_size,)
    dev_flat  = deviation_matrix(flat_lons)        # (n_perm * sample_size, S)
    S = len(COARSE_SPACINGS)
    dev_3d    = dev_flat.reshape(n_perm, sample_size, S)  # (n_perm, sample_size, S)
    hits_2d   = (dev_3d <= TIER_APLUS).sum(axis=1)        # (n_perm, S)

    # Apply binomtest over (n_perm, S) — vectorized via broadcasting
    # p-values: shape (n_perm, S)
    n = sample_size
    # Use vectorized binomial CDF: P(X >= k) = 1 - binom.cdf(k-1, n, p)
    from scipy.stats import binom as _binom
    # hits_2d: (n_perm, S), NULL_RATES: (S,)
    pvals_2d = 1.0 - _binom.cdf(hits_2d - 1, n, NULL_RATES[np.newaxis, :])

    nsig = (pvals_2d < ALPHA).sum(axis=1)          # (n_perm,)
    return nsig.astype(int)


# ── Non-3° overlap helpers (used only on observed data, not in perm loop) ────

def lcm_deg(a: float, b: float) -> float:
    a_frac = Fraction(str(a)).limit_denominator()
    b_frac = Fraction(str(b)).limit_denominator()
    num = math.lcm(a_frac.numerator, b_frac.numerator)
    den = math.gcd(a_frac.denominator, b_frac.denominator)
    return float(Fraction(num, den))

def is_multiple_of_three(sp_deg: float, tol: float = 1e-9) -> bool:
    return abs(round(sp_deg / 3.0) * 3.0 - sp_deg) < tol

def node_hits_mask_vec(lons: np.ndarray) -> np.ndarray:
    """Boolean mask: sites within TIER_APLUS of a 0.1-beru node."""
    node_sp = 0.1  # beru
    arc = np.abs(lons - GERIZIM) / BERU
    nearest = np.round(arc / node_sp) * node_sp
    return np.abs(arc - nearest) <= TIER_APLUS

def significant_nonthree_overlap(lons: np.ndarray, hits_arr, pvals_arr):
    node_mask = node_hits_mask_vec(lons)
    out = []
    n = len(lons)
    for i, (sp_deg, sp_beru, nr, hits, p) in enumerate(
            zip(CANDIDATE_DEG, COARSE_SPACINGS, NULL_RATES, hits_arr, pvals_arr)):
        if abs(sp_deg - 3.0) < 1e-9 or p >= ALPHA:
            continue
        dev = deviation_matrix(lons, np.array([sp_beru]))[:, 0]
        spacing_mask = dev <= TIER_APLUS
        overlap_hits = int((spacing_mask & node_mask).sum())
        off_hits     = int((spacing_mask & ~node_mask).sum())
        off_p = (binomtest(off_hits, int(hits), float(nr), alternative="greater").pvalue
                 if hits > 0 else 1.0)
        out.append({
            "spacing": sp_deg,
            "overlap_period": lcm_deg(3.0, sp_deg),
            "total_hits": int(hits),
            "overlap_hits": overlap_hits,
            "off_hits": off_hits,
            "p": float(p),
            "off_p": off_p,
        })
    return out


def collapsed_independent_count(pvals_arr, nonthree_results):
    cnt = 0
    for sp_deg, p in zip(CANDIDATE_DEG, pvals_arr):
        if abs(sp_deg - 3.0) < 1e-9 and p < ALPHA:
            cnt += 1
            break
    for res in nonthree_results:
        if res.get("off_p", 1.0) < ALPHA:
            cnt += 1
    return cnt


# ── Vectorized collapsed-count permutation ───────────────────────────────────

def perm_collapsed_vectorized(lons_pool: np.ndarray, sample_size: int,
                              n_perm: int, replace: bool = False) -> np.ndarray:
    """
    Vectorized collapsed-independent-count permutation.
    For each perm sample, compute collapsed_independent_count fully in NumPy.

    The 'off-node' test for non-3° spacings is approximated vectorially:
    a non-3° spacing counts as independent only if its off-node hits
    (hits not within TIER_APLUS of a 0.1-beru node) form a significant
    binomial excess.
    """
    from scipy.stats import binom as _binom

    N_pool = len(lons_pool)
    S = len(COARSE_SPACINGS)
    node_sp = 0.1  # beru

    if replace:
        idx = RNG.integers(0, N_pool, size=(n_perm, sample_size))
    else:
        idx = np.array([RNG.choice(N_pool, size=sample_size, replace=False)
                        for _ in range(n_perm)])

    flat_lons = lons_pool[idx.ravel()]  # (n_perm * sample_size,)

    # Deviation matrix for all spacings: (n_perm*ss, S)
    dev_all = deviation_matrix(flat_lons, COARSE_SPACINGS)
    hit_all = (dev_all <= TIER_APLUS).reshape(n_perm, sample_size, S)  # (P, ss, S)
    hits_2d = hit_all.sum(axis=1)  # (P, S)

    # Node mask for every site in every perm
    arc_flat = np.abs(flat_lons - GERIZIM) / BERU
    nearest_node = np.round(arc_flat / node_sp) * node_sp
    on_node_flat = (np.abs(arc_flat - nearest_node) <= TIER_APLUS)  # (P*ss,)
    on_node_2d = on_node_flat.reshape(n_perm, sample_size)  # (P, ss)

    # For each spacing, compute off-node hits: sites in spacing band but NOT on node
    # hit_all[:,:,s] & ~on_node_2d
    on_node_3d  = on_node_2d[:, :, np.newaxis]  # (P, ss, 1)
    off_hits_2d = (hit_all & ~on_node_3d).sum(axis=1)  # (P, S)

    # p-values for all spacings: (P, S)
    n = sample_size
    pvals_2d  = 1.0 - _binom.cdf(hits_2d - 1,  n, NULL_RATES[np.newaxis, :])

    # p-values for off-node hits (using same null rate — conservative)
    # off_hits_2d[p, s] out of hits_2d[p, s] total hits
    # off_p[p,s] = P(X >= off_hits | n=hits, p=nr)
    # To stay vectorized, use total_hits as the "n" for off-node test
    off_pvals_2d = np.ones((n_perm, S))
    nonzero = hits_2d > 0
    off_pvals_2d[nonzero] = (
        1.0 - _binom.cdf(
            off_hits_2d[nonzero] - 1,
            hits_2d[nonzero],
            np.broadcast_to(NULL_RATES, (n_perm, S))[nonzero],
        )
    )

    # 3° index
    idx_three = next(i for i, d in enumerate(CANDIDATE_DEG) if abs(d - 3.0) < 1e-9)
    three_sig = pvals_2d[:, idx_three] < ALPHA  # (P,)

    # Non-3° independent count per perm
    three_mult = np.array([is_multiple_of_three(d) for d in CANDIDATE_DEG])  # (S,)
    non3_mask  = ~three_mult  # (S,)
    non3_sig   = (pvals_2d < ALPHA) & non3_mask[np.newaxis, :]  # (P, S)
    non3_indep = (non3_sig & (off_pvals_2d < ALPHA)).sum(axis=1)  # (P,)

    collapsed = three_sig.astype(int) + non3_indep
    return collapsed.astype(int)


# ── Observed ──────────────────────────────────────────────────────────────────
obs_dome_hits, obs_dome_pvals = hits_and_pvals(dome_lons)
obs_full_hits, obs_full_pvals = hits_and_pvals(all_lons)

obs_dome_nsig = n_sig_vec(obs_dome_pvals)
obs_full_nsig = n_sig_vec(obs_full_pvals)

nonthree_overlap_results      = significant_nonthree_overlap(dome_lons, obs_dome_hits, obs_dome_pvals)
nonthree_overlap_results_full = significant_nonthree_overlap(all_lons,  obs_full_hits, obs_full_pvals)

for label, hits_arr, pvals_arr, nsig in [
        ("DOME",        obs_dome_hits, obs_dome_pvals, obs_dome_nsig),
        ("FULL CORPUS", obs_full_hits, obs_full_pvals, obs_full_nsig)]:
    print(f"\n  OBSERVED ({label}):")
    print(f"  {'Spacing':>8}  {'Null%':>6}  {'Hits':>5}  {'p':>10}  Sig")
    print("  " + "-" * 45)
    for sp_deg, nr, hits, p in zip(CANDIDATE_DEG, NULL_RATES, hits_arr, pvals_arr):
        print(f"  {sp_deg:>7.1f}°  {100*nr:>5.1f}%  {hits:>5}  {p:>10.4f}  {sig(p)}")
    print(f"  Significant spacings (p<{ALPHA}): {nsig} / {len(COARSE_SPACINGS)}")

if nonthree_overlap_results:
    print("\n  NON-3° SIGNIFICANT SPACINGS NOT DRIVEN BY 3° NODE OVERLAP (DOME):")
    print(f"  {'Spacing':>8}  {'Overlap':>8}  {'TotalHits':>9}  {'OvlHits':>8}  {'OffHits':>8}  {'p(total)':>10}  {'p(off)':>12}  Sig")
    print("  " + "-" * 92)
    for res in nonthree_overlap_results:
        print(f"  {res['spacing']:>7.1f}°  {res['overlap_period']:>7.1f}°  {res['total_hits']:>9}  "
              f"{res['overlap_hits']:>8}  {res['off_hits']:>8}  {res['p']:>10.4f}  "
              f"{res['off_p']:>12.4f}  {sig(res['off_p'])}")

if nonthree_overlap_results_full:
    print("\n  NON-3° SIGNIFICANT SPACINGS NOT DRIVEN BY 3° NODE OVERLAP (FULL CORPUS):")
    print(f"  {'Spacing':>8}  {'Overlap':>8}  {'TotalHits':>9}  {'OvlHits':>8}  {'OffHits':>8}  {'p(total)':>10}  {'p(off)':>12}  Sig")
    print("  " + "-" * 92)
    for res in nonthree_overlap_results_full:
        print(f"  {res['spacing']:>7.1f}°  {res['overlap_period']:>7.1f}°  {res['total_hits']:>9}  "
              f"{res['overlap_hits']:>8}  {res['off_hits']:>8}  {res['p']:>10.4f}  "
              f"{res['off_p']:>12.4f}  {sig(res['off_p'])}")

# ── Permutation null — dome label shuffle ─────────────────────────────────────
print(f"\n  Running {N_PERM} permutations (dome label shuffle) — vectorized...")
perm_dome_nsig = perm_nsig_vectorized(all_lons, N_dome, N_PERM, replace=False)

def summarise(label, obs_nsig, perm_nsig):
    mean_p = perm_nsig.mean()
    sd_p   = perm_nsig.std()
    p_zero = (perm_nsig == 0).mean()
    p_ge   = (perm_nsig >= obs_nsig).mean()
    print(f"\n  NULL DISTRIBUTION ({label}, N_perm={N_PERM}):")
    print(f"  E[N_sig] = {mean_p:.2f}  SD = {sd_p:.2f}")
    print(f"  P(N_sig=0) = {p_zero:.4f}   P(N_sig>={obs_nsig}) = {p_ge:.4f}  {sig(p_ge)}")
    print(f"  Observed N_sig = {obs_nsig}")
    for k in range(int(perm_nsig.max()) + 2):
        frac = (perm_nsig == k).mean()
        if frac > 0.001:
            print(f"    N_sig={k}: {100*frac:5.1f}%  {'█'*int(frac*40)}")
    return mean_p, sd_p, p_zero, p_ge

dome_mean, dome_sd, dome_pzero, dome_pge = summarise("DOME", obs_dome_nsig, perm_dome_nsig)

# ── Collapsed (independent) significant-spacings count ───────────────────────
obs_dome_collapsed = collapsed_independent_count(obs_dome_pvals, nonthree_overlap_results)
obs_full_collapsed = collapsed_independent_count(obs_full_pvals, nonthree_overlap_results_full)

print(f"\n  COLLAPSED (INDEPENDENT) SIGNIFICANT SPACINGS:")
print(f"  DOME: {obs_dome_collapsed}   |   FULL: {obs_full_collapsed}")

# ── Permutation null for collapsed count — dome ───────────────────────────────
print(f"\n  Running {N_PERM} permutations (collapsed, dome) — vectorized...")
perm_dome_collapsed = perm_collapsed_vectorized(all_lons, N_dome, N_PERM, replace=False)
p_ge_collapsed  = (perm_dome_collapsed >= obs_dome_collapsed).mean()
p_eq1_collapsed = (perm_dome_collapsed == 1).mean()
print(f"\n  NULL DISTRIBUTION (collapsed, DOME, N_perm={N_PERM}):")
print(f"  E[collapsed] = {perm_dome_collapsed.mean():.2f}  SD = {perm_dome_collapsed.std():.2f}")
print(f"  P(collapsed>={obs_dome_collapsed}) = {p_ge_collapsed:.4f}  {sig(p_ge_collapsed)}")
print(f"  P(collapsed=1) = {p_eq1_collapsed:.4f}")
print(f"  Observed collapsed = {obs_dome_collapsed}")

# ── Permutation null for collapsed count — full corpus bootstrap ──────────────
print(f"\n  Running {N_PERM} bootstrap permutations (full corpus) — vectorized...")
perm_full_collapsed = perm_collapsed_vectorized(all_lons, N_all, N_PERM, replace=True)
p_ge_full  = (perm_full_collapsed >= obs_full_collapsed).mean()
p_eq1_full = (perm_full_collapsed == 1).mean()
print(f"\n  NULL DISTRIBUTION (collapsed, FULL, N_perm={N_PERM}):")
print(f"  E[collapsed] = {perm_full_collapsed.mean():.2f}  SD = {perm_full_collapsed.std():.2f}")
print(f"  P(collapsed>={obs_full_collapsed}) = {p_ge_full:.4f}  {sig(p_ge_full)}")
print(f"  P(collapsed=1) = {p_eq1_full:.4f}")
print(f"  Observed collapsed = {obs_full_collapsed}")
for k in range(int(perm_full_collapsed.max()) + 2):
    frac = (perm_full_collapsed == k).mean()
    if frac > 0.001:
        print(f"    N_sig={k}: {100*frac:5.1f}%  {'█'*int(frac*40)}")

# ── LaTeX macros ──────────────────────────────────────────────────────────────
n_candidates = len(CANDIDATE_DEG)

idx_four  = next(i for i, d in enumerate(CANDIDATE_DEG) if abs(d - 4.0) < 0.01)
idx_eight = next(i for i, d in enumerate(CANDIDATE_DEG) if abs(d - 8.0) < 0.01)
p_four_dome  = float(obs_dome_pvals[idx_four])
p_eight_dome = float(obs_dome_pvals[idx_eight])

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
print(f"  \\newcommand{{\\mcDomeNsigCollapsed}}{{{obs_dome_collapsed}}}  % dome: collapsed independent sig spacings")
print(f"  \\newcommand{{\\mcFullNsigCollapsed}}{{{obs_full_collapsed}}}  % full: collapsed independent sig spacings")
print(f"  \\newcommand{{\\mcDomeCollapsedPEqOne}}{{{p_eq1_collapsed:.4f}}}  % dome: P(collapsed=1)")
print(f"  \\newcommand{{\\mcFullCollapsedPEqOne}}{{{p_eq1_full:.4f}}}  % full: P(collapsed=1)")
print(f"  \\newcommand{{\\mcFullCollapsedPObsGe}}{{{p_ge_full:.4f}}}  % full: P(collapsed>=obs)")
print(f"  \\newcommand{{\\mcFullCollapsedMean}}{{{perm_full_collapsed.mean():.2f}}}  % full: E[collapsed] under null")
print(f"  \\newcommand{{\\mcFullCollapsedSD}}{{{perm_full_collapsed.std():.2f}}}  % full: SD collapsed under null")

ResultsStore().write_many({
    "mcDomeNsigObs":          obs_dome_nsig,
    "mcDomeNsigMean":         round(float(dome_mean), 2),
    "mcDomePZero":            round(float(dome_pzero), 4),
    "mcDomePObsGe":           round(float(dome_pge), 4),
    "mcNcandidates":          n_candidates,
    "mcFullNsigObs":          obs_full_nsig,
    "mcFourDegDomeP":         round(p_four_dome, 4),
    "mcEightDegDomeP":        round(p_eight_dome, 4),
    "mcNonThreeOverlapCount": len(nonthree_overlap_results),
    "mcDomeNsigCollapsed":    obs_dome_collapsed,
    "mcFullNsigCollapsed":    obs_full_collapsed,
    "mcDomeCollapsedPEqOne":  round(float(p_eq1_collapsed), 4),
    "mcFullCollapsedPEqOne":  round(float(p_eq1_full), 4),
    "mcFullCollapsedPObsGe":  round(float(p_ge_full), 4),
    "mcFullCollapsedMean":    round(float(perm_full_collapsed.mean()), 2),
    "mcFullCollapsedSD":      round(float(perm_full_collapsed.std()), 2),
})
