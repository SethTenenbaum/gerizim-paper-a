"""
sensitivity_slope_specificity_test.py
======================================
A stronger companion to sensitivity_slope_permutation_test.py.

MOTIVATION
----------
The existing test (sensitivity_slope_permutation_test.py) asks whether the
*general* pattern of a sharp peak anywhere in the fine window is unusual under
geographic concentration.  That test is non-significant (p = permSlopePsharp, ns),
showing geographic concentration is sufficient to reproduce the sharpness.

The STRONGER question is:

    "In what fraction of permuted datasets does the canonical 0.1000-beru spacing
     specifically produce the best joint significance in the fine window — i.e.,
     simultaneously the lowest combined p-value in both populations AND all
     off-canonical spacings non-significant in at least one population?"

This tests not just sharpness but SPECIFICITY TO 0.1000 beru.  Under geographic
concentration, a sharp peak can appear at any spacing in the fine window, not
just the metrologically attested one.  If 0.1000 is the only spacing that produces
joint significance in the real data while permuted data produces sharp peaks at
arbitrary spacings, then the *location* of the peak at the canonical unit is
informative, even if the sharpness alone is not unusual.

THREE STATISTICS ARE COMPUTED
------------------------------
(1) permSlopeCanonBestJoint:
    Fraction of permutations where 0.1000 produces joint significance AND is the
    unique best spacing (smallest combined -log p summed over dome + full).

(2) permSlopeCanonRank:
    In the real data, what is the rank of 0.1000 by combined -log p among all
    fine spacings?  (rank 1 = best, 7 = worst).  Compared with the distribution
    of ranks across permutations.

(3) permSlopeCanonBestGivenSharp:
    Of the permutations that DO show a sharp-peak pattern anywhere in the window,
    what fraction have the peak specifically at 0.1000?  This conditional fraction
    addresses the question: given that a sharp peak appears, how often does it
    happen to land exactly on the metrologically attested spacing?

MACROS EMITTED
--------------
    \\permSlopeSpecN             — number of permutations (same as permSlopeNperms)
    \\permSlopeCanonBestJoint    — p-value for statistic (1) above
    \\permSlopeCanonRankObs      — observed rank of 0.1000 in the real data (integer)
    \\permSlopeCanonRankPct      — percentile of that rank among permuted datasets
    \\permSlopeCanonBestGivenSharp — fraction of sharp-peak permutations where peak is at 0.1000

CODE
----
    analysis/unesco/sensitivity_slope_specificity_test.py

DEPENDENCIES
------------
    lib/beru.py, lib/dome_filter.py, lib/results_store.py
    data/unesco_corpus.py

Writes all results to the ResultsStore.
"""

import sys
import numpy as np
from pathlib import Path
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import TIER_APLUS
from lib.dome_filter import is_dome_site
from lib.results_store import ResultsStore

# ── Configuration ──────────────────────────────────────────────────────────────
N_PERMS = 10_000
SEED = 42

# Fine spacings: same set as sensitivity_slope_permutation_test.py
FINE_SPACINGS = [0.0990, 0.0993, 0.0995, 0.1000, 0.1001, 0.1003, 0.1010]
CANONICAL = 0.1000

# Off-canonical: spacings ≥ 0.0003 beru from canonical (the "collapse" threshold)
OFF_CANONICAL = [s for s in FINE_SPACINGS if abs(s - CANONICAL) >= 0.0003]

ALPHA = 0.05  # significance threshold for joint-significance tests


# ── Core helpers ───────────────────────────────────────────────────────────────

def null_rate(spacing: float) -> float:
    """Geometric null rate: window_width / harmonic_spacing."""
    return 2 * TIER_APLUS / spacing


def count_hits(lons: np.ndarray, spacing: float) -> int:
    """Count sites within TIER_APLUS of the nearest harmonic at given spacing."""
    beru_vals = np.abs(lons) / 30.0          # distance from anchor in beru
    # deviation from nearest harmonic node
    nearest = np.round(beru_vals / spacing) * spacing
    devs = np.abs(beru_vals - nearest)
    return int(np.sum(devs <= TIER_APLUS))


def sweep_scores(dome_lons: np.ndarray, all_lons: np.ndarray):
    """
    For each spacing in FINE_SPACINGS, compute dome_p, full_p, and
    the combined score -log10(dome_p) - log10(full_p) (higher = more significant).

    Returns
    -------
    scores    : list of (spacing, dome_p, full_p, combined_score)
    """
    scores = []
    for sp in FINE_SPACINGS:
        nr = null_rate(sp)
        d_hits = count_hits(dome_lons, sp)
        f_hits = count_hits(all_lons, sp)
        d_p = binomtest(d_hits, len(dome_lons), nr, alternative='greater').pvalue
        f_p = binomtest(f_hits, len(all_lons), nr, alternative='greater').pvalue
        # Clamp to avoid log(0); p is always > 0 for binomtest
        combined = -np.log10(max(d_p, 1e-15)) + -np.log10(max(f_p, 1e-15))
        scores.append((sp, d_p, f_p, combined))
    return scores


def is_jointly_significant(d_p: float, f_p: float) -> bool:
    return d_p < ALPHA and f_p < ALPHA


def sharp_peak_spacing(scores):
    """
    Return the spacing of the sharp-peak pattern if it exists, else None.
    Sharp-peak: exactly one spacing is jointly significant AND it is jointly
    significant while all spacings >= 0.0003 from canonical are not
    jointly significant.  This generalises the original test to return
    WHICH spacing is the peak (not just whether one exists at canonical).
    """
    jointly_sig_spacings = [sp for (sp, d_p, f_p, _) in scores
                             if is_jointly_significant(d_p, f_p)]
    if len(jointly_sig_spacings) == 0:
        return None
    # Any off-canonical spacing jointly significant? → not a clean sharp peak
    off_canon_joint = any(sp in OFF_CANONICAL for sp in jointly_sig_spacings)
    if off_canon_joint:
        return None
    # The only jointly-significant spacing must be near-canonical (<0.0003 from 0.1000)
    near_canon = [sp for sp in jointly_sig_spacings if abs(sp - CANONICAL) < 0.0003]
    if near_canon:
        # Return the one that has the best combined score among near-canonical
        near_scores = [(sp, comb) for (sp, d_p, f_p, comb) in scores
                       if sp in near_canon]
        return max(near_scores, key=lambda x: x[1])[0]
    return None


def rank_of_canonical(scores) -> int:
    """
    Rank of 0.1000 by combined -log p score.
    Rank 1 = highest combined score (most significant joint result).
    """
    sorted_scores = sorted(scores, key=lambda x: x[3], reverse=True)
    for rank, (sp, _, _, _) in enumerate(sorted_scores, start=1):
        if sp == CANONICAL:
            return rank
    return len(FINE_SPACINGS)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    rng = np.random.default_rng(SEED)

    # Load data
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    all_lons_abs = np.array([abs(s.longitude - 35.272)  # deviation from Gerizim anchor
                              for s in cultural])
    # Note: we work with absolute angular distances from Gerizim,
    # matching the beru deviation computation in lib/beru.py.
    # The sign of (lon - anchor) does not matter for count_hits.
    all_lons_full = np.array([s.longitude for s in cultural])
    dome_mask = np.array([is_dome_site(s) for s in cultural])

    # Recompute using the same approach as sensitivity_slope_permutation_test.py:
    # work with raw longitudes, not pre-subtracted deviations, since deviation_at_spacing
    # uses lib/beru constants internally.
    from lib.beru import GERIZIM, BERU, deviation_at_spacing

    def _dev_at_spacing_vec(lons_arr: np.ndarray, spacing: float) -> np.ndarray:
        """Vectorized deviation_at_spacing for an array of raw longitudes."""
        arc      = np.abs(lons_arr - GERIZIM) / BERU        # beru units
        nearest  = np.round(arc / spacing) * spacing
        return np.abs(arc - nearest)

    def count_hits_correct(lons: np.ndarray, spacing: float) -> int:
        """Vectorized version — identical result to the old scalar loop."""
        return int(np.sum(_dev_at_spacing_vec(np.asarray(lons), spacing) <= TIER_APLUS))

    N_all = len(all_lons_full)
    N_dome = int(dome_mask.sum())
    print(f"Loaded {N_all} Cultural/Mixed sites, {N_dome} dome/spherical sites")
    print(f"Running {N_PERMS:,} permutations (seed={SEED})\n")

    # ── Observed sweep ─────────────────────────────────────────────────────────
    dome_lons_obs = all_lons_full[dome_mask]

    def sweep_correct(dome_lons, all_lons):
        scores = []
        for sp in FINE_SPACINGS:
            nr = null_rate(sp)
            d_hits = count_hits_correct(dome_lons, sp)
            f_hits = count_hits_correct(all_lons, sp)
            d_p = binomtest(d_hits, len(dome_lons), nr, alternative='greater').pvalue
            f_p = binomtest(f_hits, len(all_lons), nr, alternative='greater').pvalue
            combined = -np.log10(max(d_p, 1e-15)) + -np.log10(max(f_p, 1e-15))
            scores.append((sp, d_p, f_p, combined))
        return scores

    print("── Observed fine sweep ──")
    obs_scores = sweep_correct(dome_lons_obs, all_lons_full)
    for sp, d_p, f_p, comb in obs_scores:
        mark = " <<<" if sp == CANONICAL else ""
        js = "JOINT-SIG" if is_jointly_significant(d_p, f_p) else ""
        print(f"  {sp:.4f}  dome_p={d_p:.4f}  full_p={f_p:.4f}  "
              f"combined={comb:.3f}  {js}{mark}")

    obs_canon_rank = rank_of_canonical(obs_scores)
    obs_sharp_peak_spacing = sharp_peak_spacing(obs_scores)
    obs_canon_is_best_joint = (obs_sharp_peak_spacing == CANONICAL)

    print(f"\nObserved:")
    print(f"  Rank of 0.1000 beru by combined -log p:  {obs_canon_rank} of {len(FINE_SPACINGS)}")
    print(f"  Sharp-peak spacing:  {obs_sharp_peak_spacing}")
    print(f"  0.1000 is unique best joint-sig spacing:  {obs_canon_is_best_joint}")

    # ── Permutation loop — vectorized over spacings ────────────────────────────
    print(f"\n── Permutation test ({N_PERMS:,} trials) — vectorized ──")

    spacings_arr   = np.array(FINE_SPACINGS)                             # (S,)
    null_rates_arr = np.array([null_rate(sp) for sp in FINE_SPACINGS])   # (S,)
    canon_idx      = FINE_SPACINGS.index(CANONICAL)
    off_mask       = np.array([sp in OFF_CANONICAL for sp in FINE_SPACINGS])  # (S,)

    def _hits_all_spacings(lons_arr: np.ndarray) -> np.ndarray:
        """Return (S,) integer hit counts for all spacings simultaneously."""
        arc     = np.abs(lons_arr[:, None] - GERIZIM) / BERU            # (n, S)
        nearest = np.round(arc / spacings_arr) * spacings_arr            # (n, S)
        dev     = np.abs(arc - nearest)                                  # (n, S)
        return (dev <= TIER_APLUS).sum(axis=0).astype(int)               # (S,)

    # stat 1: fraction where canonical is unique best joint-sig spacing
    n_canon_best_joint = 0
    # stat 2: distribution of canonical rank
    perm_canon_ranks = []
    # stat 3: conditional — of sharp-peak perms, how many have peak AT canonical?
    n_sharp_peak_total = 0
    n_sharp_peak_at_canon = 0

    for i in range(N_PERMS):
        perm_lons      = rng.permutation(all_lons_full)
        perm_dome_lons = perm_lons[dome_mask]

        d_hits = _hits_all_spacings(perm_dome_lons)   # (S,)
        f_hits = _hits_all_spacings(perm_lons)         # (S,)

        # Compute p-values for all spacings at once
        d_pvals = np.array([binomtest(int(d_hits[j]), N_dome, null_rates_arr[j],
                                      alternative='greater').pvalue
                            for j in range(len(FINE_SPACINGS))])
        f_pvals = np.array([binomtest(int(f_hits[j]), N_all, null_rates_arr[j],
                                      alternative='greater').pvalue
                            for j in range(len(FINE_SPACINGS))])
        combined   = (-np.log10(np.maximum(d_pvals, 1e-15))
                      - np.log10(np.maximum(f_pvals, 1e-15)))            # (S,)
        joint_sig  = (d_pvals < ALPHA) & (f_pvals < ALPHA)              # (S,)

        # Build scores list for compatibility with helpers
        perm_scores = [(FINE_SPACINGS[j], d_pvals[j], f_pvals[j], combined[j])
                       for j in range(len(FINE_SPACINGS))]

        # stat 1
        perm_peak_sp = sharp_peak_spacing(perm_scores)
        if perm_peak_sp == CANONICAL:
            n_canon_best_joint += 1

        # stat 2
        perm_canon_ranks.append(rank_of_canonical(perm_scores))

        # stat 3
        if perm_peak_sp is not None:
            n_sharp_peak_total += 1
            if perm_peak_sp == CANONICAL:
                n_sharp_peak_at_canon += 1

        if (i + 1) % 1000 == 0:
            print(f"  {i+1:>6}/{N_PERMS}  "
                  f"canon_best={n_canon_best_joint}  "
                  f"sharp_peaks={n_sharp_peak_total}  "
                  f"sharp_at_canon={n_sharp_peak_at_canon}")

    # ── Compute statistics ─────────────────────────────────────────────────────
    # stat 1: p-value = fraction of perms where canonical is unique best joint-sig
    p_canon_best_joint = n_canon_best_joint / N_PERMS

    # stat 2: percentile of observed rank among permuted ranks
    perm_canon_ranks = np.array(perm_canon_ranks)
    # lower rank = better; what fraction of perms have rank <= observed rank?
    # (fraction as good or better than observed)
    pct_as_good = float(np.mean(perm_canon_ranks <= obs_canon_rank))
    # p-value for "canonical ranks this well or better" = fraction with rank >= obs
    # (unusual = obs rank is unusually low = unusually high combined score)
    p_canon_rank = float(np.mean(perm_canon_ranks >= obs_canon_rank))

    # stat 3: conditional fraction
    cond_fraction = (n_sharp_peak_at_canon / n_sharp_peak_total
                     if n_sharp_peak_total > 0 else 0.0)
    n_fine_spacings = len(FINE_SPACINGS)
    # Under the null, the sharp peak (if it appears) is equally likely at any
    # of the ~2 near-canonical spacings (0.1000, 0.1001).  The null expectation
    # for the conditional fraction is approximately 1/2 (since the sharp-peak
    # criterion requires a near-canonical peak by definition).
    # We report the raw fraction for transparency.

    # ── Print results ──────────────────────────────────────────────────────────
    print(f"\n══ Results ══")
    print(f"N permutations:  {N_PERMS:,}")
    print()
    print(f"Statistic 1 — Canonical spacing (0.1000) is unique best joint-sig peak:")
    print(f"  Observed data:   {'YES' if obs_canon_is_best_joint else 'NO'}")
    print(f"  Perm fraction:   {n_canon_best_joint}/{N_PERMS}  "
          f"(p = {p_canon_best_joint:.4f})")
    print()
    print(f"Statistic 2 — Rank of 0.1000 by combined -log p among {n_fine_spacings} fine spacings:")
    print(f"  Observed rank:   {obs_canon_rank} of {n_fine_spacings}  "
          f"(rank 1 = most jointly significant)")
    print(f"  Perm distribution:  "
          f"mean={perm_canon_ranks.mean():.2f}  "
          f"median={np.median(perm_canon_ranks):.1f}  "
          f"std={perm_canon_ranks.std():.2f}")
    print(f"  P(perm rank <= obs rank):  {pct_as_good:.4f}")
    print(f"  P-value (rank as good or better):  {p_canon_rank:.4f}")
    print()
    print(f"Statistic 3 — Of {n_sharp_peak_total} sharp-peak permutations, "
          f"{n_sharp_peak_at_canon} had peak AT 0.1000:")
    print(f"  Conditional fraction:  "
          f"{n_sharp_peak_at_canon}/{n_sharp_peak_total} = {cond_fraction:.4f}")
    print(f"  (Under the null, expected ~{1/2:.2f} since only 0.1000 and 0.1001 "
          f"can produce a sharp peak by the collapse criterion)")
    print()

    # ── Format macro values ────────────────────────────────────────────────────
    def fmt_p(p: float) -> str:
        if p == 0.0:
            return f"$< {1/N_PERMS:.4f}$"
        return f"{p:.4f}"

    def sig_label(p: float) -> str:
        if p < 0.001:
            return "***"
        elif p < 0.01:
            return "**"
        elif p < 0.05:
            return "*"
        elif p < 0.10:
            return "$\\sim$"
        else:
            return "ns"

    print("── LaTeX macros ──")
    macros = {
        "permSlopeSpecN":              N_PERMS,
        "permSlopeCanonRankObs":       obs_canon_rank,
        "permSlopeCanonBestJoint":     round(p_canon_best_joint, 4),
        "permSlopeCanonRankP":         round(p_canon_rank, 4),
        "permSlopeCanonRankPctile":    round(100 * pct_as_good, 1),
        "permSlopeCanonCondFraction":  round(cond_fraction, 3),
        "permSlopeSharpPeakTotal":     n_sharp_peak_total,
        "permSlopeSharpAtCanon":       n_sharp_peak_at_canon,
    }

    for key, val in macros.items():
        if isinstance(val, float):
            print(f"  \\newcommand{{\\{key}}}{{{val:.4f}}}")
        else:
            print(f"  \\newcommand{{\\{key}}}{{{val}}}")

    print()
    print(f"  % Significance labels for key statistics:")
    print(f"  \\newcommand{{\\permSlopeCanonBestJointSig}}{{{sig_label(p_canon_best_joint)}}}")
    print(f"  \\newcommand{{\\permSlopeCanonRankSig}}{{{sig_label(p_canon_rank)}}}")

    # ── Write to results store ─────────────────────────────────────────────────
    store_vals = {k: (float(v) if not isinstance(v, str) else v)
                  for k, v in macros.items()}
    store_vals["permSlopeCanonBestJointSig"] = sig_label(p_canon_best_joint)
    store_vals["permSlopeCanonRankSig"] = sig_label(p_canon_rank)
    ResultsStore().write_many(store_vals)
    print("\nResults written to data/store/results.json")


if __name__ == "__main__":
    main()
