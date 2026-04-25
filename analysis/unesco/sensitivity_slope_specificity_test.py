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

THREE STATISTICS ARE COMPUTED (per tier: A++, A+, A)
------------------------------------------------------
(1) permSlopeCanonBestJoint<Tier>:
    Fraction of permutations where 0.1000 produces joint significance AND is the
    unique best spacing (smallest combined -log p summed over dome + full).

(2) permSlopeCanonRank<Tier>:
    In the real data, what is the rank of 0.1000 by combined -log p among all
    fine spacings?  (rank 1 = best, 7 = worst).  Compared with the distribution
    of ranks across permutations.

(3) permSlopeCanonBestGivenSharp<Tier>:
    Of the permutations that DO show a sharp-peak pattern anywhere in the window,
    what fraction have the peak specifically at 0.1000?

MACROS EMITTED (suffix APP = A++, AP = A+, A = A)
--------------------------------------------------
    \\permSlopeSpecN                     — number of permutations
    \\permSlopeCanonBestJoint{APP,AP,A}  — p-value for statistic (1)
    \\permSlopeCanonRankObs{APP,AP,A}    — observed rank of 0.1000
    \\permSlopeCanonRankPct{APP,AP,A}    — percentile of that rank
    \\permSlopeCanonBestGivenSharp{APP,AP,A} — cond. fraction

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
from lib.beru import TIER_APP, TIER_APLUS, TIER_A_MAX, GERIZIM, BERU
from lib.dome_filter import is_dome_site, UNAMBIGUOUS_KEYWORDS, FORM_KEYWORD_RES
from lib.results_store import ResultsStore

# Import evo corpus loader — defined as a function in the evolution test script.
# We import the function directly rather than executing the script.
import importlib.util as _ilu
_evo_spec = _ilu.spec_from_file_location(
    "tumulus_dome_evolution_test",
    Path(__file__).parent / "tumulus_dome_evolution_test.py",
)
_evo_mod = _ilu.module_from_spec(_evo_spec)
_evo_spec.loader.exec_module(_evo_mod)
_load_validated_corpus = _evo_mod.load_validated_corpus

# ── Configuration ──────────────────────────────────────────────────────────────
N_PERMS = 10_000
SEED = 42

FINE_SPACINGS = [0.0990, 0.0993, 0.0995, 0.1000, 0.1001, 0.1003, 0.1010]
CANONICAL = 0.1000
OFF_CANONICAL = [s for s in FINE_SPACINGS if abs(s - CANONICAL) >= 0.0003]
ALPHA = 0.05

TIERS = [
    ("APP", "A++", TIER_APP),
    ("AP",  "A+",  TIER_APLUS),
    ("A",   "A",   TIER_A_MAX),
]

# Filter variants: (macro_suffix, label)
# Masks are built in main() after corpus load.
FILTER_KEYS = [
    ("Dome", "dome filter (is_dome_site)"),
    ("Evo",  "dome-evolution filter (validated mound/stupa/dome corpus)"),
]


# ── Core helpers ───────────────────────────────────────────────────────────────

def null_rate(spacing: float, threshold: float) -> float:
    """Geometric null rate: window_width / harmonic_spacing."""
    return 2 * threshold / spacing


def _dev_at_spacing_vec(lons_arr: np.ndarray, spacing: float) -> np.ndarray:
    arc     = np.abs(lons_arr - GERIZIM) / BERU
    nearest = np.round(arc / spacing) * spacing
    return np.abs(arc - nearest)


def count_hits(lons: np.ndarray, spacing: float, threshold: float) -> int:
    return int(np.sum(_dev_at_spacing_vec(np.asarray(lons), spacing) <= threshold))


def sweep_scores(dome_lons: np.ndarray, all_lons: np.ndarray, threshold: float):
    scores = []
    for sp in FINE_SPACINGS:
        nr  = null_rate(sp, threshold)
        d_p = binomtest(count_hits(dome_lons, sp, threshold), len(dome_lons), nr,
                        alternative='greater').pvalue
        f_p = binomtest(count_hits(all_lons,  sp, threshold), len(all_lons),  nr,
                        alternative='greater').pvalue
        combined = -np.log10(max(d_p, 1e-15)) + -np.log10(max(f_p, 1e-15))
        scores.append((sp, d_p, f_p, combined))
    return scores


def is_jointly_significant(d_p: float, f_p: float) -> bool:
    return d_p < ALPHA and f_p < ALPHA


def sharp_peak_spacing(scores):
    jointly_sig = [sp for (sp, d_p, f_p, _) in scores if is_jointly_significant(d_p, f_p)]
    if not jointly_sig:
        return None
    if any(sp in OFF_CANONICAL for sp in jointly_sig):
        return None
    near_canon = [sp for sp in jointly_sig if abs(sp - CANONICAL) < 0.0003]
    if near_canon:
        near_scores = [(sp, comb) for (sp, d_p, f_p, comb) in scores if sp in near_canon]
        return max(near_scores, key=lambda x: x[1])[0]
    return None


def rank_of_canonical(scores) -> int:
    sorted_scores = sorted(scores, key=lambda x: x[3], reverse=True)
    for rank, (sp, _, _, _) in enumerate(sorted_scores, start=1):
        if sp == CANONICAL:
            return rank
    return len(FINE_SPACINGS)


def sweep_scores_single(lons: np.ndarray, threshold: float):
    """Single-population sweep — for full-corpus geographic bootstrap."""
    scores = []
    for sp in FINE_SPACINGS:
        nr  = null_rate(sp, threshold)
        p   = binomtest(count_hits(lons, sp, threshold), len(lons), nr,
                        alternative='greater').pvalue
        comb = -np.log10(max(p, 1e-15))
        scores.append((sp, p, comb))
    return scores


def sharp_peak_spacing_single(scores):
    """Sharp-peak criterion for a single population."""
    sig = [sp for (sp, p, _) in scores if p < ALPHA]
    if not sig:
        return None
    if any(sp in OFF_CANONICAL for sp in sig):
        return None
    near_canon = [sp for sp in sig if abs(sp - CANONICAL) < 0.0003]
    if near_canon:
        near_scores = [(sp, comb) for (sp, p, comb) in scores if sp in near_canon]
        return max(near_scores, key=lambda x: x[1])[0]
    return None


def run_full_corpus_tier(tier_key: str, tier_label: str, threshold: float,
                         all_lons_full: np.ndarray,
                         rng: np.random.Generator) -> dict:
    """Geographic bootstrap for the full cultural corpus.

    Permutes all longitudes (destroys geographic clustering) and asks:
    how often does 0.1000 beru remain the unique best spacing?
    """
    N_all = len(all_lons_full)

    print(f"\n{'═'*60}")
    print(f"  Filter: Full corpus  |  Tier {tier_label}  (threshold = {threshold} beru)")
    print(f"{'═'*60}")

    obs_scores  = sweep_scores_single(all_lons_full, threshold)
    obs_peak_sp = sharp_peak_spacing_single(obs_scores)
    obs_canon_is_best = (obs_peak_sp == CANONICAL)

    print("── Observed fine sweep (full corpus) ──")
    for sp, p, comb in obs_scores:
        mark = " <<<" if sp == CANONICAL else ""
        sig  = "SIG" if p < ALPHA else ""
        print(f"  {sp:.4f}  p={p:.4f}  combined={comb:.3f}  {sig}{mark}")
    print(f"\n  Sharp-peak spacing: {obs_peak_sp}")
    print(f"  0.1000 is unique best: {obs_canon_is_best}")

    spacings_arr   = np.array(FINE_SPACINGS)
    null_rates_arr = np.array([null_rate(sp, threshold) for sp in FINE_SPACINGS])

    def _hits_all_spacings(lons_arr: np.ndarray) -> np.ndarray:
        arc     = np.abs(lons_arr[:, None] - GERIZIM) / BERU
        nearest = np.round(arc / spacings_arr) * spacings_arr
        dev     = np.abs(arc - nearest)
        return (dev <= threshold).sum(axis=0).astype(int)

    n_canon_best       = 0
    n_sharp_peak_total = 0
    n_sharp_peak_at_canon = 0

    print(f"\n── Geographic bootstrap ({N_PERMS:,} trials) ──")
    for i in range(N_PERMS):
        perm_lons = rng.permutation(all_lons_full)
        f_hits    = _hits_all_spacings(perm_lons)

        f_pvals  = np.array([binomtest(int(f_hits[j]), N_all, null_rates_arr[j],
                                       alternative='greater').pvalue
                             for j in range(len(FINE_SPACINGS))])
        combined = -np.log10(np.maximum(f_pvals, 1e-15))
        perm_scores = [(FINE_SPACINGS[j], f_pvals[j], combined[j])
                       for j in range(len(FINE_SPACINGS))]

        perm_peak_sp = sharp_peak_spacing_single(perm_scores)
        if perm_peak_sp == CANONICAL:
            n_canon_best += 1
        if perm_peak_sp is not None:
            n_sharp_peak_total += 1
            if perm_peak_sp == CANONICAL:
                n_sharp_peak_at_canon += 1

        if (i + 1) % 1000 == 0:
            print(f"  {i+1:>6}/{N_PERMS}  canon_best={n_canon_best}  "
                  f"sharp_peaks={n_sharp_peak_total}")

    p_canon_best   = n_canon_best / N_PERMS
    cond_fraction  = (n_sharp_peak_at_canon / n_sharp_peak_total
                      if n_sharp_peak_total > 0 else 0.0)

    suffix = f"Full{tier_key}"
    print(f"\n══ Results — Full corpus / Tier {tier_label} ══")
    print(f"Stat 1 — 0.1000 unique best:  {n_canon_best}/{N_PERMS}  (p = {p_canon_best:.4f})")
    print(f"Stat 3 — cond fraction:  {n_sharp_peak_at_canon}/{n_sharp_peak_total} = {cond_fraction:.4f}")

    return {
        f"permSlopeCanonBestJoint{suffix}":    round(p_canon_best, 4),
        f"permSlopeCanonBestJointSig{suffix}": sig_label(p_canon_best),
        f"permSlopeSharpPeakTotal{suffix}":    n_sharp_peak_total,
        f"permSlopeSharpAtCanon{suffix}":      n_sharp_peak_at_canon,
        f"permSlopeCanonCondFraction{suffix}": round(cond_fraction, 3),
    }


def sig_label(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "$\\sim$"
    return "ns"


def fmt_p(p: float, n_perms: int) -> str:
    if p == 0.0:
        return f"$< {1/n_perms:.4f}$"
    return f"{p:.4f}"


# ── Per-tier analysis ──────────────────────────────────────────────────────────

def run_tier(tier_key: str, tier_label: str, threshold: float,
             dome_lons_obs: np.ndarray, all_lons_full: np.ndarray,
             dome_mask: np.ndarray, rng: np.random.Generator,
             filter_key: str = "Dome") -> dict:

    N_all  = len(all_lons_full)
    N_dome = int(dome_mask.sum())

    print(f"\n{'═'*60}")
    print(f"  Filter: {filter_key}  |  Tier {tier_label}  (threshold = {threshold} beru)")
    print(f"{'═'*60}")

    # Observed
    obs_scores           = sweep_scores(dome_lons_obs, all_lons_full, threshold)
    obs_canon_rank       = rank_of_canonical(obs_scores)
    obs_sharp_peak_sp    = sharp_peak_spacing(obs_scores)
    obs_canon_is_best    = (obs_sharp_peak_sp == CANONICAL)

    print("── Observed fine sweep ──")
    for sp, d_p, f_p, comb in obs_scores:
        mark = " <<<" if sp == CANONICAL else ""
        js   = "JOINT-SIG" if is_jointly_significant(d_p, f_p) else ""
        print(f"  {sp:.4f}  dome_p={d_p:.4f}  full_p={f_p:.4f}  "
              f"combined={comb:.3f}  {js}{mark}")
    print(f"\n  Rank of 0.1000: {obs_canon_rank} of {len(FINE_SPACINGS)}")
    print(f"  Sharp-peak spacing: {obs_sharp_peak_sp}")
    print(f"  0.1000 is unique best joint-sig: {obs_canon_is_best}")

    # Vectorised permutation loop
    spacings_arr   = np.array(FINE_SPACINGS)
    null_rates_arr = np.array([null_rate(sp, threshold) for sp in FINE_SPACINGS])

    def _hits_all_spacings(lons_arr: np.ndarray) -> np.ndarray:
        arc     = np.abs(lons_arr[:, None] - GERIZIM) / BERU
        nearest = np.round(arc / spacings_arr) * spacings_arr
        dev     = np.abs(arc - nearest)
        return (dev <= threshold).sum(axis=0).astype(int)

    n_canon_best_joint  = 0
    perm_canon_ranks    = []
    n_sharp_peak_total  = 0
    n_sharp_peak_at_canon = 0

    print(f"\n── Permutation test ({N_PERMS:,} trials) ──")
    for i in range(N_PERMS):
        perm_lons      = rng.permutation(all_lons_full)
        perm_dome_lons = perm_lons[dome_mask]

        d_hits = _hits_all_spacings(perm_dome_lons)
        f_hits = _hits_all_spacings(perm_lons)

        d_pvals = np.array([binomtest(int(d_hits[j]), N_dome, null_rates_arr[j],
                                      alternative='greater').pvalue
                            for j in range(len(FINE_SPACINGS))])
        f_pvals = np.array([binomtest(int(f_hits[j]), N_all,  null_rates_arr[j],
                                      alternative='greater').pvalue
                            for j in range(len(FINE_SPACINGS))])
        combined   = (-np.log10(np.maximum(d_pvals, 1e-15))
                      - np.log10(np.maximum(f_pvals, 1e-15)))
        perm_scores = [(FINE_SPACINGS[j], d_pvals[j], f_pvals[j], combined[j])
                       for j in range(len(FINE_SPACINGS))]

        perm_peak_sp = sharp_peak_spacing(perm_scores)
        if perm_peak_sp == CANONICAL:
            n_canon_best_joint += 1

        perm_canon_ranks.append(rank_of_canonical(perm_scores))

        if perm_peak_sp is not None:
            n_sharp_peak_total += 1
            if perm_peak_sp == CANONICAL:
                n_sharp_peak_at_canon += 1

        if (i + 1) % 1000 == 0:
            print(f"  {i+1:>6}/{N_PERMS}  "
                  f"canon_best={n_canon_best_joint}  "
                  f"sharp_peaks={n_sharp_peak_total}  "
                  f"sharp_at_canon={n_sharp_peak_at_canon}")

    # Statistics
    p_canon_best_joint = n_canon_best_joint / N_PERMS
    perm_canon_ranks   = np.array(perm_canon_ranks)
    pct_as_good        = float(np.mean(perm_canon_ranks <= obs_canon_rank))
    p_canon_rank       = float(np.mean(perm_canon_ranks >= obs_canon_rank))
    cond_fraction      = (n_sharp_peak_at_canon / n_sharp_peak_total
                          if n_sharp_peak_total > 0 else 0.0)

    suffix = f"{filter_key}{tier_key}"
    print(f"\n══ Results — {filter_key} / Tier {tier_label} ══")
    print(f"Stat 1 — 0.1000 unique best joint-sig:  "
          f"{n_canon_best_joint}/{N_PERMS}  (p = {p_canon_best_joint:.4f})")
    print(f"Stat 2 — rank of 0.1000:  observed={obs_canon_rank}  "
          f"perm mean={perm_canon_ranks.mean():.2f}  p={p_canon_rank:.4f}")
    print(f"Stat 3 — cond fraction:  "
          f"{n_sharp_peak_at_canon}/{n_sharp_peak_total} = {cond_fraction:.4f}")

    return {
        f"permSlopeCanonBestJoint{suffix}":      round(p_canon_best_joint, 4),
        f"permSlopeCanonBestJointSig{suffix}":   sig_label(p_canon_best_joint),
        f"permSlopeCanonRankObs{suffix}":        obs_canon_rank,
        f"permSlopeCanonRankP{suffix}":          round(p_canon_rank, 4),
        f"permSlopeCanonRankSig{suffix}":        sig_label(p_canon_rank),
        f"permSlopeCanonRankPctile{suffix}":     round(100 * pct_as_good, 1),
        f"permSlopeCanonCondFraction{suffix}":   round(cond_fraction, 3),
        f"permSlopeSharpPeakTotal{suffix}":      n_sharp_peak_total,
        f"permSlopeSharpAtCanon{suffix}":        n_sharp_peak_at_canon,
    }


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    corpus   = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    all_lons_full = np.array([s.longitude for s in cultural])

    # ── Build filter masks ─────────────────────────────────────────────────────
    dome_mask = np.array([is_dome_site(s) for s in cultural])

    # Kw mask: unambiguous keywords only (stupa/stupas/tholos) — no ambiguous
    # dome/domed/domes/spherical — the most conservative focal population.
    def _is_unamb(s):
        ft = s.full_text or ""
        return any(FORM_KEYWORD_RES[k].search(ft) for k in UNAMBIGUOUS_KEYWORDS)
    kw_mask = np.array([_is_unamb(s) for s in cultural])

    # Evo mask: sites present in the context-validated mound/stupa/dome corpus.
    evo_selected, _, _ = _load_validated_corpus()
    evo_site_ids = {(e["name"], round(e["lon"], 4)) for e in evo_selected}
    evo_mask = np.array([
        (s.site, round(s.longitude, 4)) in evo_site_ids for s in cultural
    ])

    FILTERS = [
        ("Dome", dome_mask),
        ("Kw",   kw_mask),
        ("Evo",  evo_mask),
    ]

    N_all = len(all_lons_full)
    for fkey, fmask in FILTERS:
        print(f"Filter {fkey}: {int(fmask.sum())} focal sites out of {N_all}")
    print(f"Running {N_PERMS:,} permutations per filter × tier (seed={SEED})\n")

    all_results: dict = {"permSlopeSpecN": N_PERMS}

    for fkey, fmask in FILTERS:
        focal_lons = all_lons_full[fmask]
        for tier_key, tier_label, threshold in TIERS:
            rng = np.random.default_rng(SEED)   # independent seed per cell
            tier_results = run_tier(tier_key, tier_label, threshold,
                                    focal_lons, all_lons_full, fmask, rng,
                                    filter_key=fkey)
            all_results.update(tier_results)

    # Full corpus geographic bootstrap
    print(f"\n\n{'━'*60}")
    print("  Full corpus geographic bootstrap")
    print(f"{'━'*60}")
    for tier_key, tier_label, threshold in TIERS:
        rng = np.random.default_rng(SEED)
        all_results.update(
            run_full_corpus_tier(tier_key, tier_label, threshold, all_lons_full, rng)
        )

    # Print LaTeX macros
    print("\n\n── LaTeX macros ──")
    print(f"  \\newcommand{{\\permSlopeSpecN}}{{{N_PERMS}}}")
    for fkey, _ in FILTERS:
        for tier_key, tier_label, _ in TIERS:
            suffix = f"{fkey}{tier_key}"
            print(f"\n  % {fkey} / Tier {tier_label}")
            keys = [
                f"permSlopeCanonBestJoint{suffix}",
                f"permSlopeCanonBestJointSig{suffix}",
                f"permSlopeCanonRankObs{suffix}",
                f"permSlopeCanonRankP{suffix}",
                f"permSlopeCanonRankSig{suffix}",
                f"permSlopeCanonRankPctile{suffix}",
                f"permSlopeCanonCondFraction{suffix}",
                f"permSlopeSharpPeakTotal{suffix}",
                f"permSlopeSharpAtCanon{suffix}",
            ]
            for k in keys:
                v = all_results[k]
                if isinstance(v, float):
                    print(f"  \\newcommand{{\\{k}}}{{{v:.4f}}}")
                else:
                    print(f"  \\newcommand{{\\{k}}}{{{v}}}")
    for tier_key, tier_label, _ in TIERS:
        suffix = f"Full{tier_key}"
        print(f"\n  % Full corpus / Tier {tier_label}")
        for k in [f"permSlopeCanonBestJoint{suffix}",
                  f"permSlopeCanonBestJointSig{suffix}",
                  f"permSlopeCanonCondFraction{suffix}",
                  f"permSlopeSharpPeakTotal{suffix}",
                  f"permSlopeSharpAtCanon{suffix}"]:
            v = all_results[k]
            if isinstance(v, float):
                print(f"  \\newcommand{{\\{k}}}{{{v:.4f}}}")
            else:
                print(f"  \\newcommand{{\\{k}}}{{{v}}}")

    # Write to results store
    store_vals = {k: (float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else v)
                  for k, v in all_results.items()}
    ResultsStore().write_many(store_vals)
    print("\nResults written to data/store/results.json")

    # Backward-compat: keep the old unprefixed A+ macros (Dome filter, A+ tier)
    # so existing manuscript references don't break until they are updated.
    compat = {
        "permSlopeCanonBestJoint":     all_results["permSlopeCanonBestJointDomeAP"],
        "permSlopeCanonBestJointSig":  all_results["permSlopeCanonBestJointSigDomeAP"],
        "permSlopeCanonRankObs":       all_results["permSlopeCanonRankObsDomeAP"],
        "permSlopeCanonRankP":         all_results["permSlopeCanonRankPDomeAP"],
        "permSlopeCanonRankSig":       all_results["permSlopeCanonRankSigDomeAP"],
        "permSlopeCanonRankPctile":    all_results["permSlopeCanonRankPctileDomeAP"],
        "permSlopeCanonCondFraction":  all_results["permSlopeCanonCondFractionDomeAP"],
        "permSlopeSharpPeakTotal":     all_results["permSlopeSharpPeakTotalDomeAP"],
        "permSlopeSharpAtCanon":       all_results["permSlopeSharpAtCanonDomeAP"],
        # Per-tier backward compat (old macros without filter prefix = Dome)
        "permSlopeCanonBestJointAPP":    all_results["permSlopeCanonBestJointDomeAPP"],
        "permSlopeCanonBestJointSigAPP": all_results["permSlopeCanonBestJointSigDomeAPP"],
        "permSlopeCanonBestJointAP":     all_results["permSlopeCanonBestJointDomeAP"],
        "permSlopeCanonBestJointSigAP":  all_results["permSlopeCanonBestJointSigDomeAP"],
        "permSlopeCanonBestJointA":      all_results["permSlopeCanonBestJointDomeA"],
        "permSlopeCanonBestJointSigA":   all_results["permSlopeCanonBestJointSigDomeA"],
    }
    ResultsStore().write_many(
        {k: (float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else v)
         for k, v in compat.items()}
    )


if __name__ == "__main__":
    main()
