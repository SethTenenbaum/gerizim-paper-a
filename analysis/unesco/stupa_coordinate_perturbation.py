"""
stupa_coordinate_perturbation.py
==================================
Coordinate-perturbation sensitivity check for the stupa sub-population.

MOTIVATION
----------
The stupa corpus is small (N = NevoStupa ≈ 18).  The Tier-A+ window is
narrow (≤ 6.7 km equatorial).  A reviewer correctly notes that individual
coordinate errors could shift a site across the tier boundary and
materially affect the result.

This script asks: is the stupa A+ result robust to small random errors
in site coordinates consistent with the precision of UNESCO XML centroids?

DESIGN
------
Displacement model:
  Each stupa site longitude is independently perturbed by drawing from
  Uniform(−δ, +δ) where δ = 5 km equatorial ≈ 0.0450° at the equator.
  (At 27°N, the median stupa latitude, δ_lon = 5 / (111 × cos(27°)) ≈ 0.0505°.
  We use the conservative equatorial approximation δ = 0.0450°, i.e.
  5 / 111.0 degrees, consistent with the paper's equatorial-approximation
  convention throughout.)

  This captures ±5 km centroid uncertainty — a realistic upper bound for
  UNESCO XML coordinates, which represent administrative or park centroids
  and may differ from the primary architectural feature by 1–3 km.

Statistic:
  For each of N_TRIALS perturbations, compute the A+ count in the perturbed
  stupa corpus and run a one-sided binomial test (H₁: rate > 4%).  Record
  whether p < 0.05.

Outputs:
  - Mean and SD of perturbed A+ count over all trials.
  - Fraction of trials retaining p < 0.05.
  - LaTeX macros for the paper.

USAGE
-----
    python3 analysis/unesco/stupa_coordinate_perturbation.py

MANUSCRIPT MACROS PRODUCED
---------------------------
    \\stupaCoordPerturbTrials    — number of perturbation trials
    \\stupaCoordPerturbMean      — mean perturbed A+ count
    \\stupaCoordPerturbStd       — SD of perturbed A+ count
    \\stupaCoordPerturbSigPct    — % of trials retaining p < 0.05
"""

from __future__ import annotations

import sys
from pathlib import Path
from scipy.stats import binomtest  # type: ignore

import numpy as np

np.random.seed(42)

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, TIER_APLUS, P_NULL_AP
from lib.beru import deviation as _beru_deviation
from lib.results_store import ResultsStore

import re
import json

# Load stupa keywords from keywords.json (same source as tumulus_dome_evolution_raw_sweep.py)
_KW_PATH = Path(__file__).parent.parent.parent / "keywords.json"
with open(_KW_PATH) as _f:
    _KW = json.load(_f)
STUPA_KEYWORDS = _KW["mound_evolution"]["stupa"]
STUPA_KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    for kw in STUPA_KEYWORDS
}

N_TRIALS = 10_000
# ±5 km displacement as an equatorial longitude offset
DELTA_DEG = 5.0 / 111.0   # ≈ 0.0450°


def beru_dev(lon: float) -> float:
    return _beru_deviation(lon)


def is_stupa_site(site) -> bool:
    text = getattr(site, "full_text", "") or ""
    return any(rx.search(text) for rx in STUPA_KEYWORD_RES.values())


def count_aplus(lons: np.ndarray) -> int:
    """Vectorized A+ counter — no Python loop over sites."""
    arr = lons if isinstance(lons, np.ndarray) else np.asarray(lons, dtype=float)
    arc = np.abs(arr - GERIZIM) / BERU
    dev = np.abs(arc - np.round(arc / 0.1) * 0.1)
    return int(np.sum(dev <= TIER_APLUS))


def binomial_p(n_ap: int, n_total: int) -> float:
    """One-sided binomial p-value (H₁: rate > P_NULL_AP)."""
    try:
        return float(binomtest(n_ap, n_total, P_NULL_AP, alternative="greater").pvalue)
    except TypeError:
        # scipy ≥ 1.11 uses binomtest
        from scipy.stats import binomtest as _bt
        return float(_bt(n_ap, n_total, P_NULL_AP, alternative="greater").pvalue)


def main():
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    stupa_lons = np.array([
        s.longitude for s in cultural if is_stupa_site(s)
    ])
    N_stupa = len(stupa_lons)
    obs_ap = count_aplus(stupa_lons)
    obs_rate = 100 * obs_ap / N_stupa
    obs_p = binomial_p(obs_ap, N_stupa)

    print("=" * 80)
    print("  STUPA COORDINATE-PERTURBATION SENSITIVITY")
    print(f"  N_stupa = {N_stupa}, observed A+ = {obs_ap} ({obs_rate:.1f}%)")
    print(f"  Observed one-sided binomial p = {obs_p:.4f}")
    print(f"  Perturbation: ±{DELTA_DEG:.4f}° per site (≈ ±5 km equatorial)")
    print(f"  N_trials = {N_TRIALS:,}, seed = 42")
    print("=" * 80)

    rng = np.random.default_rng(42)
    perturbed_ap = np.zeros(N_TRIALS, dtype=int)
    sig_count = 0

    for i in range(N_TRIALS):
        offsets = rng.uniform(-DELTA_DEG, DELTA_DEG, size=N_stupa)
        perturbed = stupa_lons + offsets
        ap = count_aplus(perturbed)
        perturbed_ap[i] = ap
        if binomial_p(ap, N_stupa) < 0.05:
            sig_count += 1

    mean_ap = float(perturbed_ap.mean())
    std_ap = float(perturbed_ap.std())
    sig_pct = round(100.0 * sig_count / N_TRIALS, 1)

    print(f"\n  Results over {N_TRIALS:,} trials:")
    print(f"    Perturbed A+ mean = {mean_ap:.2f} ± {std_ap:.2f}")
    print(f"    Trials retaining p < 0.05: {sig_count:,} / {N_TRIALS:,} = {sig_pct}%")
    print(f"\n  INTERPRETATION:")
    if sig_pct >= 80:
        verdict = "robust"
        print(f"    {sig_pct}% of perturbation trials retain significance (p < 0.05).")
        print(f"    The stupa result is robust to ±5 km coordinate uncertainty.")
    elif sig_pct >= 50:
        verdict = "moderately robust"
        print(f"    {sig_pct}% of perturbation trials retain significance.")
        print(f"    The result is moderately robust; some site coordinates are")
        print(f"    close to the tier boundary.")
    else:
        verdict = "sensitive"
        print(f"    Only {sig_pct}% of perturbation trials retain significance.")
        print(f"    The result is sensitive to coordinate precision at this N.")

    # ── LaTeX macros ─────────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("  LATEX MACROS (stupa coordinate-perturbation sensitivity):")
    print("=" * 80)
    print(f"  \\newcommand{{\\stupaCoordPerturbTrials}}{{{N_TRIALS:,}}}")
    print(f"  \\newcommand{{\\stupaCoordPerturbMean}}{{{mean_ap:.2f}}}")
    print(f"  \\newcommand{{\\stupaCoordPerturbStd}}{{{std_ap:.2f}}}")
    print(f"  \\newcommand{{\\stupaCoordPerturbSigPct}}{{{sig_pct:.0f}}}")

    # ── Store ─────────────────────────────────────────────────────────────────
    ResultsStore().write_many({
        "stupaCoordPerturbTrials": float(N_TRIALS),
        "stupaCoordPerturbMean":   round(mean_ap, 2),
        "stupaCoordPerturbStd":    round(std_ap, 2),
        "stupaCoordPerturbSigPct": sig_pct,
    })
    print("\nResults written to data/store/results.json")


if __name__ == "__main__":
    main()
