"""
dome_footprint_window_sensitivity.py
======================================
Sensitivity analysis for the ±5° footprint window used in Null C of
dome_geographic_concentration_test.py.

MOTIVATION
----------
The ±5° footprint radius in Null C is a reasonable default but is not
uniquely motivated.  A reviewer may ask whether the Null C result depends
on this choice.  This script repeats the Null C restricted draw at
±2° and ±10° window radii, answering:

    "Does the dome morphological-specificity finding survive across
     different footprint window definitions?"

DESIGN
------
For each window radius in {2°, 10°}:
  1. Build a longitude pool from all full-corpus sites whose longitude
     falls within ±radius of any dome site longitude.
  2. Draw N_dome (context-validated) longitudes with replacement from
     the pool, 100,000 times.
  3. Count the fraction of draws that produce >= obs A+.

The ±5° result is already computed by dome_geographic_concentration_test.py
and is reported in the paper.  This script adds the ±2° and ±10° variants.

USAGE
-----
    python3 analysis/unesco/dome_footprint_window_sensitivity.py

MANUSCRIPT MACROS PRODUCED
---------------------------
    \\geoNullDomeTwoDegN        — pool size at ±2°
    \\geoNullDomeTwoDegMean     — bootstrap mean A+ at ±2°
    \\geoNullDomeTwoDegP        — p-value at ±2°
    \\geoNullDomeTwoDegStar     — significance label at ±2°
    \\geoNullDomeTenDegN        — pool size at ±10°
    \\geoNullDomeTenDegMean     — bootstrap mean A+ at ±10°
    \\geoNullDomeTenDegP        — p-value at ±10°
    \\geoNullDomeTenDegStar     — significance label at ±10°
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

np.random.seed(42)

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import TIER_APLUS, P_NULL_AP
from lib.beru import deviation as _beru_deviation
from lib.dome_filter import is_dome_site
from lib.results_store import ResultsStore

N_PERMS = 100_000
WINDOWS = [2.0, 10.0]


def beru_dev(lon: float) -> float:
    return _beru_deviation(lon)


def count_aplus(lons) -> int:
    return sum(1 for lon in lons if beru_dev(lon) <= TIER_APLUS)


def sig_label(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    if p < 0.10:
        return "~"
    return "ns"


def run_null_c(all_lons: np.ndarray, dome_lons: np.ndarray,
               radius: float, obs_ap: int, rng: np.random.Generator):
    """Run a single restricted-draw null at a given footprint radius."""
    pool = np.array([
        lon for lon in all_lons
        if any(abs(lon - d) <= radius for d in dome_lons)
    ])
    N_pool = len(pool)
    N_dome = len(dome_lons)

    counts = np.array([
        count_aplus(rng.choice(pool, size=N_dome, replace=True))
        for _ in range(N_PERMS)
    ])
    p = float(np.mean(counts >= obs_ap))
    mean = float(counts.mean())
    std = float(counts.std())
    z = (obs_ap - mean) / std if std > 0 else float("nan")
    return N_pool, mean, std, z, p


def main():
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    all_lons = np.array([s.longitude for s in cultural])
    dome_lons = np.array([s.longitude for s in cultural if is_dome_site(s)])
    N_dome = len(dome_lons)
    obs_dome_ap = count_aplus(dome_lons)
    rate_obs = 100 * obs_dome_ap / N_dome

    rng = np.random.default_rng(42)

    print("=" * 80)
    print("  DOME NULL C — FOOTPRINT WINDOW SENSITIVITY")
    print(f"  Dome corpus N = {N_dome}, observed A+ = {obs_dome_ap} ({rate_obs:.1f}%)")
    print(f"  Testing windows: {WINDOWS}")
    print("=" * 80)

    results: dict[float, tuple] = {}
    for radius in WINDOWS:
        print(f"\n  Running ±{radius:.0f}° window ({N_PERMS:,} trials)...")
        N_pool, mean, std, z, p = run_null_c(
            all_lons, dome_lons, radius, obs_dome_ap, rng
        )
        results[radius] = (N_pool, mean, std, z, p)
        label = sig_label(p)
        print(f"    Pool N = {N_pool:,}, mean A+ = {mean:.2f} ± {std:.2f}, "
              f"Z = {z:.2f}, p = {p:.4f} {label}")

    # ── LaTeX macros ─────────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("  LATEX MACROS (dome footprint window sensitivity):")
    print("=" * 80)

    r2_N, r2_mean, r2_std, r2_z, r2_p = results[2.0]
    r10_N, r10_mean, r10_std, r10_z, r10_p = results[10.0]

    macros = {
        "geoNullDomeTwoDegN":    r2_N,
        "geoNullDomeTwoDegMean": round(r2_mean, 2),
        "geoNullDomeTwoDegP":    round(r2_p, 4),
        "geoNullDomeTenDegN":    r10_N,
        "geoNullDomeTenDegMean": round(r10_mean, 2),
        "geoNullDomeTenDegP":    round(r10_p, 4),
    }

    for k, v in macros.items():
        print(f"  \\newcommand{{\\{k}}}{{{v}}}")
    print(f"  \\newcommand{{\\geoNullDomeTwoDegStar}}{{{sig_label(r2_p)}}}")
    print(f"  \\newcommand{{\\geoNullDomeTenDegStar}}{{{sig_label(r10_p)}}}")

    # ── Store ─────────────────────────────────────────────────────────────────
    store_data = {k: float(v) for k, v in macros.items()}
    ResultsStore().write_many(store_data)
    print("\nResults written to data/store/results.json")


if __name__ == "__main__":
    main()
