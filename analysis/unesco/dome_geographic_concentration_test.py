"""
dome_geographic_concentration_test.py
======================================
Geographic-concentration null for the dome/spherical sub-population.

MOTIVATION
----------
The primary permutation null (simulation_null_model.py, Approach 1) shuffles
all 1011 corpus longitudes across all sites and reads off the slots that
happen to be dome sites.  This correctly tests whether dome sites are special
*relative to the full corpus* longitude distribution, and it is significant
(p = simDomePermP).  But it does not directly answer the reviewer question:

    "If you randomly draw N sites from the same longitude distribution as
     the dome sites themselves, how often do you get >= obs A+ hits?"

That question asks whether the dome enrichment can be explained purely by the
geographic concentration of domed monuments in the Eurasian corridor.

TWO NULL MODELS — each tested at three tier thresholds (A, A+, A++)
----------------------------------------------------------------------
Null A — Within-dome bootstrap (resampling):
    Draw N_dome longitudes with replacement from the observed dome site
    longitude list.  Tested at all three tiers.

Null C — Restricted geographic draw (most conservative):
    Build a longitude pool from ALL full-corpus sites whose longitude falls
    within ±W° of any dome site longitude.  Draw N_dome from this pool.
    Tested at all three tiers and three window widths (±2°, ±5°, ±10°).

USAGE
-----
    python3 analysis/unesco/dome_geographic_concentration_test.py

IMPORTANT: ~2–3 minutes at N_PERMS=100_000.  Writes to disc.
"""

import sys
from pathlib import Path

import numpy as np

np.random.seed(42)

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (GERIZIM, BERU,
                      TIER_APP, P_NULL_APP,
                      TIER_APLUS, P_NULL_AP,
                      TIER_A_MAX, P_NULL_A)
from lib.beru import deviation as _beru_deviation
from lib.dome_filter import is_dome_site_raw
from lib.results_store import ResultsStore

N_PERMS = 100_000

# Tier definitions: label, bēru threshold, geometric null rate
TIERS = [
    ("A",   TIER_A_MAX, P_NULL_A),
    ("A+",  TIER_APLUS, P_NULL_AP),
    ("A++", TIER_APP,   P_NULL_APP),
]


# ── Helpers ──────────────────────────────────────────────────────────────────

def beru_dev(lon: float) -> float:
    return _beru_deviation(lon)


def count_tier(lons, threshold) -> int:
    """Count sites within bēru threshold of a harmonic. Vectorized."""
    arr = lons if isinstance(lons, np.ndarray) else np.asarray(lons)
    arc = np.abs(arr - GERIZIM) / BERU
    dev = np.abs(arc - np.round(arc / 0.1) * 0.1)
    return int(np.sum(dev <= threshold))


def boot_counts(mat, threshold):
    """Given (N_PERMS, N) longitude matrix, return per-draw tier counts."""
    arc = np.abs(mat - GERIZIM) / BERU
    dev = np.abs(arc - np.round(arc / 0.1) * 0.1)
    return (dev <= threshold).sum(axis=1).astype(int)


# legacy aliases
def count_aplusplus(lons): return count_tier(lons, TIER_APP)
def count_aplus(lons):     return count_tier(lons, TIER_APLUS)

_sig = lambda p: "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "~" if p < 0.10 else "ns"


# ── Load corpus ───────────────────────────────────────────────────────────────

def load_all():
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    all_lons = np.array([s.longitude for s in cultural])
    dome_lons = np.array([s.longitude for s in cultural if is_dome_site_raw(s)])
    return all_lons, dome_lons


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    all_lons, dome_lons = load_all()
    N_dome = len(dome_lons)
    rng = np.random.default_rng(42)

    # Pre-draw a single bootstrap matrix for Null A (shared across tiers)
    mat_a = rng.choice(dome_lons, size=(N_PERMS, N_dome), replace=True)

    print("=" * 80)
    print("  DOME GEOGRAPHIC-CONCENTRATION NULL TEST  (all tiers: A, A+, A++)")
    print(f"  Anchor: {GERIZIM}°E  |  BERU: {BERU}°  |  N_perms: {N_PERMS:,}")
    print(f"  Corpus: raw-sweep dome/stupa/tholos sites (is_dome_site_raw)")
    print(f"  N_dome = {N_dome}")
    print("=" * 80)

    eurasian = float(np.sum((dome_lons >= 20) & (dome_lons <= 120)) / N_dome * 100)
    full_e   = float(np.sum((all_lons  >= 20) & (all_lons  <= 120)) / len(all_lons) * 100)
    print(f"\n  Dome sites in 20°–120°E Eurasian corridor: {eurasian:.0f}%  "
          f"(full corpus: {full_e:.0f}%)")

    # ── Null A: within-dome bootstrap — all tiers ─────────────────────────────
    print("\n" + "═" * 80)
    print("  NULL A: WITHIN-DOME BOOTSTRAP  (resample dome longitudes, all tiers)")
    print("═" * 80)
    print(f"  {'Tier':<6}  {'Obs':>4}  {'Null%':>6}  {'Exp':>6}  "
          f"{'Mean':>6}  {'SD':>5}  {'Z':>6}  {'p':>9}  {'Sig':>4}")
    print("  " + "─" * 68)

    null_a_results = {}
    for tlabel, thresh, p_null_rate in TIERS:
        obs   = count_tier(dome_lons, thresh)
        exp   = N_dome * p_null_rate
        boot  = boot_counts(mat_a, thresh)
        p_val = float(np.mean(boot >= obs))
        mean  = float(boot.mean())
        std   = float(boot.std())
        z     = (obs - mean) / std if std > 0 else float("nan")
        null_a_results[tlabel] = dict(obs=obs, exp=exp, mean=mean, std=std, z=z, p=p_val)
        rate  = 100 * obs / N_dome
        print(f"  {tlabel:<6}  {obs:>4} ({rate:.0f}%)  {100*p_null_rate:.0f}%  "
              f"{exp:>6.2f}  {mean:>6.2f}  {std:>5.2f}  {z:>6.2f}  {p_val:>9.4f}  {_sig(p_val):>4}")

    # ── Null C: restricted geographic draw — all tiers × all windows ──────────
    print("\n" + "═" * 80)
    print("  NULL C: RESTRICTED GEOGRAPHIC DRAW  (all tiers × ±2°, ±5°, ±10°)")
    print("  Pool: all full-corpus sites within ±W° of any dome site longitude.")
    print("  Tests: do dome sites outperform other monuments from the same footprint?")
    print("═" * 80)

    PRIMARY_W = 5.0
    null_c_results = {}   # keyed by (tlabel, w)
    pools = {}            # keyed by w

    for w in [2.0, 5.0, 10.0]:
        dists = np.abs(all_lons[:, None] - dome_lons[None, :])
        pool  = all_lons[dists.min(axis=1) <= w]
        pools[w] = pool
        mat_w = rng.choice(pool, size=(N_PERMS, N_dome), replace=True)

        print(f"\n  Window ±{w:.0f}°  (pool N = {len(pool)})")
        print(f"  {'Tier':<6}  {'Obs':>4}  {'Pool%':>6}  {'Mean':>6}  "
              f"{'SD':>5}  {'Z':>6}  {'p':>9}  {'Sig':>4}")
        print("  " + "─" * 60)

        for tlabel, thresh, p_null_rate in TIERS:
            obs        = count_tier(dome_lons, thresh)
            pool_rate  = count_tier(pool, thresh) / len(pool) * 100
            boot       = boot_counts(mat_w, thresh)
            p_val      = float(np.mean(boot >= obs))
            mean       = float(boot.mean())
            std        = float(boot.std())
            z          = (obs - mean) / std if std > 0 else float("nan")
            null_c_results[(tlabel, w)] = dict(
                obs=obs, pool_rate=pool_rate, mean=mean, std=std, z=z, p=p_val,
                n_pool=len(pool)
            )
            print(f"  {tlabel:<6}  {obs:>4}        {pool_rate:>5.1f}%  {mean:>6.2f}  "
                  f"{std:>5.2f}  {z:>6.2f}  {p_val:>9.4f}  {_sig(p_val):>4}")

    # ── Summary table — primary tier A++ ──────────────────────────────────────
    print("\n" + "═" * 80)
    print("  SUMMARY  (primary tier A++)")
    print("═" * 80)
    ra = null_a_results["A++"]
    rc = null_c_results[("A++", PRIMARY_W)]
    print(f"  {'Null':<52}  {'Mean':>6}  {'Z':>6}  {'p':>9}  {'Sig':>4}")
    print(f"  {'─'*52}  {'─'*6}  {'─'*6}  {'─'*9}  {'─'*4}")
    print(f"  {'A: within-dome bootstrap':<52}  {ra['mean']:>6.2f}  {ra['z']:>6.2f}  "
          f"{ra['p']:>9.4f}  {_sig(ra['p']):>4}")
    c_label = f"C: restricted pool ±{PRIMARY_W:.0f}° (N={rc['n_pool']})"
    print(f"  {c_label:<52}  {rc['mean']:>6.2f}  {rc['z']:>6.2f}  "
          f"{rc['p']:>9.4f}  {_sig(rc['p']):>4}")
    print(f"\n  Observed A++ = {ra['obs']}/{N_dome}  "
          f"(geometric null {100*P_NULL_APP:.0f}%, expected {N_dome*P_NULL_APP:.2f})")

    # ── LaTeX macros (GROUP 26) ───────────────────────────────────────────────
    # Primary macros use A++ (most stringent tier) for Null A and Null C ±5°
    print("\n" + "═" * 80)
    print("  LATEX MACROS (GROUP 26 — geographic-concentration null, primary = A++)")
    print("═" * 80)

    ra_app = null_a_results["A++"]
    rc_app_5 = null_c_results[("A++", 5.0)]
    _tag_map = {2.0: "two", 5.0: "five", 10.0: "ten"}

    # Null A macros (primary = A++)
    for k, v in [
        ("geoNullDomeBootP",    ra_app["p"]),
        ("geoNullDomeBootZ",    ra_app["z"]),
        ("geoNullDomeBootMean", ra_app["mean"]),
    ]:
        print(f"  \\newcommand{{\\{k}}}{{{v:.4f}}}")

    # Null C primary (±5°, A++)
    for k, v in [
        ("geoNullDomeRestrictedP",    rc_app_5["p"]),
        ("geoNullDomeRestrictedZ",    rc_app_5["z"]),
        ("geoNullDomeRestrictedMean", rc_app_5["mean"]),
        ("geoNullDomeRestrictedN",    float(rc_app_5["n_pool"])),
        ("geoNullDomeRestrictedRate", rc_app_5["pool_rate"]),
    ]:
        print(f"  \\newcommand{{\\{k}}}{{{v:.4f}}}")

    # Null C sensitivity sweep macros (A++)
    for w, r in [(w, null_c_results[("A++", w)]) for w in [2.0, 5.0, 10.0]]:
        tag = _tag_map[w]
        print(f"  \\newcommand{{\\geoNullDomeRestrictedN{tag}}}{{{r['n_pool']}}}")
        print(f"  \\newcommand{{\\geoNullDomeRestrictedP{tag}}}{{{r['p']:.4f}}}")
        print(f"  \\newcommand{{\\geoNullDomeRestrictedMean{tag}}}{{{r['mean']:.4f}}}")

    # Corridor diagnostic
    print(f"  \\newcommand{{\\domeEurasianFraction}}{{{eurasian:.0f}}}")
    print(f"  \\newcommand{{\\fullCorpusEurasianFraction}}{{{full_e:.0f}}}")

    # ── Write to results store ────────────────────────────────────────────────
    store_data = {
        "geoNullDomeBootP":           ra_app["p"],
        "geoNullDomeBootZ":           ra_app["z"],
        "geoNullDomeBootMean":        ra_app["mean"],
        "geoNullDomeRestrictedP":     rc_app_5["p"],
        "geoNullDomeRestrictedZ":     rc_app_5["z"],
        "geoNullDomeRestrictedMean":  rc_app_5["mean"],
        "geoNullDomeRestrictedN":     float(rc_app_5["n_pool"]),
        "domeEurasianFraction":       eurasian,
        "fullCorpusEurasianFraction": full_e,
    }
    for w in [2.0, 5.0, 10.0]:
        tag = _tag_map[w]
        r   = null_c_results[("A++", w)]
        store_data[f"geoNullDomeRestrictedN{tag}"]    = float(r["n_pool"])
        store_data[f"geoNullDomeRestrictedP{tag}"]    = r["p"]
        store_data[f"geoNullDomeRestrictedMean{tag}"] = r["mean"]
    ResultsStore().write_many(store_data)
    print("\nResults written to data/store/results.json")

    return null_a_results, null_c_results, dome_lons, all_lons, N_dome


if __name__ == "__main__":
    main()
