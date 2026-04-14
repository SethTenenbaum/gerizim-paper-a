"""
stupa_geographic_concentration_test.py
=======================================
Geographic-concentration null for the stupa sub-population (Test 2b stage).

MOTIVATION
----------
Table 3 of the paper (Test 2b stage breakdown) shows that stupas drive the
strongest per-stage enrichment (5.56×, N = NevoStupa, p = pEvoStupa **).
A reviewer correctly notes: are stupas enriched because they are
morphologically special, or because they cluster in South/Southeast Asian
longitude bands that coincide with the beru-harmonic structure?

This script runs the same Null A / Null C design used for the full dome
sub-population (dome_geographic_concentration_test.py) against the stupa
sub-population specifically.

  NULL A — Within-stupa bootstrap:
    Draw N_stupa longitudes with replacement from the stupa sites' own
    longitude list.  If p >= 0.05: geographic concentration of stupas
    in their actual locations explains the A+ rate.

  NULL B (new, stupa-specific) — South/Southeast Asia draw:
    Build a longitude pool from all full-corpus sites whose longitude
    falls in the 70°–110°E band (the primary stupa heartland: India,
    Nepal, Myanmar, Sri Lanka, Indonesia, Cambodia, Laos, Thailand).
    Draw N_stupa sites from this pool.  Tests: do stupa sites outperform
    OTHER monuments from their own regional heartland?

  NULL C — Restricted geographic draw:
    Same as in dome_geographic_concentration_test.py: build a pool from
    all full-corpus sites within ±5° of any stupa site longitude.
    Draw N_stupa sites from this pool.

The three-null design directly addresses the reviewer's concern:
  - Null A shows whether stupa longitudes encode the enrichment.
  - Null B tests specificity within the South/SE Asian heartland.
  - Null C tests specificity within the full geographic footprint.

STUPA KEYWORD MATCHING
----------------------
Uses the same keyword logic as tumulus_dome_evolution_raw_sweep.py:
keywords from keywords.json → mound_evolution → stupa list.
No context gating (raw sweep, consistent with the primary Test 2b result).

USAGE
-----
    cd /path/to/gerizim-paper-a
    python3 analysis/unesco/stupa_geographic_concentration_test.py

IMPORTANT: ~2 minutes at N_PERMS=100,000.  Writes results to store.

MANUSCRIPT MACROS PRODUCED
---------------------------
    \\stupaGeoBootP        — Null A: within-stupa bootstrap p
    \\stupaGeoBootZ        — Null A: Z-score
    \\stupaGeoBootMean     — Null A: bootstrap mean A+
    \\stupaRegionP         — Null B: South/SE Asia draw p
    \\stupaRegionZ         — Null B: Z-score
    \\stupaRegionMean      — Null B: pool mean A+
    \\stupaRegionN         — Null B: pool size
    \\stupaGeoRestrictedP  — Null C: restricted pool p
    \\stupaGeoRestrictedZ  — Null C: Z-score
    \\stupaGeoRestrictedMean — Null C: pool mean A+
    \\stupaGeoRestrictedN  — Null C: pool size
    \\stupaHeartlandFrac   — % stupa sites in 70°–110°E band
    \\fullCorpusHeartlandFrac — % full corpus in 70°–110°E band

GROUP: 29b (runs alongside dome_geographic_concentration_test.py)
"""

import sys
import json
import re
from pathlib import Path

import numpy as np

np.random.seed(42)

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, TIER_APLUS, P_NULL_AP
from lib.beru import deviation as _beru_deviation
from lib.results_store import ResultsStore

N_PERMS = 100_000

# Stupa heartland longitude band (South/Southeast Asia)
HEARTLAND_LO = 70.0
HEARTLAND_HI = 110.0


# ── Load stupa keyword list from keywords.json ────────────────────────────────
_KW_PATH = Path(__file__).parent.parent.parent / "keywords.json"
with open(_KW_PATH) as _f:
    _KW = json.load(_f)

STUPA_KEYWORDS = _KW["mound_evolution"]["stupa"]

STUPA_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    for kw in STUPA_KEYWORDS
}


def is_stupa_site(site_obj) -> bool:
    """Return True if any stupa keyword matches the site's full text."""
    full_text = site_obj.full_text
    return any(rx.search(full_text) for rx in STUPA_RES.values())


# ── Helpers ───────────────────────────────────────────────────────────────────

def beru_dev(lon: float) -> float:
    return _beru_deviation(lon)


def count_aplus(lons) -> int:
    """Count A+ sites. Accepts list or ndarray; always vectorized."""
    arr = lons if isinstance(lons, np.ndarray) else np.asarray(lons)
    arc = np.abs(arr - GERIZIM)
    bv  = arc / BERU
    dev = np.abs(bv - np.round(bv / 0.1) * 0.1)
    return int(np.sum(dev <= TIER_APLUS))


# ── Load corpus ───────────────────────────────────────────────────────────────

def load_all():
    corpus   = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    all_lons   = np.array([s.longitude for s in cultural])
    stupa_lons = np.array([s.longitude for s in cultural if is_stupa_site(s)])
    return all_lons, stupa_lons


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    all_lons, stupa_lons = load_all()
    N_stupa       = len(stupa_lons)
    obs_stupa_ap  = count_aplus(stupa_lons)
    exp_null      = N_stupa * P_NULL_AP
    rate_obs      = 100 * obs_stupa_ap / N_stupa

    print("=" * 80)
    print("  STUPA GEOGRAPHIC-CONCENTRATION NULL TEST")
    print(f"  Anchor: {GERIZIM}°E  |  BERU: {BERU}°  |  N_perms: {N_PERMS:,}")
    print(f"  Stupa sites (raw keyword sweep): N = {N_stupa}")
    print(f"  Keywords: {STUPA_KEYWORDS}")
    print(f"  Observed A+ = {obs_stupa_ap}  ({rate_obs:.1f}%)")
    print(f"  Geometric null: {P_NULL_AP:.0%} → expected A+ = {exp_null:.2f}")
    print("=" * 80)

    # Geographic diagnostics
    heartland_stupa = float(
        np.sum((stupa_lons >= HEARTLAND_LO) & (stupa_lons <= HEARTLAND_HI)) / N_stupa * 100
    )
    heartland_full = float(
        np.sum((all_lons >= HEARTLAND_LO) & (all_lons <= HEARTLAND_HI)) / len(all_lons) * 100
    )
    print(f"\n  Stupa sites in {HEARTLAND_LO:.0f}°–{HEARTLAND_HI:.0f}°E (S/SE Asia): "
          f"{heartland_stupa:.0f}%  (full corpus: {heartland_full:.0f}%)")
    print(f"  → Stupa sites are concentrated in the S/SE Asian heartland at "
          f"{heartland_stupa/heartland_full:.1f}× the full-corpus rate.")

    rng = np.random.default_rng(42)
    _sig = lambda p: "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "~" if p < 0.10 else "ns"

    # ── Null A: within-stupa bootstrap ───────────────────────────────────────
    print("\n" + "─" * 80)
    print("  NULL A: WITHIN-STUPA BOOTSTRAP")
    print(f"  Draw N={N_stupa} from the stupa sites' own longitudes (with replacement).")
    print("─" * 80)

    # Fully batched: draw (N_PERMS, N_stupa) matrix, compute A+ count per row
    mat_a = rng.choice(stupa_lons, size=(N_PERMS, N_stupa), replace=True)
    arc_a = np.abs(mat_a - GERIZIM) / BERU
    dev_a = np.abs(arc_a - np.round(arc_a / 0.1) * 0.1)
    boot_a = (dev_a <= TIER_APLUS).sum(axis=1).astype(int)

    p_a    = float(np.mean(boot_a >= obs_stupa_ap))
    mean_a = float(boot_a.mean())
    std_a  = float(boot_a.std())
    z_a    = (obs_stupa_ap - mean_a) / std_a if std_a > 0 else float("nan")

    print(f"\n  Null A: within-stupa bootstrap (N={N_stupa}, {N_PERMS:,} trials)")
    print(f"    Bootstrap mean A+ = {mean_a:.2f} ± {std_a:.2f}")
    print(f"    Observed A+       = {obs_stupa_ap}")
    print(f"    Z = {z_a:.2f},  p (>= {obs_stupa_ap}) = {p_a:.4f}  {_sig(p_a)}")

    # ── Null B: South/SE Asia heartland draw ─────────────────────────────────
    print("\n" + "─" * 80)
    print(f"  NULL B: SOUTH/SE ASIA HEARTLAND DRAW  ({HEARTLAND_LO:.0f}°–{HEARTLAND_HI:.0f}°E)")
    print(f"  Pool: all full-corpus sites in the stupa heartland band.")
    print("  Tests whether stupa sites are more enriched than OTHER monuments")
    print("  from the same regional heartland.")
    print("─" * 80)

    heartland_pool = all_lons[
        (all_lons >= HEARTLAND_LO) & (all_lons <= HEARTLAND_HI)
    ]
    N_pool_b      = len(heartland_pool)
    pool_b_ap_rate = count_aplus(heartland_pool) / N_pool_b * 100

    print(f"\n  Heartland pool: N = {N_pool_b}  (A+ rate in pool = {pool_b_ap_rate:.1f}%)")

    mat_b = rng.choice(heartland_pool, size=(N_PERMS, N_stupa), replace=True)
    arc_b = np.abs(mat_b - GERIZIM) / BERU
    dev_b = np.abs(arc_b - np.round(arc_b / 0.1) * 0.1)
    boot_b = (dev_b <= TIER_APLUS).sum(axis=1).astype(int)

    p_b    = float(np.mean(boot_b >= obs_stupa_ap))
    mean_b = float(boot_b.mean())
    std_b  = float(boot_b.std())
    z_b    = (obs_stupa_ap - mean_b) / std_b if std_b > 0 else float("nan")

    print(f"\n  Null B: heartland pool (N_pool={N_pool_b}, {N_PERMS:,} draws of N={N_stupa})")
    print(f"    Pool mean A+ = {mean_b:.2f} ± {std_b:.2f}")
    print(f"    Observed A+  = {obs_stupa_ap}")
    print(f"    Z = {z_b:.2f},  p (>= {obs_stupa_ap}) = {p_b:.4f}  {_sig(p_b)}")

    # ── Null C: restricted geographic draw (±5° footprint) ───────────────────
    print("\n" + "─" * 80)
    print("  NULL C: RESTRICTED GEOGRAPHIC DRAW  (±5° of any stupa site)")
    print("  Pool: all full-corpus sites within ±5° of any stupa site longitude.")
    print("  Tests whether stupa sites outperform OTHER monuments from the same")
    print("  broad longitude footprint.")
    print("─" * 80)

    # Vectorized pool construction: broadcast (N_all,) vs (N_stupa,) → (N_all, N_stupa)
    dists = np.abs(all_lons[:, None] - stupa_lons[None, :])  # (N_all, N_stupa)
    in_pool = dists.min(axis=1) <= 5.0                        # (N_all,)
    restricted_pool = all_lons[in_pool]
    N_pool_c      = len(restricted_pool)
    pool_c_ap_rate = count_aplus(restricted_pool) / N_pool_c * 100

    print(f"\n  Restricted pool: N = {N_pool_c}  (A+ rate in pool = {pool_c_ap_rate:.1f}%)")

    mat_c = rng.choice(restricted_pool, size=(N_PERMS, N_stupa), replace=True)
    arc_c = np.abs(mat_c - GERIZIM) / BERU
    dev_c = np.abs(arc_c - np.round(arc_c / 0.1) * 0.1)
    boot_c = (dev_c <= TIER_APLUS).sum(axis=1).astype(int)

    p_c    = float(np.mean(boot_c >= obs_stupa_ap))
    mean_c = float(boot_c.mean())
    std_c  = float(boot_c.std())
    z_c    = (obs_stupa_ap - mean_c) / std_c if std_c > 0 else float("nan")

    print(f"\n  Null C: restricted pool (N_pool={N_pool_c}, {N_PERMS:,} draws of N={N_stupa})")
    print(f"    Pool mean A+ = {mean_c:.2f} ± {std_c:.2f}")
    print(f"    Observed A+  = {obs_stupa_ap}")
    print(f"    Z = {z_c:.2f},  p (>= {obs_stupa_ap}) = {p_c:.4f}  {_sig(p_c)}")

    # ── Combined interpretation ───────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("  COMBINED INTERPRETATION")
    print("=" * 80)
    print(f"""
  Null A (within-stupa bootstrap, p = {p_a:.4f}, {_sig(p_a)}):
    Resampling from the stupa sites' own longitudes yields ~{mean_a:.1f} A+ on
    average (vs. observed {obs_stupa_ap}).
    {'Geographic concentration of stupa sites in their actual locations explains' if p_a >= 0.05 else 'The stupa sites are MORE enriched than their own geographic distribution predicts.'}
    {'their A+ rate.  Null B and Null C below test specificity.' if p_a >= 0.05 else ''}

  Null B (S/SE Asian heartland draw, p = {p_b:.4f}, {_sig(p_b)}):
    Drawing N={N_stupa} from the {N_pool_b}-site S/SE Asian heartland pool (A+ rate
    {pool_b_ap_rate:.1f}%) yields ~{mean_b:.1f} A+ on average.
    {'Stupa sites are MORE enriched than other monuments from the same regional heartland.' if p_b < 0.05 else 'Other monuments from the S/SE Asian heartland show similar enrichment.'}

  Null C (±5° restricted pool, p = {p_c:.4f}, {_sig(p_c)}):
    Drawing N={N_stupa} from the {N_pool_c}-site restricted pool (A+ rate
    {pool_c_ap_rate:.1f}%) yields ~{mean_c:.1f} A+ on average.
    {'Stupa sites significantly exceed the expectation from their broader longitude footprint.' if p_c < 0.05 else 'Other monuments from the same longitude footprint show comparable enrichment.'}
""")

    # ── Summary table ─────────────────────────────────────────────────────────
    print("─" * 80)
    print("  SUMMARY TABLE  (stupa sub-population geographic-concentration nulls)")
    print("─" * 80)
    print(f"  {'Null':<52}  {'mean A+':>8}  {'Z':>6}  {'p':>9}  {'Sig':>5}")
    print(f"  {'-'*52}  {'-'*8}  {'-'*6}  {'-'*9}  {'-'*5}")
    rows = [
        ("A: within-stupa bootstrap",           mean_a, z_a, p_a),
        (f"B: S/SE Asia heartland (N={N_pool_b})",  mean_b, z_b, p_b),
        (f"C: restricted ±5° pool (N={N_pool_c})",  mean_c, z_c, p_c),
    ]
    for label, mean, z, p in rows:
        print(f"  {label:<52}  {mean:>8.2f}  {z:>6.2f}  {p:>9.4f}  {_sig(p):>5}")
    print(f"\n  Observed: {obs_stupa_ap}/{N_stupa} A+ = {rate_obs:.1f}%  "
          f"(geometric null {100*P_NULL_AP:.0f}%, expected {exp_null:.2f})")

    # ── LaTeX macros ──────────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("  LATEX MACROS (GROUP 29b — stupa geographic-concentration null):")
    print("=" * 80)
    macro_pairs = [
        ("stupaGeoBootP",           f"{p_a:.4f}"),
        ("stupaGeoBootZ",           f"{z_a:.2f}"),
        ("stupaGeoBootMean",        f"{mean_a:.2f}"),
        ("stupaRegionP",            f"{p_b:.4f}"),
        ("stupaRegionZ",            f"{z_b:.2f}"),
        ("stupaRegionMean",         f"{mean_b:.2f}"),
        ("stupaRegionN",            str(N_pool_b)),
        ("stupaGeoRestrictedP",     f"{p_c:.4f}"),
        ("stupaGeoRestrictedZ",     f"{z_c:.2f}"),
        ("stupaGeoRestrictedMean",  f"{mean_c:.2f}"),
        ("stupaGeoRestrictedN",     str(N_pool_c)),
        ("stupaHeartlandFrac",      f"{heartland_stupa:.0f}"),
        ("fullCorpusHeartlandFrac", f"{heartland_full:.0f}"),
    ]
    for name, val in macro_pairs:
        print(f"\\newcommand{{\\{name}}}{{{val}}}")

    # ── Write to results store ────────────────────────────────────────────────
    ResultsStore().write_many({
        "stupaGeoBootP":           p_a,
        "stupaGeoBootZ":           z_a,
        "stupaGeoBootMean":        mean_a,
        "stupaRegionP":            p_b,
        "stupaRegionZ":            z_b,
        "stupaRegionMean":         mean_b,
        "stupaRegionN":            float(N_pool_b),
        "stupaGeoRestrictedP":     p_c,
        "stupaGeoRestrictedZ":     z_c,
        "stupaGeoRestrictedMean":  mean_c,
        "stupaGeoRestrictedN":     float(N_pool_c),
        "stupaHeartlandFrac":      heartland_stupa,
        "fullCorpusHeartlandFrac": heartland_full,
    })
    print("\nResults written to data/store/results.json")
    print("=" * 80)


if __name__ == "__main__":
    main()
