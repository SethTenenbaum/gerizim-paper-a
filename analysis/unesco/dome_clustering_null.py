"""
dome_clustering_null.py
=======================
Geographic-clustering null for the UNESCO dome/stupa corpus.

MOTIVATION
----------
Null A (within-dome bootstrap) tests whether harmonic-proximity hit counts
are reproducible by resampling from dome-site longitudes — it is non-
significant, localising the signal to where domed monuments are built.

Null B (uniform draw) tests whether the *hit count* exceeds what a
globally uniform placement produces — it is significant, confirming the
observed proximity enrichment is globally atypical.

Neither null directly answers the underlying spatial question:

    "Is the geographic concentration of dome/stupa sites itself improbable
     under a null model of uniform random placement?  Or is the observed
     degree of clustering typical of N sites drawn uniformly from [-180°, 180°]?"

This script (Null D) answers that question by measuring the geographic
*clustering* of dome sites via two statistics:

  (1) Longitude standard deviation (σ).
      Under uniform[-180, 180], E[σ] ≈ 103.9°.
      A smaller σ indicates geographic concentration.
      p = P(null σ ≤ observed σ) — left-tail, smaller = more clustered.

  (2) Mean nearest-neighbour longitude distance (d̄_NN).
      For each site, find its closest neighbour in longitude.
      Lower d̄_NN = tighter clustering.
      p = P(null d̄_NN ≤ observed d̄_NN) — left-tail.

Both statistics are computed for N_PERMS draws of N_dome longitudes
from Uniform[-180°, 180°].  Empirical p-values are one-sided (left-tail),
measuring how rarely a uniform draw produces at least as much clustering
as observed.

Run from repo root:
    python3 analysis/unesco/dome_clustering_null.py

Emits LaTeX macros and writes to ResultsStore (GROUP 27).
~15 seconds at N_PERMS = 100_000.
"""

from __future__ import annotations
import sys
import json
import numpy as np
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT))

from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.dome_filter import is_dome_site_raw
from lib.results_store import ResultsStore

_CFG    = json.loads((_ROOT / "config.json").read_text())
N_PERMS = _CFG["simulation"]["n_permutations"]
SEED    = _CFG["simulation"]["random_seed"]

LON_MIN, LON_MAX = -180.0, 180.0

MACRO_OUT = _ROOT / "analysis" / "unesco" / "dome_clustering_null_macros.tex"


# ── Clustering statistics ─────────────────────────────────────────────────────

def lon_std(lons: np.ndarray) -> float:
    """Longitude standard deviation (scalar)."""
    return float(np.std(lons, ddof=0))


def mean_nn_distance(lons: np.ndarray) -> float:
    """Mean 1-D nearest-neighbour longitude distance."""
    s = np.sort(lons)
    # For each point, NN is either the previous or next in sorted order.
    # Edge cases: first and last only have one neighbour in sorted list.
    d_right = np.diff(s)   # distance to right neighbour (length N-1)
    d_right_full = np.append(d_right, np.inf)   # last point: no right
    d_left_full  = np.append(np.inf, d_right)   # first point: no left
    nn = np.minimum(d_right_full, d_left_full)
    return float(np.mean(nn[np.isfinite(nn)]))


def batch_std(draws: np.ndarray) -> np.ndarray:
    """(N_PERMS, N) → (N_PERMS,) standard deviations."""
    return draws.std(axis=1, ddof=0)


def batch_mean_nn(draws: np.ndarray) -> np.ndarray:
    """(N_PERMS, N) → (N_PERMS,) mean NN distances. Vectorised."""
    sorted_draws = np.sort(draws, axis=1)          # (N_PERMS, N)
    d_right = np.diff(sorted_draws, axis=1)        # (N_PERMS, N-1)
    # Pad with inf on both sides for edge cases
    inf_col = np.full((draws.shape[0], 1), np.inf)
    d_right_full = np.concatenate([d_right, inf_col], axis=1)   # (N_PERMS, N)
    d_left_full  = np.concatenate([inf_col, d_right], axis=1)   # (N_PERMS, N)
    nn = np.minimum(d_right_full, d_left_full)                   # (N_PERMS, N)
    # Replace inf (only at N=1, never in practice) with nan
    nn[~np.isfinite(nn)] = np.nan
    return np.nanmean(nn, axis=1)


def sig(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"


def fmt_p(p: float) -> str:
    if p >= 0.10:  return f"{p:.2f}"
    if p >= 0.001: return f"{p:.3f}"
    return f"{p:.4f}"


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    corpus   = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    dome_lons = np.array([s.longitude for s in cultural if is_dome_site_raw(s)])
    N = len(dome_lons)

    rng = np.random.default_rng(SEED)

    # Observed clustering statistics
    obs_std  = lon_std(dome_lons)
    obs_nn   = mean_nn_distance(dome_lons)

    print()
    print("=" * 80)
    print("  NULL D: GEOGRAPHIC CLUSTERING NULL  (uniform random placement)")
    print(f"  UNESCO dome/stupa corpus  N = {N}")
    print(f"  Null draw: Uniform[{LON_MIN}°, {LON_MAX}°]  |  N_PERMS = {N_PERMS:,}  |  seed = {SEED}")
    print(f"  Under Uniform[-180,180]:  E[σ] ≈ {LON_MAX/np.sqrt(3):.1f}°")
    print(f"  Observed longitude σ  = {obs_std:.2f}°  (lower = more clustered)")
    print(f"  Observed mean NN dist = {obs_nn:.2f}°  (lower = more clustered)")
    print("=" * 80)
    print()
    print("  Both p-values are LEFT-TAILED:")
    print("  p = P(null statistic ≤ observed)  — probability that a uniform")
    print("  random draw is AT LEAST as clustered as the observed corpus.")
    print()

    # Draw null samples
    draws = rng.uniform(LON_MIN, LON_MAX, size=(N_PERMS, N))

    null_stds = batch_std(draws)
    null_nns  = batch_mean_nn(draws)

    # Left-tail p-values (how often is the null as clustered as observed?)
    p_std = float(np.mean(null_stds <= obs_std))
    p_nn  = float(np.mean(null_nns  <= obs_nn))

    z_std = (obs_std - null_stds.mean()) / null_stds.std()
    z_nn  = (obs_nn  - null_nns.mean())  / null_nns.std()

    print(f"  {'Statistic':<30}  {'Observed':>10}  {'Null mean':>10}  {'Null SD':>8}  "
          f"{'Z':>7}  {'p (left)':>10}  {'Sig':>4}")
    print("  " + "─" * 82)
    print(f"  {'Longitude σ (°)':<30}  {obs_std:>10.2f}  "
          f"{null_stds.mean():>10.2f}  {null_stds.std():>8.2f}  "
          f"{z_std:>7.2f}  {p_std:>10.4f}  {sig(p_std):>4}")
    print(f"  {'Mean NN distance (°)':<30}  {obs_nn:>10.2f}  "
          f"{null_nns.mean():>10.2f}  {null_nns.std():>8.2f}  "
          f"{z_nn:>7.2f}  {p_nn:>10.4f}  {sig(p_nn):>4}")

    print()
    print("  INTERPRETATION")
    print(f"  Observed σ = {obs_std:.2f}° vs uniform null mean {null_stds.mean():.1f}° "
          f"(Z = {z_std:.2f}, p = {fmt_p(p_std)}, {sig(p_std)})")
    print(f"  Observed mean NN = {obs_nn:.2f}° vs uniform null mean {null_nns.mean():.1f}° "
          f"(Z = {z_nn:.2f}, p = {fmt_p(p_nn)}, {sig(p_nn)})")

    # ── LaTeX macros (GROUP 27) ───────────────────────────────────────────────
    macros = [
        ("clusterNullObsSigma",     f"{obs_std:.2f}"),
        ("clusterNullObsNN",        f"{obs_nn:.2f}"),
        ("clusterNullNullMeanSigma",f"{null_stds.mean():.1f}"),
        ("clusterNullNullMeanNN",   f"{null_nns.mean():.2f}"),
        ("clusterNullNullSDSigma",  f"{null_stds.std():.2f}"),
        ("clusterNullNullSDNN",     f"{null_nns.std():.2f}"),
        ("clusterNullZSigma",       f"{z_std:.2f}"),
        ("clusterNullZNN",          f"{z_nn:.2f}"),
        ("clusterNullPSigma",       fmt_p(p_std)),
        ("clusterNullPNN",          fmt_p(p_nn)),
        ("clusterNullSigSigma",     sig(p_std)),
        ("clusterNullSigNN",        sig(p_nn)),
        ("clusterNullN",            str(N)),
        ("clusterNullNperms",       str(N_PERMS)),
    ]

    print()
    print("=" * 80)
    print("  LATEX MACROS (GROUP 27 — geographic clustering null)")
    print("=" * 80)
    macro_lines = []
    for name, val in macros:
        line = f"\\newcommand{{\\{name}}}{{{val}}}"
        print(f"  {line}")
        macro_lines.append(line)

    MACRO_OUT.write_text("\n".join(macro_lines) + "\n")
    print(f"\n  Local macro file written to {MACRO_OUT.relative_to(_ROOT)}")

    # Write to results store
    store_data = {
        "clusterNullObsSigma":      obs_std,
        "clusterNullObsNN":         obs_nn,
        "clusterNullNullMeanSigma": float(null_stds.mean()),
        "clusterNullNullMeanNN":    float(null_nns.mean()),
        "clusterNullZSigma":        float(z_std),
        "clusterNullZNN":           float(z_nn),
        "clusterNullPSigma":        p_std,
        "clusterNullPNN":           p_nn,
        "clusterNullN":             float(N),
    }
    ResultsStore().write_many(store_data)
    print("  Results written to data/store/results.json")


if __name__ == "__main__":
    main()
