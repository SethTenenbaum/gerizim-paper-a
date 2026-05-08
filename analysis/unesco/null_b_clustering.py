"""
null_b_clustering.py
====================
DEPRECATED — not used in the manuscript.

Comparing dome-site longitude σ and mean NN distance against any placement
null (global uniform or corpus-conditioned) operates at the ~100° geographic
scale, which has no projection onto the 3° harmonic period. The 3° grid tiles
uniformly over Eurasia, so a site cluster spanning 20° crosses ~7 full periods
and each member independently faces the same geometric null rate p₀ = 2d/T
regardless of where in Eurasia it sits. Null A (within-dome bootstrap) already
fully captures this: it resamples from actual dome positions and returns ns,
meaning the hit count is fully consistent with dome-site geography at the 3°
scale. A clustering null at the 100° scale adds no information. Retained for
reference only.
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

MACRO_OUT = _ROOT / "analysis" / "unesco" / "null_b_clustering_macros.tex"


# ---------------------------------------------------------------------------
# Clustering statistics
# ---------------------------------------------------------------------------

def lon_std(lons: np.ndarray) -> float:
    return float(np.std(lons, ddof=0))


def mean_nn_distance(lons: np.ndarray) -> float:
    """Mean 1-D nearest-neighbour longitude distance (scalar)."""
    s = np.sort(lons)
    d = np.diff(s)
    d_right = np.append(d, np.inf)
    d_left  = np.append(np.inf, d)
    nn = np.minimum(d_right, d_left)
    return float(np.mean(nn[np.isfinite(nn)]))


def batch_std(draws: np.ndarray) -> np.ndarray:
    """(N_PERMS, N) -> (N_PERMS,) std devs."""
    return draws.std(axis=1, ddof=0)


def batch_mean_nn(draws: np.ndarray) -> np.ndarray:
    """(N_PERMS, N) -> (N_PERMS,) mean NN distances. Vectorised."""
    s = np.sort(draws, axis=1)
    d = np.diff(s, axis=1)
    inf_c = np.full((draws.shape[0], 1), np.inf)
    dr = np.concatenate([d, inf_c], axis=1)
    dl = np.concatenate([inf_c, d], axis=1)
    nn = np.minimum(dr, dl)
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


def run_null(draws: np.ndarray, obs_std: float, obs_nn: float,
             label: str) -> dict:
    null_stds = batch_std(draws)
    null_nns  = batch_mean_nn(draws)

    p_std = float(np.mean(null_stds <= obs_std))
    p_nn  = float(np.mean(null_nns  <= obs_nn))
    z_std = float((obs_std - null_stds.mean()) / null_stds.std())
    z_nn  = float((obs_nn  - null_nns.mean())  / null_nns.std())

    W = 80
    print(f"\n  {label}")
    print("  " + "-" * W)
    print(f"  {'Statistic':<28}  {'Observed':>9}  {'Null mean':>9}  {'Null SD':>8}"
          f"  {'Z':>7}  {'p':>8}  {'Sig':>4}")
    print("  " + "-" * W)
    print(f"  {'Longitude sigma (deg)':<28}  {obs_std:>9.2f}  "
          f"{null_stds.mean():>9.2f}  {null_stds.std():>8.2f}  "
          f"{z_std:>7.2f}  {fmt_p(p_std):>8}  {sig(p_std):>4}")
    print(f"  {'Mean NN distance (deg)':<28}  {obs_nn:>9.2f}  "
          f"{null_nns.mean():>9.2f}  {null_nns.std():>8.2f}  "
          f"{z_nn:>7.2f}  {fmt_p(p_nn):>8}  {sig(p_nn):>4}")

    return dict(
        null_mean_std=float(null_stds.mean()), null_sd_std=float(null_stds.std()),
        null_mean_nn=float(null_nns.mean()),   null_sd_nn=float(null_nns.std()),
        p_std=p_std, p_nn=p_nn, z_std=z_std, z_nn=z_nn,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    corpus    = load_corpus()
    cultural  = cultural_sites_with_coords(corpus)
    all_lons  = np.array([s.longitude for s in cultural])
    dome_lons = np.array([s.longitude for s in cultural if is_dome_site_raw(s)])
    N        = len(dome_lons)
    N_corpus = len(all_lons)

    obs_std = lon_std(dome_lons)
    obs_nn  = mean_nn_distance(dome_lons)

    rng = np.random.default_rng(SEED)

    print()
    print("=" * 80)
    print("  NULL B: GEOGRAPHIC CLUSTERING NULL")
    print(f"  Dome N = {N}  |  full corpus N = {N_corpus}  |  N_PERMS = {N_PERMS:,}  |  seed = {SEED}")
    print(f"  Observed dome sigma  = {obs_std:.2f}deg   (full corpus sigma = {all_lons.std():.2f}deg)")
    print(f"  Observed dome mean NN = {obs_nn:.2f}deg  (full corpus mean NN = {mean_nn_distance(all_lons):.2f}deg)")
    print("  p-values LEFT-TAILED: p = P(null <= observed), lower value = more clustered.")
    print("=" * 80)

    # Null B1: global uniform
    draws_b1 = rng.uniform(LON_MIN, LON_MAX, size=(N_PERMS, N))
    r_b1 = run_null(draws_b1, obs_std, obs_nn,
                    "NULL B1 -- GLOBAL UNIFORM  draw: Uniform[-180, +180]")

    # Null B2: corpus-conditioned  <-- PRIMARY
    draws_b2 = rng.choice(all_lons, size=(N_PERMS, N), replace=True)
    r_b2 = run_null(draws_b2, obs_std, obs_nn,
                    f"NULL B2 -- CORPUS-CONDITIONED (PRIMARY)  draw: resample {N_corpus} UNESCO sites")

    print()
    print("=" * 80)
    print("  INTERPRETATION")
    print("=" * 80)
    print(f"  B1 (uniform, easy bar): dome sites far more clustered than {N} random global")
    print(f"    sites (sigma Z={r_b1['z_std']:.2f}, p={fmt_p(r_b1['p_std'])}, {sig(r_b1['p_std'])}).")
    print(f"  B2 (corpus-conditioned, conservative): controlling for UNESCO Eurasian bias,")
    print(f"    dome sites are {'still more' if r_b2['p_std'] < 0.10 else 'not significantly more'}"
          f" clustered than a random same-sized draw from the full corpus")
    print(f"    (sigma: Z={r_b2['z_std']:.2f}, p={fmt_p(r_b2['p_std'])}, {sig(r_b2['p_std'])};  "
          f"NN: Z={r_b2['z_nn']:.2f}, p={fmt_p(r_b2['p_nn'])}, {sig(r_b2['p_nn'])}).")

    # ---- LaTeX macros (GROUP 27) ------------------------------------------
    macros: list[tuple[str, str]] = [
        # Observed
        ("clusterNullObsSigma",      f"{obs_std:.2f}"),
        ("clusterNullObsNN",         f"{obs_nn:.2f}"),
        ("clusterNullN",             str(N)),
        ("clusterNullNperms",        str(N_PERMS)),
        ("clusterNullCorpusN",       str(N_corpus)),
        # B1 uniform
        ("clusterNullB1MeanSigma",   f"{r_b1['null_mean_std']:.1f}"),
        ("clusterNullB1MeanNN",      f"{r_b1['null_mean_nn']:.2f}"),
        ("clusterNullB1ZSigma",      f"{r_b1['z_std']:.2f}"),
        ("clusterNullB1ZNN",         f"{r_b1['z_nn']:.2f}"),
        ("clusterNullB1PSigma",      fmt_p(r_b1['p_std'])),
        ("clusterNullB1PNN",         fmt_p(r_b1['p_nn'])),
        ("clusterNullB1SigSigma",    sig(r_b1['p_std'])),
        ("clusterNullB1SigNN",       sig(r_b1['p_nn'])),
        # B2 corpus-conditioned (PRIMARY)
        ("clusterNullB2MeanSigma",   f"{r_b2['null_mean_std']:.2f}"),
        ("clusterNullB2MeanNN",      f"{r_b2['null_mean_nn']:.2f}"),
        ("clusterNullB2ZSigma",      f"{r_b2['z_std']:.2f}"),
        ("clusterNullB2ZNN",         f"{r_b2['z_nn']:.2f}"),
        ("clusterNullB2PSigma",      fmt_p(r_b2['p_std'])),
        ("clusterNullB2PNN",         fmt_p(r_b2['p_nn'])),
        ("clusterNullB2SigSigma",    sig(r_b2['p_std'])),
        ("clusterNullB2SigNN",       sig(r_b2['p_nn'])),
        # backward-compat aliases pointing at B2 (primary)
        ("clusterNullZSigma",        f"{r_b2['z_std']:.2f}"),
        ("clusterNullZNN",           f"{r_b2['z_nn']:.2f}"),
        ("clusterNullPSigma",        fmt_p(r_b2['p_std'])),
        ("clusterNullPNN",           fmt_p(r_b2['p_nn'])),
        ("clusterNullSigSigma",      sig(r_b2['p_std'])),
        ("clusterNullSigNN",         sig(r_b2['p_nn'])),
        ("clusterNullNullMeanSigma", f"{r_b2['null_mean_std']:.2f}"),
        ("clusterNullNullMeanNN",    f"{r_b2['null_mean_nn']:.2f}"),
    ]

    print()
    print("=" * 80)
    print("  LATEX MACROS (GROUP 27)")
    print("=" * 80)
    macro_lines = []
    for name, val in macros:
        line = f"\\newcommand{{\\{name}}}{{{val}}}"
        print(f"  {line}")
        macro_lines.append(line)

    MACRO_OUT.write_text(
        "% Null B clustering macros -- auto-generated by null_b_clustering.py\n"
        + "\n".join(macro_lines) + "\n"
    )
    print(f"\n  Local macros -> {MACRO_OUT.relative_to(_ROOT)}")

    # ---- Results store -------------------------------------------------------
    ResultsStore().write_many({
        "clusterNullObsSigma":  obs_std,
        "clusterNullObsNN":     obs_nn,
        "clusterNullN":         float(N),
        "clusterNullB1PSigma":  r_b1["p_std"],  "clusterNullB1PNN": r_b1["p_nn"],
        "clusterNullB1ZSigma":  r_b1["z_std"],  "clusterNullB1ZNN": r_b1["z_nn"],
        "clusterNullB2PSigma":  r_b2["p_std"],  "clusterNullB2PNN": r_b2["p_nn"],
        "clusterNullB2ZSigma":  r_b2["z_std"],  "clusterNullB2ZNN": r_b2["z_nn"],
    })
    print("  Results -> data/store/results.json")


if __name__ == "__main__":
    main()
