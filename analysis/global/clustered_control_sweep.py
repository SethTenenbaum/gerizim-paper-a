"""
clustered_control_sweep.py
==========================
"Random Clustered Control" experiment — dome sub-population.

Addresses the critique: "You tuned a statistic to find a peak — any clustered
dataset would produce a similar optimal frame."

Scope: dome sites (N~90), the sub-population in which the (3 deg, ~35 E)
score-maximizing frame was identified.  Running this experiment on the full
corpus is uninformative because the full-corpus global maximiser sits at a
completely different (omega, phi), confirming that the 3 deg/35 E result is
specific to the dome sub-population.

Method:
  1. Filter UNESCO corpus to dome/stupa/monument sites (same keywords used in
     dome_periodicity_audit.py and simulation_null_model.py).
  2. Fit a Gaussian KDE to dome longitudes (preserving clustering structure).
  3. Resample N_SIM synthetic dome datasets of the same size from the KDE.
  4. For each synthetic dataset run the (omega, phi) parameter sweep:
     sweep candidate spacings omega and anchor phases phi, compute the
     enrichment ratio, record T = max S(omega, phi).
  5. Compare the distribution of T_sim to T_obs (the real dome data's maximum).

Interpretation:
  If T_obs falls well above the synthetic distribution -> the dome pattern is
  unusual even given dome-site clustering; the result is not a clustering
  artifact.
  If T_obs is typical of the synthetic distribution -> the statistic can
  manufacture an "optimal frame" from clustering alone; strong caveat needed.

Performance notes (vectorized):
  Inner loop over phases is fully vectorized using numpy broadcasting.
  Loop over omega values is sequential to bound memory to O(N x M) per step.

Run from repo root:
    python3 analysis/global/clustered_control_sweep.py

Outputs to ResultsStore:
    clusterCtrlPermP        p-value: P(T_sim >= T_obs)
    clusterCtrlTobs         observed maximum enrichment ratio (dome corpus)
    clusterCtrlTsimMean     mean T_sim across simulations
    clusterCtrlTsimPct95    95th percentile of T_sim
    clusterCtrlOmegaMatchPct  % of sims where best omega is within 0.5 deg of 3 deg
    clusterCtrlPhaseMatchPct  % of sims where best phi is within 2 deg of 35 E
"""

import sys
import time
import numpy as np
from pathlib import Path
from scipy.stats import gaussian_kde

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from data.unesco_corpus import load_corpus
from lib.dome_filter import FORM_KEYWORDS, FORM_KEYWORD_RES
from lib.beru import GERIZIM, TIER_APLUS_DEG
from lib.results_store import ResultsStore

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

N_SIM        = 2_000          # synthetic datasets
SEED         = 42
KDE_BW       = 0.05           # KDE bandwidth (same as simulation_null_model.py)
THRESHOLD    = TIER_APLUS_DEG # 0.15 deg — A+ absolute threshold (fixed across omega)

# omega sweep: candidate spacings in degrees
OMEGA_VALUES = np.concatenate([
    np.arange(0.5, 3.0, 0.25),    # fine grid near small spacings
    np.array([3.0]),               # exact 3 deg
    np.arange(3.25, 15.5, 0.25),  # above 3 deg
])

# phi sweep: candidate anchor phases (1 deg resolution, 0-360 deg)
PHASE_VALUES = np.arange(0.0, 360.0, 1.0)   # 360 phases

OMEGA_MATCH_TOL = 0.5    # omega counts as "near 3 deg" if within this many degrees
PHASE_MATCH_TOL = 2.0    # phi counts as "near 35 E" if within this many degrees


# ---------------------------------------------------------------------------
# Score function — fully vectorized over phases
# ---------------------------------------------------------------------------

def enrichment_by_phase(lons: np.ndarray, omega: float,
                         phases: np.ndarray) -> np.ndarray:
    """
    For one (omega, all phases) pair, return the enrichment ratio per phase.

    Parameters
    ----------
    lons   : (N,) array of site longitudes
    omega  : spacing in degrees
    phases : (M,) array of candidate anchor phases

    Returns
    -------
    enrichment : (M,) enrichment ratios — observed_count / expected_count
                 expected_count = N * 2 * THRESHOLD / omega
    """
    N = len(lons)
    # arc[i, j] = shortest circular distance from lons[i] to phases[j]
    arc = np.abs(lons[:, None] - phases[None, :])        # (N, M)
    arc = np.minimum(arc, 360.0 - arc)                   # wrap to [0, 180]

    # deviation from nearest harmonic node at this spacing
    dev = arc % omega                                      # (N, M)
    dev = np.minimum(dev, omega - dev)                    # fold to [0, omega/2]

    count = (dev <= THRESHOLD).sum(axis=0).astype(float)  # (M,)

    # expected count under uniform null (null rate = 2*threshold/omega)
    expected = N * 2.0 * THRESHOLD / omega                # scalar
    return count / expected if expected > 0 else count


def sweep_best(lons: np.ndarray,
               omega_values: np.ndarray,
               phases: np.ndarray) -> tuple[float, float, float]:
    """
    Run the full (omega, phi) sweep and return (best_enrichment, best_omega, best_phase).

    Memory per step: O(N x M) — bounded and manageable.
    """
    best_enrich = -1.0
    best_omega  = float("nan")
    best_phase  = float("nan")

    for omega in omega_values:
        enrich = enrichment_by_phase(lons, omega, phases)   # (M,)
        idx    = int(np.argmax(enrich))
        if enrich[idx] > best_enrich:
            best_enrich = float(enrich[idx])
            best_omega  = float(omega)
            best_phase  = float(phases[idx])

    return best_enrich, best_omega, best_phase


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("  CLUSTERED CONTROL SWEEP  (dome sub-population)")
    print(f"  N_SIM={N_SIM:,}  |  omega candidates={len(OMEGA_VALUES)}  "
          f"|  phi candidates={len(PHASE_VALUES)}")
    print("=" * 80)

    # -- Load corpus and filter to dome sites --------------------------------
    corpus = load_corpus()
    dome_lons = []
    for site in corpus:
        if site.category == "Natural" or not site.has_coords:
            continue
        ft = site.full_text
        if any(FORM_KEYWORD_RES[kw].search(ft) for kw in FORM_KEYWORDS):
            dome_lons.append(site.longitude)

    lons = np.array(dome_lons, dtype=float)
    N    = len(lons)
    print(f"\n  Dome sub-population: N = {N} sites")
    print(f"  Longitude range: {lons.min():.1f} to {lons.max():.1f} deg")

    # -- Observed sweep ------------------------------------------------------
    print("\n  Running sweep on real dome data...", end=" ", flush=True)
    t0 = time.perf_counter()
    T_obs, obs_omega, obs_phase = sweep_best(lons, OMEGA_VALUES, PHASE_VALUES)
    print(f"done ({time.perf_counter()-t0:.1f}s)")
    print(f"  T_obs = {T_obs:.4f}  (best omega={obs_omega:.2f} deg, best phi={obs_phase:.1f} deg)")

    # Also report enrichment at the fixed (3 deg, Gerizim) target frame
    enrich_at_target = enrichment_by_phase(lons, 3.0, PHASE_VALUES)
    target_phase_idx = int(np.argmin(np.abs(PHASE_VALUES - GERIZIM)))
    S_target = float(enrich_at_target[target_phase_idx])
    print(f"  S at (3 deg, {GERIZIM:.1f} E) = {S_target:.4f}  "
          f"[target frame; T_obs is the global max over all (omega, phi)]")

    # -- Fit KDE to dome longitudes ------------------------------------------
    print("\n  Fitting KDE to dome longitudes...", end=" ", flush=True)
    kde = gaussian_kde(lons, bw_method=KDE_BW)
    print("done")

    # -- Simulate ------------------------------------------------------------
    rng = np.random.default_rng(SEED)
    T_sim       = np.empty(N_SIM)
    omega_sim   = np.empty(N_SIM)
    phase_sim   = np.empty(N_SIM)

    print(f"\n  Running {N_SIM:,} synthetic sweeps...", end=" ", flush=True)
    t0 = time.perf_counter()

    for k in range(N_SIM):
        lons_k = kde.resample(N, seed=rng).flatten()
        T_sim[k], omega_sim[k], phase_sim[k] = sweep_best(
            lons_k, OMEGA_VALUES, PHASE_VALUES
        )
        if (k + 1) % 500 == 0:
            elapsed = time.perf_counter() - t0
            rate = (k + 1) / elapsed
            print(f"\n    [{k+1:,}/{N_SIM:,}] {elapsed:.0f}s elapsed, "
                  f"{rate:.0f} sims/s, ETA ~{(N_SIM-k-1)/rate:.0f}s",
                  end="", flush=True)

    elapsed = time.perf_counter() - t0
    print(f"\n  Done — {elapsed:.1f}s total")

    # -- Results -------------------------------------------------------------
    p_val = float(np.mean(T_sim >= T_obs))
    t_mean = float(T_sim.mean())
    t_p95  = float(np.percentile(T_sim, 95))
    t_p99  = float(np.percentile(T_sim, 99))

    omega_match_pct = float(np.mean(np.abs(omega_sim - 3.0) <= OMEGA_MATCH_TOL) * 100)
    phase_match_pct = float(np.mean(np.abs(phase_sim - GERIZIM) <= PHASE_MATCH_TOL) * 100)

    print("\n" + "-" * 80)
    print("  RESULTS")
    print("-" * 80)
    print(f"  T_obs             = {T_obs:.4f}")
    print(f"  T_sim mean        = {t_mean:.4f}")
    print(f"  T_sim 95th pct    = {t_p95:.4f}")
    print(f"  T_sim 99th pct    = {t_p99:.4f}")
    print(f"  P(T_sim >= T_obs) = {p_val:.4f}   (clusterCtrlPermP)")
    print(f"\n  Observed sweep maximiser:    omega={obs_omega:.2f} deg,  phi={obs_phase:.1f} deg")
    print(f"  Sim: omega near 3 deg (+-{OMEGA_MATCH_TOL} deg):  {omega_match_pct:.1f}% of sims")
    print(f"  Sim: phi near 35 E (+-{PHASE_MATCH_TOL} deg): {phase_match_pct:.1f}% of sims")
    print("-" * 80)

    # -- Write to ResultsStore -----------------------------------------------
    store = ResultsStore()
    store.write_many({
        "clusterCtrlPermP":         round(p_val, 4),
        "clusterCtrlTobs":          round(T_obs, 4),
        "clusterCtrlTsimMean":      round(t_mean, 4),
        "clusterCtrlTsimPct95":     round(t_p95, 4),
        "clusterCtrlNsim":          N_SIM,
        "clusterCtrlOmegaMatchPct": round(omega_match_pct, 1),
        "clusterCtrlPhaseMatchPct": round(phase_match_pct, 1),
        "clusterCtrlObsOmega":      round(obs_omega, 2),
        "clusterCtrlObsPhase":      round(obs_phase, 1),
        "clusterCtrlSatTarget":     round(S_target, 4),
    })
    print(f"\n  Wrote 10 macros to ResultsStore.")

    # -- Emit LaTeX \newcommand lines (captured by reproduce_all_macros.sh) ---
    print(f"  \\newcommand{{\\clusterCtrlPermP}}{{{round(p_val, 4)}}}        % permutation p-value (clustered control sweep)")
    print(f"  \\newcommand{{\\clusterCtrlNsim}}{{{N_SIM}}}         % number of simulations (clustered control sweep)")
    print(f"  \\newcommand{{\\clusterCtrlTobs}}{{{round(T_obs, 4)}}}        % observed maximum enrichment ratio")
    print(f"  \\newcommand{{\\clusterCtrlTsimMean}}{{{round(t_mean, 4)}}}        % mean T_sim across simulations")
    print(f"  \\newcommand{{\\clusterCtrlTsimPctNF}}{{{round(t_p95, 4)}}}        % 95th percentile of T_sim")
    print(f"  \\newcommand{{\\clusterCtrlOmegaMatchPct}}{{{round(omega_match_pct, 1)}}}        % pct of sims where best omega is near 3 deg")
    print(f"  \\newcommand{{\\clusterCtrlPhaseMatchPct}}{{{round(phase_match_pct, 1)}}}        % pct of sims where best phi is near 35 E")
    print(f"  \\newcommand{{\\clusterCtrlObsOmega}}{{{round(obs_omega, 2)}}}        % observed best-fit period (deg)")
    print(f"  \\newcommand{{\\clusterCtrlObsPhase}}{{{round(obs_phase, 1)}}}        % observed best-fit anchor longitude (deg)")
    print(f"  \\newcommand{{\\clusterCtrlSatTarget}}{{{round(S_target, 4)}}}        % saturation target S")


if __name__ == "__main__":
    main()
