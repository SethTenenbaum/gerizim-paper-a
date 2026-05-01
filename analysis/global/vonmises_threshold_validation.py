"""
vonmises_threshold_validation.py
=================================
Fits a von Mises mixture model to the harmonic phase distributions of two
external corpora (OWTRAD Silk Road nodes and Wikidata Q180987 stupas) and
compares the inferred angular dispersion σ = κ^{-1/2} to the metrologically
derived tier thresholds.

The model:
    f(θ | μ, κ) = w · VM(μ=0, κ) + (1-w) · Uniform(0, 2π)

where θ is the harmonic phase angle, 0 = on harmonic, π = midpoint.

Results are written back to data/store/results.json under the keys:
    kappaOwtrad, kappaStupa, sigmaBeruOwtrad, sigmaBeruStupa,
    sigmaBeruPooled, vmMixtureWeightOwtrad, vmMixtureWeightStupa

Run from repo root:
    python3 analysis/global/vonmises_threshold_validation.py
"""

from __future__ import annotations

import json
import sys
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import i0

# ── Repo paths ────────────────────────────────────────────────────────────────
_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(_ROOT))

OWTRAD_CSV  = _ROOT / "data/store/silk_road/owtrad_nodes.csv"
STUPA_CSV   = _ROOT / "data/store/unesco/wikidata_stupas_q180987.csv"

# ── Constants (must match config.json / lib/beru.py) ─────────────────────────
ANCHOR  = 35.269   # degrees E
BERU    = 30.0     # degrees per beru
SPACING = 0.1      # beru per harmonic interval


# ── Angular mapping ───────────────────────────────────────────────────────────

def harmonic_deviation_beru(lon: float) -> float:
    """Absolute deviation from nearest 0.1-beru harmonic, in beru."""
    b    = (lon - ANCHOR) / BERU
    frac = b % SPACING
    if frac > SPACING / 2:
        frac = SPACING - frac
    return frac


def to_phase(lons: np.ndarray) -> np.ndarray:
    """
    Map longitudes to phase angles θ ∈ [0, π] where
      θ = 0   → exactly on harmonic
      θ = π   → midpoint between harmonics
    """
    devs  = np.array([harmonic_deviation_beru(lon) for lon in lons])
    theta = devs * (np.pi / (SPACING / 2))   # [0, π]
    return theta


# ── Von Mises mixture likelihood ──────────────────────────────────────────────

def neg_llk(params: np.ndarray, thetas: np.ndarray) -> float:
    """
    Negative log-likelihood for:
        f(θ) = w · VM(μ=0, κ) + (1-w) · Uniform(0, 2π)

    The deviation distribution is symmetric, so we use both +θ and -θ.
    """
    kappa, w = params
    if kappa <= 0 or not (0 < w < 1):
        return 1e10

    full_theta = np.concatenate([thetas, -thetas])
    log_vm     = kappa * np.cos(full_theta) - np.log(2 * np.pi * i0(kappa))
    log_unif   = -np.log(2 * np.pi)
    log_mix    = np.log(w * np.exp(log_vm) + (1 - w) * np.exp(log_unif))
    return -float(np.sum(log_mix))


def fit_mixture(lons: np.ndarray) -> tuple[float, float]:
    """
    Returns (kappa_hat, w_hat) for the best-fit von Mises mixture.
    Grid initialisation avoids local minima.
    """
    thetas = to_phase(lons)
    obj    = partial(neg_llk, thetas=thetas)

    best_nll, best_params = np.inf, None
    for k0 in [1.0, 5.0, 10.0, 20.0, 50.0, 100.0]:
        for w0 in [0.02, 0.05, 0.10, 0.20]:
            res = minimize(
                obj,
                x0     = [k0, w0],
                bounds = [(0.01, 500.0), (0.001, 0.999)],
                method = "L-BFGS-B",
            )
            if res.success and res.fun < best_nll:
                best_nll    = res.fun
                best_params = res.x

    if best_params is None:
        raise RuntimeError("von Mises mixture optimisation failed to converge.")

    return float(best_params[0]), float(best_params[1])


# ── Threshold comparison ──────────────────────────────────────────────────────

METROLOGICAL = {
    "A++": 0.00160,
    "A+":  0.00321,
    "A":   0.00642,
}


def sigma_beru(kappa: float) -> float:
    """Effective angular dispersion in beru: σ ≈ κ^{-1/2} (large-κ approximation)."""
    sigma_rad = 1.0 / np.sqrt(kappa)
    return sigma_rad * (SPACING / 2) / np.pi


def report(label: str, lons: np.ndarray) -> dict:
    kappa, w = fit_mixture(lons)
    sig      = sigma_beru(kappa)
    n        = len(lons)

    print(f"\n  {'='*58}")
    print(f"  {label}  (N={n})")
    print(f"  κ̂ = {kappa:.4f}   ŵ = {w:.4f}  ({100*w:.1f}% of mass in VM component)")
    print(f"  σ̂ = κ̂^{{-1/2}} = {sig:.5f} beru  ({sig*BERU:.4f}°)")
    print()
    print(f"  Implied thresholds (σ multiples):")
    for mult, name in [(0.5, "0.5σ"), (1.0, "1.0σ"), (2.0, "2.0σ")]:
        t = mult * sig
        print(f"    {name}  = {t:.5f} beru  ({t*BERU:.4f}°  ≈ {t*BERU*111:.1f} km)")
    print()
    print(f"  Metrological tiers for comparison:")
    for tier, t in METROLOGICAL.items():
        diff = abs(0.5 * sig - t) / t * 100
        print(f"    {tier}  = {t:.5f} beru   |  |0.5σ − {tier}| / {tier} = {diff:.1f}%")

    return {"kappa": kappa, "w": w, "sigma_beru": sig, "n": n}


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    print("\n" + "=" * 62)
    print("  VON MISES THRESHOLD VALIDATION")
    print("  Mixture:  f(θ|μ,κ) = w·VM(0,κ) + (1−w)·Uniform(0,2π)")
    print("  θ ∈ [0,π]: 0 = on-harmonic, π = midpoint")
    print("=" * 62)

    owtrad = pd.read_csv(OWTRAD_CSV, comment="#")
    stupas = pd.read_csv(STUPA_CSV,  comment="#")

    r_o = report("OWTRAD Silk Road nodes",      owtrad["lon"].values)
    r_s = report("Wikidata Q180987 stupas",      stupas["lon"].values)

    # N-weighted pooled σ
    pooled = (r_o["sigma_beru"] * r_o["n"] + r_s["sigma_beru"] * r_s["n"]) / (
        r_o["n"] + r_s["n"]
    )
    ratio  = max(r_o["sigma_beru"], r_s["sigma_beru"]) / min(
        r_o["sigma_beru"], r_s["sigma_beru"]
    )

    print(f"\n  {'='*58}")
    print(f"  Cross-corpus summary")
    print(f"  OWTRAD  κ̂ = {r_o['kappa']:.4f}   σ̂ = {r_o['sigma_beru']:.5f} beru")
    print(f"  Stupas  κ̂ = {r_s['kappa']:.4f}   σ̂ = {r_s['sigma_beru']:.5f} beru")
    print(f"  σ ratio = {ratio:.3f}  (κ ratio = {max(r_o['kappa'],r_s['kappa'])/min(r_o['kappa'],r_s['kappa']):.3f})")
    print(f"  N-weighted pooled σ̂ = {pooled:.5f} beru  ({pooled*BERU:.4f}°)")
    print()

    # ── Derived threshold values from stupa κ (primary reference corpus) ──────
    # Naming: half-sigma → A++ proxy, 1-sigma → A+ proxy, 2-sigma → A proxy
    sig_s      = r_s["sigma_beru"]
    sig_o      = r_o["sigma_beru"]
    kappa_s    = r_s["kappa"]
    kappa_o    = r_o["kappa"]
    kappa_ratio = max(kappa_o, kappa_s) / min(kappa_o, kappa_s)

    vm_app  = round(0.25 * sig_s, 5)  # 0.25σ stupa (≈ A++)
    vm_ap   = round(0.50 * sig_s, 5)  # 0.50σ stupa (≈ A+)
    vm_a    = round(1.00 * sig_s, 5)  # 1.00σ stupa (≈ A)

    # ── Write to results.json ─────────────────────────────────────────────────
    from lib.results_store import ResultsStore
    ResultsStore().write_many({
        "kappaOwtrad":             round(kappa_o,      2),
        "kappaStupa":              round(kappa_s,      2),
        "kappaRatio":              round(kappa_ratio,  2),
        "sigmaBeruOwtrad":         round(sig_o,        5),
        "sigmaBeruStupa":          round(sig_s,        5),
        "sigmaBeruPooled":         round(pooled,       5),
        "sigmaDegOwtrad":          round(sig_o * BERU, 4),
        "sigmaDegStupa":           round(sig_s * BERU, 4),
        "sigmaDegPooled":          round(pooled * BERU, 4),
        "vmMixtureWeightOwtrad":   round(r_o["w"],     4),
        "vmMixtureWeightStupa":    round(r_s["w"],     4),
        "vmHalfSigmaStupa":        vm_app,   # 0.25σ → A++ proxy (beru)
        "vmOneSigmaStupa":         vm_ap,    # 0.50σ → A+  proxy (beru)
        "vmTwoSigmaStupa":         vm_a,     # 1.00σ → A   proxy (beru)
        "vmHalfSigmaDegStupa":     round(vm_app * BERU, 4),  # 0.25σ → A++ proxy (deg)
        "vmOneSigmaDegStupa":      round(vm_ap  * BERU, 4),  # 0.50σ → A+  proxy (deg)
        "vmTwoSigmaDegStupa":      round(vm_a   * BERU, 4),  # 1.00σ → A   proxy (deg)
        "vmNOwtrad":               r_o["n"],
        "vmNStupa":                r_s["n"],
    })
    print(f"  Results written → data/store/results.json")

    # ── Emit LaTeX macros ─────────────────────────────────────────────────────
    print()
    print("  % LaTeX macros (von Mises threshold validation):")
    print(f"  \\newcommand{{\\vmKappaOwtrad}}{{{round(kappa_o, 2)}}}  % von Mises kappa, OWTRAD corpus")
    print(f"  \\newcommand{{\\vmKappaStupa}}{{{round(kappa_s, 2)}}}   % von Mises kappa, Q180987 stupa corpus")
    print(f"  \\newcommand{{\\vmKappaRatio}}{{{round(kappa_ratio, 1)}}}  % ratio max/min kappa across corpora")
    print(f"  \\newcommand{{\\vmNOwtrad}}{{{r_o['n']}}}     % N, OWTRAD corpus")
    print(f"  \\newcommand{{\\vmNStupa}}{{{r_s['n']}}}      % N, Q180987 stupa corpus")
    print(f"  \\newcommand{{\\vmSigmaBeruOwtrad}}{{{round(sig_o, 5)}}}   % sigma (beru), OWTRAD — appendix use")
    print(f"  \\newcommand{{\\vmSigmaBeruStupa}}{{{round(sig_s, 5)}}}    % sigma (beru), stupa — appendix use")
    print(f"  \\newcommand{{\\vmSigmaBeruPooled}}{{{round(pooled, 5)}}}  % N-weighted pooled sigma (beru) — appendix use")
    print(f"  \\newcommand{{\\vmSigmaDegOwtrad}}{{{round(sig_o * BERU, 4)}}}   % sigma (degrees), OWTRAD")
    print(f"  \\newcommand{{\\vmSigmaDegStupa}}{{{round(sig_s * BERU, 4)}}}    % sigma (degrees), stupa")
    print(f"  \\newcommand{{\\vmSigmaDegPooled}}{{{round(pooled * BERU, 4)}}}  % N-weighted pooled sigma (degrees)")
    print(f"  \\newcommand{{\\vmHalfSigmaStupa}}{{{vm_app}}}  % 0.25*sigma_stupa (beru) — A++ proxy — appendix use")
    print(f"  \\newcommand{{\\vmOneSigmaStupa}}{{{vm_ap}}}   % 0.50*sigma_stupa (beru) — A+  proxy — appendix use")
    print(f"  \\newcommand{{\\vmTwoSigmaStupa}}{{{vm_a}}}   % 1.00*sigma_stupa (beru) — A   proxy — appendix use")
    print(f"  \\newcommand{{\\vmHalfSigmaDegStupa}}{{{round(vm_app * BERU, 4)}}}  % 0.25*sigma_stupa (degrees) — A++ proxy")
    print(f"  \\newcommand{{\\vmOneSigmaDegStupa}}{{{round(vm_ap * BERU, 4)}}}   % 0.50*sigma_stupa (degrees) — A+ proxy")
    print(f"  \\newcommand{{\\vmTwoSigmaDegStupa}}{{{round(vm_a * BERU, 4)}}}   % 1.00*sigma_stupa (degrees) — A proxy")


if __name__ == "__main__":
    main()
