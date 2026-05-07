"""
fortification_negative_control.py
==================================
NEGATIVE ITERATIVE CONTROL — Fortification-type UNESCO sites.

PURPOSE
-------
If the 3° harmonic enrichment seen in dome/stupa monuments reflects a
genuine architectural-type signal rather than a statistical artefact of
the UNESCO corpus or its geographic distribution, then morphologically
distinct monument types should NOT show comparable enrichment.

Fortification-type sites (castles, citadels, fortresses, ramparts,
bastions) represent a class of monumental architecture that:
  • Is geographically widespread (Europe, Middle East, Asia, Americas)
  • Is well-documented in UNESCO metadata
  • Has no architectural or religious relationship to the beru/stupa
    tradition
  • Is large enough (N > 50 expected) to provide meaningful statistical
    power

If this corpus shows comparable Rayleigh R or A+/A++ enrichment, it
weakens the type-specificity claim. If it shows no enrichment, it
strengthens it.

KEYWORDS SEARCHED
-----------------
  Unambiguous (no context validation needed):
    citadel, citadels, fortification, fortifications, fortified,
    fortress, fortresses, bastion, bastions, rampart, ramparts,
    castle, castles

TESTS
-----
1. Corpus extraction: keyword match against full UNESCO descriptions
2. Tier distribution: A, A+, A++ counts vs geometric null (binomial)
3. Rayleigh phase-concentration permutation (phase randomisation,
   same method as circConcStupaP / primary result)
4. Circular-shift permutation (same method as Null A / LCO tests)
5. Comparison table vs dome corpus

OUTPUT
------
  Prints self-contained report + LaTeX macros to stdout.
  Macros written to ResultsStore.

  Key macros emitted (GROUP 30 — fortification negative control):
    \\fortN              — fortification corpus size
    \\fortTierAppN       — A++ count
    \\fortTierApN        — A+ count
    \\fortTierAN         — A count
    \\fortAppRate        — A++ rate (%)
    \\fortApRate         — A+ rate (%)
    \\fortARate          — A rate (%)
    \\fortBinomAppP      — binomial p vs A++ null
    \\fortBinomApP       — binomial p vs A+ null
    \\fortBinomAP        — binomial p vs A null
    \\fortRayleighR      — Rayleigh R (phase-randomisation)
    \\fortRayleighP      — Rayleigh p (phase-randomisation)
    \\fortCircShiftP     — circular-shift permutation p
    \\fortKeywords       — number of keywords used
"""

from __future__ import annotations

import sys
import re
import json
import numpy as np
from pathlib import Path
from scipy.stats import binomtest

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT))

from data.unesco_corpus import load_corpus, cultural_sites_with_coords_extended
from lib.beru import (
    GERIZIM, BERU,
    deviation as beru_dev,
    tier_label,
    TIER_APP, TIER_APLUS, TIER_A_MAX,
    P_NULL_APP, P_NULL_AP, P_NULL_A,
)
from lib.results_store import ResultsStore

_CFG        = json.loads((_ROOT / "config.json").read_text())
N_PERMS     = _CFG["simulation"]["n_permutations"]   # 100,000
SEED        = _CFG["simulation"]["random_seed"]
PERIOD_DEG  = 3.0   # 1 beru

# ── Fortification keywords ────────────────────────────────────────────────────
FORT_KEYWORDS = [
    "fortification", "fortifications", "fortified",
    "citadel", "citadels",
    "fortress", "fortresses",
    "bastion", "bastions",
    "rampart", "ramparts",
    "castle", "castles",
]

FORT_KW_RES = [
    re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    for kw in FORT_KEYWORDS
]

def is_fort_site(site) -> bool:
    """Return True if any fortification keyword appears in the site's text."""
    text = " ".join([
        getattr(site, "name", "") or "",
        getattr(site, "short_description", "") or "",
        getattr(site, "long_description", "") or "",
    ])
    return any(rx.search(text) for rx in FORT_KW_RES)


def rayleigh_phase_rand(lons: np.ndarray, n_perms: int, seed: int):
    """Phase-randomisation Rayleigh test (same as circConcStupaP)."""
    phases = (lons / BERU) * 2 * np.pi
    R_obs  = float(np.abs(np.mean(np.exp(1j * phases))))
    rng    = np.random.default_rng(seed)
    theta  = rng.uniform(0.0, 2 * np.pi, size=(n_perms, len(lons)))
    R_null = np.abs(np.mean(np.exp(1j * (phases[None, :] + theta)), axis=1))
    p_val  = float((R_null >= R_obs).mean())
    return round(R_obs, 4), round(p_val, 4)


def rayleigh_circ_shift(lons: np.ndarray, n_perms: int, seed: int):
    """Circular-shift Rayleigh test (same as Null A / LCO)."""
    phases    = (lons / BERU) * 2 * np.pi
    R_obs     = float(np.abs(np.mean(np.exp(1j * phases))))
    rng       = np.random.default_rng(seed + 1)
    shifts    = rng.uniform(0.0, PERIOD_DEG, n_perms)
    shift_rad = (shifts / BERU) * 2 * np.pi
    R_null    = np.abs(np.mean(
        np.exp(1j * (phases[None, :] + shift_rad[:, None])), axis=1
    ))
    p_val = float((R_null >= R_obs).mean())
    return round(R_obs, 4), round(p_val, 4)


def sig_label(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


def main():
    print("Loading UNESCO corpus (extended descriptions)…")
    all_sites = load_corpus()
    sites = cultural_sites_with_coords_extended(all_sites)
    print(f"  {len(sites)} cultural/mixed sites with coordinates")

    # ── Extract fortification sub-corpus ─────────────────────────────────────
    fort_sites = [s for s in sites if is_fort_site(s)]
    N = len(fort_sites)
    lons = np.array([s.longitude for s in fort_sites])

    print(f"  Fortification keyword hits: {N}")
    print()

    # ── Tier counts ───────────────────────────────────────────────────────────
    devs = np.array([beru_dev(lon) for lon in lons])
    n_app = int(np.sum(devs <= TIER_APP))
    n_ap  = int(np.sum(devs <= TIER_APLUS))
    n_a   = int(np.sum(devs <= TIER_A_MAX))
    n_b   = int(np.sum((devs > TIER_A_MAX) & (devs <= 0.5 - TIER_A_MAX)))

    rate_app = 100 * n_app / N if N else 0.0
    rate_ap  = 100 * n_ap  / N if N else 0.0
    rate_a   = 100 * n_a   / N if N else 0.0

    # Expected counts under geometric null
    exp_app = round(P_NULL_APP * N, 2)
    exp_ap  = round(P_NULL_AP  * N, 2)
    exp_a   = round(P_NULL_A   * N, 2)

    enrich_app = (n_app / N) / P_NULL_APP if N else 0.0
    enrich_ap  = (n_ap  / N) / P_NULL_AP  if N else 0.0
    enrich_a   = (n_a   / N) / P_NULL_A   if N else 0.0

    binom_app = binomtest(n_app, N, P_NULL_APP, alternative="greater")
    binom_ap  = binomtest(n_ap,  N, P_NULL_AP,  alternative="greater")
    binom_a   = binomtest(n_a,   N, P_NULL_A,   alternative="greater")

    # ── Rayleigh tests ────────────────────────────────────────────────────────
    R_pr, p_pr  = rayleigh_phase_rand(lons, N_PERMS, SEED)
    R_cs, p_cs  = rayleigh_circ_shift(lons, N_PERMS, SEED)

    # ── Print report ──────────────────────────────────────────────────────────
    sep  = "=" * 78
    dash = "─" * 78
    print(sep)
    print("  FORTIFICATION NEGATIVE CONTROL")
    print(f"  Keywords: {', '.join(FORT_KEYWORDS)}")
    print(f"  N = {N}  |  Anchor: {GERIZIM}°E  |  BERU: {BERU}°")
    print(sep)
    print()

    print(dash)
    print("  TIER DISTRIBUTION")
    print(dash)
    print(f"  {'Tier':<6}  {'Obs':>5}  {'Rate':>7}  {'Exp':>6}  {'Enrich':>7}  {'Binom-p':>9}  Sig")
    print(f"  {'────':<6}  {'───':>5}  {'────':>7}  {'───':>6}  {'──────':>7}  {'───────':>9}  ───")
    for label, obs, rate, exp, enr, bp in [
        ("A++", n_app, rate_app, exp_app, enrich_app, binom_app.pvalue),
        ("A+",  n_ap,  rate_ap,  exp_ap,  enrich_ap,  binom_ap.pvalue),
        ("A",   n_a,   rate_a,   exp_a,   enrich_a,   binom_a.pvalue),
    ]:
        print(f"  {label:<6}  {obs:>5}  {rate:>6.1f}%  {exp:>6.1f}  {enr:>6.2f}x  {bp:>9.4f}  {sig_label(bp)}")
    print()

    print(dash)
    print("  RAYLEIGH PHASE-CONCENTRATION TESTS")
    print(dash)
    print(f"  Phase-randomisation:  R = {R_pr:.4f},  p = {p_pr:.4f}  {sig_label(p_pr)}")
    print(f"  Circular-shift:       R = {R_cs:.4f},  p = {p_cs:.4f}  {sig_label(p_cs)}")
    print()

    print(dash)
    print("  COMPARISON: FORTIFICATION vs DOME (A++ primary)")
    print(dash)
    print(f"  {'Corpus':<20}  {'N':>5}  {'A++ obs':>8}  {'A++ rate':>9}  {'Binom-p':>9}  {'Rayleigh-p':>11}")
    print(f"  {'──────':<20}  {'─':>5}  {'───────':>8}  {'────────':>9}  {'───────':>9}  {'──────────':>11}")
    print(f"  {'Fortification':<20}  {N:>5}  {n_app:>8}  {rate_app:>8.1f}%  {binom_app.pvalue:>9.4f}  {p_pr:>11.4f}")
    print(f"  {'Dome/stupa':<20}  (see dome_geographic_concentration_test.py)")
    print()

    print(dash)
    print("  VERDICT")
    print(dash)
    if binom_app.pvalue >= 0.05 and p_pr >= 0.05:
        verdict = ("No significant enrichment or Rayleigh concentration detected "
                   "in fortification corpus. Supports type-specificity of dome signal.")
    elif binom_app.pvalue < 0.05 or p_pr < 0.05:
        verdict = ("WARNING: Significant result in fortification corpus. "
                   "This weakens the type-specificity argument. Investigate.")
    else:
        verdict = "Mixed results — interpret with caution."
    print(f"  {verdict}")
    print()

    # ── LaTeX macros (GROUP 30) ───────────────────────────────────────────────
    print(dash)
    print("  LATEX MACROS (GROUP 30 — fortification negative control)")
    print(dash)
    macros = [
        ("fortN",            str(N)),
        ("fortKeywords",     str(len(FORT_KEYWORDS))),
        ("fortTierAppN",     str(n_app)),
        ("fortTierApN",      str(n_ap)),
        ("fortTierAN",       str(n_a)),
        ("fortAppRate",      f"{rate_app:.1f}"),
        ("fortApRate",       f"{rate_ap:.1f}"),
        ("fortARate",        f"{rate_a:.1f}"),
        ("fortBinomAppP",    f"{binom_app.pvalue:.4f}"),
        ("fortBinomApP",     f"{binom_ap.pvalue:.4f}"),
        ("fortBinomAP",      f"{binom_a.pvalue:.4f}"),
        ("fortRayleighR",    f"{R_pr:.4f}"),
        ("fortRayleighP",    f"{p_pr:.4f}"),
        ("fortCircShiftP",   f"{p_cs:.4f}"),
    ]
    for name, val in macros:
        print(f"  \\newcommand{{\\{name}}}{{{val}}}")
    print()

    # ── Write to results store ────────────────────────────────────────────────
    store_data = {
        "fortN":          float(N),
        "fortKeywords":   float(len(FORT_KEYWORDS)),
        "fortTierAppN":   float(n_app),
        "fortTierApN":    float(n_ap),
        "fortTierAN":     float(n_a),
        "fortAppRate":    round(rate_app, 1),
        "fortApRate":     round(rate_ap,  1),
        "fortARate":      round(rate_a,   1),
        "fortBinomAppP":  round(binom_app.pvalue, 4),
        "fortBinomApP":   round(binom_ap.pvalue,  4),
        "fortBinomAP":    round(binom_a.pvalue,   4),
        "fortRayleighR":  R_pr,
        "fortRayleighP":  p_pr,
        "fortCircShiftP": p_cs,
    }
    ResultsStore().write_many(store_data)
    print("Results written to data/store/results.json")


if __name__ == "__main__":
    main()
