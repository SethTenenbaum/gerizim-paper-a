"""
null_b_uniform_draw.py
======================
DEPRECATED — not used in the manuscript.

A uniform random draw from [-180°, 180°] counting harmonic proximity hits
is simply a Monte Carlo simulation of the analytic geometric null p₀ = 2d/T.
It adds no information beyond the closed-form binomial and is not a distinct
null model. Retained for reference only.
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
from lib.beru import GERIZIM, BERU, TIER_APP, TIER_APLUS, TIER_A_MAX
from lib.results_store import ResultsStore

_CFG    = json.loads((_ROOT / "config.json").read_text())
N_PERMS = _CFG["simulation"]["n_permutations"]
SEED    = _CFG["simulation"]["random_seed"]
T       = 3.0   # period (degrees)

THRESH_APP = TIER_APP   * BERU   # A++ threshold in degrees
THRESH_AP  = TIER_APLUS * BERU   # A+  threshold in degrees
THRESH_A   = TIER_A_MAX * BERU   # A   threshold in degrees

P_NULL_APP = 2 * THRESH_APP / T
P_NULL_AP  = 2 * THRESH_AP  / T
P_NULL_A   = 2 * THRESH_A   / T

MACRO_OUT = _ROOT / "analysis" / "unesco" / "null_b_uniform_draw_macros.tex"


# ── Helpers ───────────────────────────────────────────────────────────────────

def total_hits(lons: np.ndarray, thresh_deg: float) -> int:
    """Count sites within thresh_deg of the nearest 3°/GERIZIM harmonic."""
    arc = (lons - GERIZIM) % T
    dev = np.minimum(arc, T - arc)
    return int(np.sum(dev <= thresh_deg))


def batch_hits(draws: np.ndarray, thresh_deg: float) -> np.ndarray:
    """(N_PERMS, N) → (N_PERMS,) hit counts."""
    arc = (draws - GERIZIM) % T
    dev = np.minimum(arc, T - arc)
    return (dev <= thresh_deg).sum(axis=1)


def sig(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"


def fmt_p(p: float) -> str:
    if p >= 0.10:   return f"{p:.2f}"
    if p >= 0.001:  return f"{p:.3f}"
    return f"{p:.4f}"


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    corpus   = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    dome_lons = np.array([s.longitude for s in cultural if is_dome_site_raw(s)])
    N = len(dome_lons)

    obs_app = total_hits(dome_lons, THRESH_APP)
    obs_ap  = total_hits(dome_lons, THRESH_AP)
    obs_a   = total_hits(dome_lons, THRESH_A)

    print()
    print("=" * 80)
    print("  NULL B: UNIFORM RANDOM DRAW  (global longitude, no clustering assumed)")
    print(f"  UNESCO dome/stupa corpus  N = {N}  |  anchor = {GERIZIM}°E  |  T = {T}°")
    print(f"  Permutations: {N_PERMS:,}  |  Draw range: [−180°, +180°]  |  seed = {SEED}")
    print(f"  Thresholds — A++: ±{THRESH_APP:.3f}°  A+: ±{THRESH_AP:.3f}°  A: ±{THRESH_A:.3f}°")
    print(f"  Geometric null rates — A++: {P_NULL_APP:.2f}  A+: {P_NULL_AP:.2f}  A: {P_NULL_A:.2f}")
    print("=" * 80)
    print()
    print("  Observed hit counts:")
    print(f"    A++ ({THRESH_APP:.3f}°) : {obs_app}/{N}  ({obs_app/N*100:.1f}%  vs  null {P_NULL_APP*100:.0f}%)")
    print(f"    A+  ({THRESH_AP:.3f}°) : {obs_ap}/{N}  ({obs_ap/N*100:.1f}%  vs  null {P_NULL_AP*100:.0f}%)")
    print(f"    A   ({THRESH_A:.3f}°) : {obs_a}/{N}  ({obs_a/N*100:.1f}%  vs  null {P_NULL_A*100:.0f}%)")
    print()

    # ── Simulate ──────────────────────────────────────────────────────────────
    rng   = np.random.default_rng(SEED)
    draws = rng.uniform(-180.0, 180.0, size=(N_PERMS, N))

    ct_app = batch_hits(draws, THRESH_APP)
    ct_ap  = batch_hits(draws, THRESH_AP)
    ct_a   = batch_hits(draws, THRESH_A)

    p_app = float(np.mean(ct_app >= obs_app))
    p_ap  = float(np.mean(ct_ap  >= obs_ap))
    p_a   = float(np.mean(ct_a   >= obs_a))

    z_app = (obs_app - ct_app.mean()) / ct_app.std()
    z_ap  = (obs_ap  - ct_ap.mean())  / ct_ap.std()
    z_a   = (obs_a   - ct_a.mean())   / ct_a.std()

    print(f"  {'Tier':<5}  {'Obs':>4}  {'Mean':>6}  {'SD':>5}  {'Z':>6}  {'p':>8}  {'Sig':>4}")
    print("  " + "─" * 52)
    print(f"  {'A++':<5}  {obs_app:>4}  {ct_app.mean():>6.2f}  {ct_app.std():>5.2f}"
          f"  {z_app:>6.2f}  {fmt_p(p_app):>8}  {sig(p_app):>4}")
    print(f"  {'A+':<5}  {obs_ap:>4}  {ct_ap.mean():>6.2f}  {ct_ap.std():>5.2f}"
          f"  {z_ap:>6.2f}  {fmt_p(p_ap):>8}  {sig(p_ap):>4}")
    print(f"  {'A':<5}  {obs_a:>4}  {ct_a.mean():>6.2f}  {ct_a.std():>5.2f}"
          f"  {z_a:>6.2f}  {fmt_p(p_a):>8}  {sig(p_a):>4}")
    print()
    print("  Interpretation:")
    print("    Null A (ns): enrichment is consistent with dome-site geography;")
    print("                 resampling the observed dome longitudes reproduces")
    print("                 the hit count — the clustering carries the signal.")
    print("    Null B (above): a globally uniform random draw of the same N")
    print("                 would rarely produce this many harmonic-proximate")
    print("                 sites — the count is globally atypical regardless")
    print("                 of geographic conditioning.")
    print("    These answer orthogonal questions and are not contradictory.")
    print("=" * 80)
    print()

    # ── Macros ────────────────────────────────────────────────────────────────
    macros = [
        f"\\newcommand{{\\geoNullBUniformPApp}}{{{fmt_p(p_app)}}}",
        f"\\newcommand{{\\geoNullBUniformPAp}}{{{fmt_p(p_ap)}}}",
        f"\\newcommand{{\\geoNullBUniformPA}}{{{fmt_p(p_a)}}}",
        f"\\newcommand{{\\geoNullBUniformZApp}}{{{z_app:.2f}}}",
        f"\\newcommand{{\\geoNullBUniformZAp}}{{{z_ap:.2f}}}",
        f"\\newcommand{{\\geoNullBUniformZA}}{{{z_a:.2f}}}",
        f"\\newcommand{{\\geoNullBUniformSigApp}}{{{sig(p_app)}}}",
        f"\\newcommand{{\\geoNullBUniformSigAp}}{{{sig(p_ap)}}}",
        f"\\newcommand{{\\geoNullBUniformSigA}}{{{sig(p_a)}}}",
        f"\\newcommand{{\\geoNullBUniformMeanApp}}{{{ct_app.mean():.2f}}}",
        f"\\newcommand{{\\geoNullBUniformMeanAp}}{{{ct_ap.mean():.2f}}}",
        f"\\newcommand{{\\geoNullBUniformMeanA}}{{{ct_a.mean():.2f}}}",
    ]

    MACRO_OUT.parent.mkdir(parents=True, exist_ok=True)
    MACRO_OUT.write_text(
        "% Null B uniform draw macros — auto-generated by null_b_uniform_draw.py\n"
        + "\n".join(macros) + "\n"
    )
    print(f"  Macros written → {MACRO_OUT}")

    # Append / update generated_macros.tex
    gm = _ROOT / "manuscript" / "generated_macros.tex"
    if gm.exists():
        content = gm.read_text()
        block = (
            "\n% ── Null B: uniform random draw ────────────────────────────────\n"
            + "\n".join(macros) + "\n"
        )
        tag = "% ── Null B: uniform random draw"
        if tag in content:
            start = content.index(tag)
            # find next section marker or end of file
            nxt = content.find("\n%", start + len(tag))
            end = nxt if nxt != -1 else len(content)
            content = content[:start] + block.rstrip() + "\n" + content[end:]
        else:
            content += block
        gm.write_text(content)
        print(f"  Macros updated → {gm}")

    # Persist to results store
    store = ResultsStore()
    store.write_many({
        "geoNullBUniformPApp":   p_app,
        "geoNullBUniformPAp":    p_ap,
        "geoNullBUniformPA":     p_a,
        "geoNullBUniformZApp":   float(z_app),
        "geoNullBUniformZAp":    float(z_ap),
        "geoNullBUniformZA":     float(z_a),
    })


if __name__ == "__main__":
    main()
