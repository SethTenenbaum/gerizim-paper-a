"""
null_b_toroidal_shift.py
========================
DEPRECATED — not used in the manuscript.

A toroidal/circular shift of dome-site residuals tests whether cluster
phase alignment with the harmonic grid is non-random. However, since
the 3° harmonic grid tiles uniformly over all longitudes, regional
clustering at the ~100° Eurasian scale has no projection onto the 3°
period — a 20° cluster spans ~7 full periods, so its members each
independently face the same geometric null rate p₀ = 2d/T. The
toroidal shift therefore does not add a test orthogonal to Null A.
Retained for reference only.
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

MACRO_OUT = _ROOT / "analysis" / "unesco" / "null_b_toroidal_shift_macros.tex"


# ── Helpers ───────────────────────────────────────────────────────────────────

def total_hits(lons: np.ndarray, thresh_deg: float) -> int:
    """Count sites within thresh_deg of the nearest harmonic (vectorised)."""
    arc = (lons - GERIZIM) % T
    dev = np.minimum(arc, T - arc)
    return int(np.sum(dev <= thresh_deg))


def shift_counts(residuals: np.ndarray, offsets: np.ndarray,
                 thresh_deg: float) -> np.ndarray:
    """
    Parameters
    ----------
    residuals : (N_dome,)  pre-computed (lon − GERIZIM) mod T
    offsets   : (N_PERMS,) random phase offsets in [0, T)
    thresh_deg: proximity threshold in degrees

    Returns
    -------
    (N_PERMS,) integer hit counts under the shifted grid.
    """
    # Broadcast: (N_PERMS, N_dome)
    shifted = (residuals[None, :] + offsets[:, None]) % T
    dev = np.minimum(shifted, T - shifted)
    return (dev <= thresh_deg).sum(axis=1).astype(int)


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
    corpus    = load_corpus()
    cultural  = cultural_sites_with_coords(corpus)
    dome_lons = np.array([s.longitude for s in cultural if is_dome_site_raw(s)])
    N = len(dome_lons)

    obs_app = total_hits(dome_lons, THRESH_APP)
    obs_ap  = total_hits(dome_lons, THRESH_AP)
    obs_a   = total_hits(dome_lons, THRESH_A)

    # Pre-compute residuals once
    residuals = (dome_lons - GERIZIM) % T   # shape (N,)

    print()
    print("=" * 80)
    print("  NULL B (REVISED): TOROIDAL / CIRCULAR SHIFT TEST")
    print(f"  UNESCO dome/stupa corpus  N = {N}  |  anchor = {GERIZIM}°E  |  T = {T}°")
    print(f"  Permutations: {N_PERMS:,}  |  Offsets: Uniform[0, {T})  |  seed = {SEED}")
    print(f"  Thresholds — A++: ±{THRESH_APP:.3f}°  A+: ±{THRESH_AP:.3f}°  A: ±{THRESH_A:.3f}°")
    print()
    print("  Design: all inter-site distances preserved; only grid alignment is")
    print("  randomised.  Tests whether cluster phase alignment is non-random.")
    print("=" * 80)
    print()
    print("  Observed hit counts:")
    print(f"    A++ ({THRESH_APP:.3f}°) : {obs_app}/{N}  ({obs_app/N*100:.1f}%  vs  null {P_NULL_APP*100:.0f}%)")
    print(f"    A+  ({THRESH_AP:.3f}°) : {obs_ap}/{N}  ({obs_ap/N*100:.1f}%  vs  null {P_NULL_AP*100:.0f}%)")
    print(f"    A   ({THRESH_A:.3f}°) : {obs_a}/{N}  ({obs_a/N*100:.1f}%  vs  null {P_NULL_A*100:.0f}%)")
    print()

    # ── Simulate ──────────────────────────────────────────────────────────────
    rng     = np.random.default_rng(SEED)
    offsets = rng.uniform(0.0, T, size=N_PERMS)   # shape (N_PERMS,)

    ct_app = shift_counts(residuals, offsets, THRESH_APP)
    ct_ap  = shift_counts(residuals, offsets, THRESH_AP)
    ct_a   = shift_counts(residuals, offsets, THRESH_A)

    p_app = float(np.mean(ct_app >= obs_app))
    p_ap  = float(np.mean(ct_ap  >= obs_ap))
    p_a   = float(np.mean(ct_a   >= obs_a))

    z_app = float((obs_app - ct_app.mean()) / ct_app.std())
    z_ap  = float((obs_ap  - ct_ap.mean())  / ct_ap.std())
    z_a   = float((obs_a   - ct_a.mean())   / ct_a.std())

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
    print("    A non-significant result indicates that the observed clusters,")
    print("    if shifted to a random phase within the harmonic period, would")
    print("    produce comparable hit counts — i.e. the alignment is not")
    print("    surprising given the clustering structure.")
    print("    A significant result indicates the phase alignment of the clusters")
    print("    with the harmonic grid is itself non-random, independent of")
    print("    how tightly the monuments cluster.")
    print("=" * 80)
    print()

    # ── LaTeX macros ──────────────────────────────────────────────────────────
    macros = [
        f"\\newcommand{{\\geoNullBShiftPApp}}{{{fmt_p(p_app)}}}",
        f"\\newcommand{{\\geoNullBShiftPAp}}{{{fmt_p(p_ap)}}}",
        f"\\newcommand{{\\geoNullBShiftPA}}{{{fmt_p(p_a)}}}",
        f"\\newcommand{{\\geoNullBShiftZApp}}{{{z_app:.2f}}}",
        f"\\newcommand{{\\geoNullBShiftZAp}}{{{z_ap:.2f}}}",
        f"\\newcommand{{\\geoNullBShiftZA}}{{{z_a:.2f}}}",
        f"\\newcommand{{\\geoNullBShiftSigApp}}{{{sig(p_app)}}}",
        f"\\newcommand{{\\geoNullBShiftSigAp}}{{{sig(p_ap)}}}",
        f"\\newcommand{{\\geoNullBShiftSigA}}{{{sig(p_a)}}}",
        f"\\newcommand{{\\geoNullBShiftMeanApp}}{{{ct_app.mean():.2f}}}",
        f"\\newcommand{{\\geoNullBShiftMeanAp}}{{{ct_ap.mean():.2f}}}",
        f"\\newcommand{{\\geoNullBShiftMeanA}}{{{ct_a.mean():.2f}}}",
    ]

    MACRO_OUT.parent.mkdir(parents=True, exist_ok=True)
    MACRO_OUT.write_text(
        "% Null B toroidal shift macros — auto-generated by null_b_toroidal_shift.py\n"
        + "\n".join(macros) + "\n"
    )
    print(f"  Macros written → {MACRO_OUT}")

    # Append / update generated_macros.tex
    gm = _ROOT / "manuscript" / "generated_macros.tex"
    if gm.exists():
        content = gm.read_text()
        block = (
            "\n% ── Null B: toroidal shift ──────────────────────────────────────\n"
            + "\n".join(macros) + "\n"
        )
        tag = "% ── Null B: toroidal shift"
        if tag in content:
            start = content.index(tag)
            nxt   = content.find("\n%", start + len(tag))
            end   = nxt if nxt != -1 else len(content)
            content = content[:start] + block.rstrip() + "\n" + content[end:]
        else:
            # Remove the old uniform-draw block if present, then append new
            old_tag = "% ── Null B: uniform random draw"
            if old_tag in content:
                start = content.index(old_tag)
                nxt   = content.find("\n%", start + len(old_tag))
                end   = nxt if nxt != -1 else len(content)
                content = content[:start] + block.rstrip() + "\n" + content[end:]
            else:
                content += block
        gm.write_text(content)
        print(f"  Macros updated → {gm}")

    # Persist to results store
    store = ResultsStore()
    store.write_many({
        "geoNullBShiftPApp":  p_app,
        "geoNullBShiftPAp":   p_ap,
        "geoNullBShiftPA":    p_a,
        "geoNullBShiftZApp":  z_app,
        "geoNullBShiftZAp":   z_ap,
        "geoNullBShiftZA":    z_a,
    })
    print(f"  Results written → data/store/results.json")


if __name__ == "__main__":
    main()
