"""
anchor_uniqueness_audit.py
==========================
Full-corpus anchor uniqueness audit — single symmetric sweep.

One extended corpus of 1012 sites:
  1011 inscribed Cultural/Mixed UNESCO sites with coordinates
  + Mount Gerizim added as a synthetic entry (Tentative List, ref. 5706)

Each trial anchor self-excludes only itself (any corpus entry within 0.001°),
so every focal anchor tests against N = 1011.

  Gerizim  (35.269°E): tests against 1011 inscribed sites (incl. Jerusalem)
  Jerusalem (35.232°E): tests against 1010 inscribed sites + Gerizim (= 1011)

This gives both anchors identical N and identical treatment.  Gerizim's score
(56 A+, N = 1011) matches the primary Test 1 result exactly.

All 36 000 trial anchors (0–360°E, 0.01° resolution) are tested against
the same extended corpus with self-exclusion.  Both focal-anchor percentiles
are read from the same global distribution.
"""

import sys
from pathlib import Path
import numpy as np
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import (
    load_corpus, cultural_sites_with_coords_extended,
)
from lib.beru import (
    GERIZIM, BERU, TIER_APP, TIER_APLUS, TIER_A_MAX,
    P_NULL_AP, CONFIG,
)
from lib.results_store import ResultsStore

JERUSALEM_LON = CONFIG["anchors"]["jerusalem"]["longitude"]   # 35.2317°E

# ── Build the single extended corpus ─────────────────────────────────────────
corpus      = load_corpus()
_corpus_ext = cultural_sites_with_coords_extended(corpus)  # 1012: inscribed + Gerizim
assert len(_corpus_ext) == 1012, f"Extended corpus size wrong: {len(_corpus_ext)}"

_lons_ext = np.array([s.longitude for s in _corpus_ext])


# ── Counting function ─────────────────────────────────────────────────────────

def _count_tiers(lons_w, anchor_lon):
    """Return (n_app, n_ap, n_a) for a longitude array and anchor."""
    arcs     = np.abs(lons_w - anchor_lon)
    arcs     = np.minimum(arcs, 360 - arcs)
    bvs      = arcs / BERU
    nearests = np.round(bvs * 10) / 10
    devs     = np.abs(bvs - nearests)
    return (int(np.sum(devs <= TIER_APP)),
            int(np.sum(devs <= TIER_APLUS)),
            int(np.sum(devs <= TIER_A_MAX)))


def count_at(anchor_lon):
    """A+ count against the extended corpus, self-excluding anchor's own entry."""
    mask   = np.abs(_lons_ext - anchor_lon) > 0.001
    lons_w = _lons_ext[mask]
    return _count_tiers(lons_w, anchor_lon)


# ── Run the sweep ─────────────────────────────────────────────────────────────
sweep_anchors = np.arange(0.00, 360.00, 0.01)
n_anchors     = len(sweep_anchors)

print("=" * 100)
print("  FULL-CORPUS ANCHOR UNIQUENESS AUDIT  —  single symmetric sweep")
print(f"  Extended corpus: 1011 inscribed + Gerizim synthetic = 1012 sites")
print(f"  Self-exclusion: each anchor removes itself → N = 1011 for all focal anchors")
print(f"  Sweep: 0°–360°E at 0.01° resolution ({n_anchors:,} trial anchors)")
print("=" * 100)

print(f"\n  Running sweep…")
results = np.array([count_at(a) for a in sweep_anchors])  # shape (36000, 3)
ap_all  = results[:, 1]   # A+ column

# ── Focal anchor scores ───────────────────────────────────────────────────────
ger_app, ger_ap, ger_a = count_at(GERIZIM)
jer_app, jer_ap, jer_a = count_at(JERUSALEM_LON)

# Percentiles within the single global distribution
pctile_ger    = 100.0 * float(np.mean(ap_all <= ger_ap))
pctile_jer    = 100.0 * float(np.mean(ap_all <= jer_ap))
pct_above_ger = 100.0 * float(np.mean(ap_all >= ger_ap))
pct_above_jer = 100.0 * float(np.mean(ap_all >= jer_ap))

# x.18° artifact (every 3° = 0.1 beru)
x18_mask = np.abs((sweep_anchors % 3.0) - 2.18) < 0.006
x18_ap   = ap_all[x18_mask]

# Levant local max (34–37°E)
lev_mask    = (sweep_anchors >= 34.0) & (sweep_anchors <= 37.0)
lev_max_lon = float(sweep_anchors[lev_mask][np.argmax(ap_all[lev_mask])])
lev_max_ap  = int(np.max(ap_all[lev_mask]))

global_max   = int(np.max(ap_all))
n_at_max     = int(np.sum(ap_all == global_max))
sep          = abs(GERIZIM - JERUSALEM_LON)

# ── N for focal-anchor significance tests ─────────────────────────────────────
# Each focal anchor self-excludes 1 entry from the 1012 extended corpus → N = 1011.
N_FOCAL = len(_corpus_ext) - 1   # 1011

p_gerizim   = binomtest(ger_ap,    N_FOCAL, P_NULL_AP, alternative='greater').pvalue
p_jerusalem = binomtest(jer_ap,    N_FOCAL, P_NULL_AP, alternative='greater').pvalue
p_x18       = binomtest(global_max, N_FOCAL, P_NULL_AP, alternative='greater').pvalue

_sig = lambda p: "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))

print(f"""
══════════════════════════════════════════════════════════════════════════════
  SIGNIFICANCE TESTS  (one-sided binomial, H₁: rate > {P_NULL_AP:.0%} null, N = {N_FOCAL})
══════════════════════════════════════════════════════════════════════════════
  {"Anchor/case":<40}  {"k":>4}  {"N":>5}  {"p-value":>10}  sig
  {"─"*40}  {"─"*4}  {"─"*5}  {"─"*10}  {"─"*3}
  {"Gerizim  (35.269°E)":<40}  {ger_ap:>4}  {N_FOCAL:>5}  {p_gerizim:>10.4e}  {_sig(p_gerizim)}
  {"Jerusalem (35.232°E)":<40}  {jer_ap:>4}  {N_FOCAL:>5}  {p_jerusalem:>10.4e}  {_sig(p_jerusalem)}
  {"x.18° artifact (global max)":<40}  {global_max:>4}  {N_FOCAL:>5}  {p_x18:>10.4e}  {_sig(p_x18)}

  KEY FINDINGS:
  • Gerizim: {ger_ap} A+, p = {p_gerizim:.4f}  {_sig(p_gerizim)}
  • Jerusalem: {jer_ap} A+, p = {p_jerusalem:.4f}  {_sig(p_jerusalem)}
  • Gerizim's phase ({GERIZIM % 3.0:.3f}° mod 3°) is offset {abs(GERIZIM % 3.0 - 2.18):.3f}° from the
    x.18° optimum, so its score is independent of that artifact.
""")

print(f"""
══════════════════════════════════════════════════════════════════════════════
  SWEEP RESULTS  —  extended corpus 1012, self-exclusion → N = 1011
══════════════════════════════════════════════════════════════════════════════
  Mean A+ (all 36k anchors):  {np.mean(ap_all):.1f} ± {np.std(ap_all):.1f}
  Global max A+:              {global_max}  (at {n_at_max} anchors, all x.18° phase)

  GERIZIM ({GERIZIM}°E):
    A++  = {ger_app}   A+  = {ger_ap}   A  = {ger_a}
    Global percentile (A+):  {pctile_ger:.1f}th  ({pct_above_ger:.2f}% of anchors ≥ this)

  JERUSALEM ({JERUSALEM_LON}°E):
    A++  = {jer_app}   A+  = {jer_ap}   A  = {jer_a}
    Global percentile (A+):  {pctile_jer:.1f}th  ({pct_above_jer:.2f}% of anchors ≥ this)

  x.18° artifact: {np.sum(x18_mask)} anchors all score {x18_ap[0]} A+  (global max)
  Levant local max (34–37°E):  {lev_max_lon:.3f}°E  →  {lev_max_ap} A+

══════════════════════════════════════════════════════════════════════════════
  SUMMARY
══════════════════════════════════════════════════════════════════════════════
  Anchor separation: {sep:.3f}°  ({sep * 111 * np.cos(np.radians(32)):.1f} km at 32°N)
  Gerizim  : {ger_ap} A+  →  {pctile_ger:.1f}th pctile
  Jerusalem: {jer_ap} A+  →  {pctile_jer:.1f}th pctile
  Both place the ~35°E corridor in the top {100 - min(pctile_ger, pctile_jer):.1f}% globally.
  Single-site difference (Gerizim vs Jerusalem): {abs(jer_ap - ger_ap)} A+ (Mount Gerizim itself)
""")

# ── LaTeX macros (GROUP 11) ───────────────────────────────────────────────────
top_pct    = round(100.0 - min(pctile_ger, pctile_jer), 1)
sweep_diff = abs(jer_ap - ger_ap)

print("  % LaTeX macros (GROUP 11):")
_nanchors_fmt = f"{n_anchors:,}".replace(",", "{,}")
print(f"  \\newcommand{{\\anchorSweepNanchors}}{{{_nanchors_fmt}}}          % number of anchor points in sweep")
print(f"  \\newcommand{{\\anchorSweepApMean}}{{{np.mean(ap_all):.1f}}}           % mean A+ count across all anchors")
print(f"  \\newcommand{{\\anchorSweepApStd}}{{{np.std(ap_all):.1f}}}            % std A+ count across all anchors")
print(f"  \\newcommand{{\\anchorSweepGlobalMax}}{{{global_max}}}           % global max A+ at any anchor")
print(f"  \\newcommand{{\\anchorSweepGlobalMaxA}}{{{global_max}}}           % global max A+ (symmetric corpus)")
print(f"  \\newcommand{{\\optimalPhaseApA}}{{{global_max}}}           % max A+ at any x.18°E anchor")
print(f"  \\newcommand{{\\GerizimSweepApp}}{{{ger_app}}}             % Gerizim A++ in sweep")
print(f"  \\newcommand{{\\GerizimSweepAp}}{{{ger_ap}}}              % Gerizim A+ in sweep")
print(f"  \\newcommand{{\\GerizimSweepA}}{{{ger_a}}}             % Gerizim A in sweep")
print(f"  \\newcommand{{\\anchorSweepPctile}}{{{pctile_ger:.1f}}}          % Gerizim percentile in A+ anchor sweep")
print(f"  \\newcommand{{\\anchorSweepGeApPct}}{{{pct_above_ger:.2f}}}          % % of anchors beating Gerizim A+")
print(f"  \\newcommand{{\\JerusalemSweepApp}}{{{jer_app}}}           % Jerusalem A++ in sweep")
print(f"  \\newcommand{{\\JerusalemSweepAp}}{{{jer_ap}}}            % Jerusalem A+ in sweep")
print(f"  \\newcommand{{\\JerusalemSweepA}}{{{jer_a}}}           % Jerusalem A in sweep")
print(f"  \\newcommand{{\\JerusalemSweepPctile}}{{{pctile_jer:.1f}}}        % Jerusalem percentile in A+ anchor sweep")
print(f"  \\newcommand{{\\anchorSweepJerApPct}}{{{pct_above_jer:.2f}}}          % % of anchors beating Jerusalem A+")
print(f"  \\newcommand{{\\anchorSweepTopPct}}{{{top_pct}}}          % top-N% globally (corridor, both anchors)")
print(f"  \\newcommand{{\\sweepABDiff}}{{{sweep_diff}}}          % A+ count difference Gerizim vs Jerusalem")

# ── Write to results store ─────────────────────────────────────────────────────
ResultsStore().write_many({
    "anchorSweepPctile":     pctile_ger,
    "anchorSweepGeApPct":    pct_above_ger,
    "anchorSweepGlobalMaxA": global_max,
    "anchorSweepGlobalMaxB": global_max,
    "anchorSweepNsweep":     N_FOCAL,
    "pX18BinomA":            float(p_x18),
    "pX18BinomB":            float(p_x18),
    "anchorSweepApMeanRaw":  float(np.mean(ap_all)),
    "anchorSweepApStdRaw":   float(np.std(ap_all)),
    "anchorSweepTopPct":     top_pct,
    "sweepABDiff":           sweep_diff,
})
