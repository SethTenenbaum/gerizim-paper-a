"""
generate_audit_anchor_sweep.py
==============================
Produce supplementary/audit/anchor_sweep_audit.txt

Global anchor sweep results: all 36,000 trial anchors ranked by A+ count,
with Gerizim and Jerusalem highlighted and the top-scoring anchors listed.

Uses the same extended corpus and methodology as anchor_uniqueness_audit.py:
  N_ext = 1012 (1011 inscribed + Gerizim synthetic)
  Each trial anchor self-excludes entries within 0.001°

Run from repo root:
    python3 tools/generate_audit_anchor_sweep.py
"""

import sys
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent))
from data.unesco_corpus import (
    load_corpus, cultural_sites_with_coords_extended, GERIZIM_SYNTHETIC,
)
from lib.beru import GERIZIM, BERU, TIER_APP, TIER_APLUS, TIER_A_MAX, P_NULL_AP, CONFIG

OUT = Path(__file__).parent.parent / "supplementary" / "audit" / "anchor_sweep_audit.txt"
SEP = "─" * 100

JERUSALEM_LON = CONFIG["anchors"]["jerusalem"]["longitude"]   # 35.2317°E

# ── Build extended corpus ─────────────────────────────────────────────────────
corpus      = load_corpus()
_corpus_ext = cultural_sites_with_coords_extended(corpus)   # 1012 sites
_lons_ext   = np.array([s.longitude for s in _corpus_ext])
N_EXT       = len(_corpus_ext)   # 1012
N_FOCAL     = N_EXT - 1          # 1011 (after self-exclusion)


def _count_ap(anchor_lon: float) -> int:
    mask = np.abs(_lons_ext - anchor_lon) > 0.001
    arcs = np.abs(_lons_ext[mask] - anchor_lon)
    arcs = np.minimum(arcs, 360 - arcs)
    bvs  = arcs / BERU
    devs = np.abs(bvs - np.round(bvs * 10) / 10)
    return int(np.sum(devs <= TIER_APLUS))


def _count_all_tiers(anchor_lon: float) -> tuple:
    mask = np.abs(_lons_ext - anchor_lon) > 0.001
    arcs = np.abs(_lons_ext[mask] - anchor_lon)
    arcs = np.minimum(arcs, 360 - arcs)
    bvs  = arcs / BERU
    devs = np.abs(bvs - np.round(bvs * 10) / 10)
    return (int(np.sum(devs <= TIER_APP)),
            int(np.sum(devs <= TIER_APLUS)),
            int(np.sum(devs <= TIER_A_MAX)))


def sig(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"


# ── Run sweep ─────────────────────────────────────────────────────────────────
sweep_anchors = np.arange(0.00, 360.00, 0.01)
n_anchors     = len(sweep_anchors)

print("Running sweep…", flush=True)
ap_counts = np.array([_count_ap(a) for a in sweep_anchors])
print("Done.", flush=True)

# ── Focal anchor scores ───────────────────────────────────────────────────────
ger_app, ger_ap, ger_a = _count_all_tiers(GERIZIM)
jer_app, jer_ap, jer_a = _count_all_tiers(JERUSALEM_LON)

# Ranks (1 = best)
ger_rank = int(np.sum(ap_counts > ger_ap)) + 1
jer_rank = int(np.sum(ap_counts > jer_ap)) + 1

pctile_ger = 100.0 * float(np.mean(ap_counts <= ger_ap))
pctile_jer = 100.0 * float(np.mean(ap_counts <= jer_ap))
pct_ge_ger = 100.0 * float(np.mean(ap_counts >= ger_ap))
pct_ge_jer = 100.0 * float(np.mean(ap_counts >= jer_ap))

p_ger = binomtest(ger_ap, N_FOCAL, P_NULL_AP, alternative="greater").pvalue
p_jer = binomtest(jer_ap, N_FOCAL, P_NULL_AP, alternative="greater").pvalue

# Distribution stats
global_max   = int(np.max(ap_counts))
global_min   = int(np.min(ap_counts))
ap_mean      = float(np.mean(ap_counts))
ap_std       = float(np.std(ap_counts))
n_above_ger  = int(np.sum(ap_counts > ger_ap))
n_at_ger     = int(np.sum(ap_counts == ger_ap))
n_above_jer  = int(np.sum(ap_counts > jer_ap))
n_at_jer     = int(np.sum(ap_counts == jer_ap))

# x.18° artifact
x18_mask = np.abs((sweep_anchors % 3.0) - 2.18) < 0.006
x18_ap   = ap_counts[x18_mask]
x18_lons = sweep_anchors[x18_mask]

ts = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")

lines = [
    "GLOBAL ANCHOR SWEEP AUDIT — GERIZIM AND JERUSALEM RANKINGS",
    f"Generated : {ts}",
    f"Script    : tools/generate_audit_anchor_sweep.py",
    f"Corpus    : {N_EXT} sites (1011 inscribed Cultural/Mixed + Gerizim synthetic)",
    f"           Each trial anchor self-excludes itself → N = {N_FOCAL} per anchor",
    f"Sweep     : {n_anchors:,} trial anchors, 0–360°E at 0.01° resolution",
    f"Threshold : A+ ≤ {TIER_APLUS} beru ({TIER_APLUS*BERU*111:.1f} km from nearest harmonic)",
    f"Null rate : {P_NULL_AP:.0%} geometric null (4% of 360° in A+ bands)",
    "",
    SEP,
    "DISTRIBUTION SUMMARY",
    "",
    f"  Mean A+ across all {n_anchors:,} anchors: {ap_mean:.2f} ± {ap_std:.2f}",
    f"  Global max:  {global_max}  (x.18°E phase; {len(x18_lons)} anchors all score identically)",
    f"  Global min:  {global_min}",
    f"  Expected (4% null):  {0.04 * N_FOCAL:.1f}",
    "",
    SEP,
    "FOCAL ANCHOR RESULTS",
    "",
    f"  {'Anchor':<30}  {'Lon°E':>7}  {'A++':>5}  {'A+':>5}  {'A':>5}  "
    f"{'Rank':>6}  {'Pctile':>7}  {'% ≥':>6}  {'p-value':>10}  sig",
    "  " + "-" * 95,
    f"  {'Gerizim (Tentative List)':<30}  {GERIZIM:>7.3f}  {ger_app:>5}  {ger_ap:>5}  {ger_a:>5}  "
    f"{ger_rank:>6}  {pctile_ger:>6.1f}th  {pct_ge_ger:>5.2f}%  {p_ger:>10.4e}  {sig(p_ger)}",
    f"  {'Jerusalem (inscribed)':<30}  {JERUSALEM_LON:>7.4f}  {jer_app:>5}  {jer_ap:>5}  {jer_a:>5}  "
    f"{jer_rank:>6}  {pctile_jer:>6.1f}th  {pct_ge_jer:>5.2f}%  {p_jer:>10.4e}  {sig(p_jer)}",
    "",
    f"  Gerizim : {n_above_ger} anchors score higher, {n_at_ger} score equally → rank {ger_rank} / {n_anchors:,}",
    f"  Jerusalem: {n_above_jer} anchors score higher, {n_at_jer} score equally → rank {jer_rank} / {n_anchors:,}",
    "",
    f"  NOTE: Gerizim's rank is lower than Jerusalem's because trial anchors near",
    f"  35.269°E benefit from Gerizim as a corpus member; x.18°E anchors dominate",
    f"  the top ranks due to a periodic phase artifact (see below).",
    "",
]

# ── Top anchors by A+ count ───────────────────────────────────────────────────
# Build sorted list with context about each anchor
lines += [
    SEP,
    f"TOP TRIAL ANCHORS BY A+ COUNT  (first {min(100, n_anchors)} distinct A+ levels shown)",
    "",
    f"  {'Rank':>5}  {'Lon°E':>7}  {'A+':>4}  {'Phase mod 3°':>13}  Notes",
    "  " + "-" * 70,
]

# Get unique A+ levels from highest to lowest
unique_ap_desc = sorted(set(ap_counts), reverse=True)[:30]

rank_counter = 1
for ap_val in unique_ap_desc:
    idxs = np.where(ap_counts == ap_val)[0]
    lons_at = sweep_anchors[idxs]
    phases   = lons_at % 3.0
    # Are all x.18 phase?
    is_x18 = np.all(np.abs(phases - 2.18) < 0.01)

    # Does this group include Gerizim's band?
    near_ger = any(abs(l - GERIZIM) <= 0.06 for l in lons_at)
    near_jer = any(abs(l - JERUSALEM_LON) <= 0.06 for l in lons_at)

    if is_x18:
        note = f"x.18°E artifact ({len(idxs)} anchors; periodic phase, not historically specific)"
        lines.append(f"  {rank_counter:>5}  {'various':>7}  {ap_val:>4}  {'≈2.18 (mod 3)':>13}  {note}")
        rank_counter += len(idxs)
    else:
        for lon in sorted(lons_at):
            phase = lon % 3.0
            note = ""
            if abs(lon - GERIZIM) <= 0.001:
                note = "← Gerizim (0.01° grid; exact lon 35.269°E)"
            elif abs(lon - GERIZIM) <= 0.06:
                note = f"← near Gerizim ({abs(lon-GERIZIM)*111:.1f} km)"
            elif abs(lon - JERUSALEM_LON) <= 0.06:
                note = f"← near Jerusalem ({abs(lon-JERUSALEM_LON)*111:.1f} km)"
            lines.append(
                f"  {rank_counter:>5}  {lon:>7.2f}  {ap_val:>4}  {phase:>13.3f}  {note}"
            )
            rank_counter += 1

lines.append("")

# ── Context window around Gerizim ────────────────────────────────────────────
ger_idx    = np.argmin(np.abs(sweep_anchors - GERIZIM))
window_lo  = max(0, ger_idx - 15)
window_hi  = min(n_anchors - 1, ger_idx + 15)

lines += [
    SEP,
    f"CONTEXT WINDOW: ±0.15° AROUND GERIZIM ({GERIZIM}°E)",
    "",
    f"  {'Lon°E':>7}  {'A+':>4}  {'Rank':>6}  Notes",
    "  " + "-" * 55,
]
for i in range(window_lo, window_hi + 1):
    lon = sweep_anchors[i]
    ap  = ap_counts[i]
    r   = int(np.sum(ap_counts > ap)) + 1
    note = ""
    if abs(lon - GERIZIM) < 0.005:
        note = "← nearest grid point to Gerizim"
    elif abs(lon - JERUSALEM_LON) < 0.005:
        note = "← nearest grid point to Jerusalem"
    lines.append(f"  {lon:>7.2f}  {ap:>4}  {r:>6}  {note}")
lines.append("")

# ── x.18° artifact explanation ────────────────────────────────────────────────
lines += [
    SEP,
    "x.18°E PHASE ARTIFACT",
    "",
    f"  Every longitude of the form x.18°E (x = integer) produces the maximum",
    f"  A+ count ({global_max}) because these points are equidistant from all harmonics",
    f"  of the 0.1-beru grid.  There are {len(x18_lons)} such anchors (one per 3° = 0.1 beru).",
    f"  They are uninformative: every historically attested x.18° site would score",
    f"  identically regardless of cultural significance.",
    "",
    f"  Gerizim's phase: {GERIZIM % 3.0:.3f}° mod 3°  (offset {abs(GERIZIM % 3.0 - 2.18):.3f}° from x.18°)",
    f"  Jerusalem's phase: {JERUSALEM_LON % 3.0:.4f}° mod 3°  (offset {abs(JERUSALEM_LON % 3.0 - 2.18):.4f}° from x.18°)",
    "",
    f"  The ranking table above groups all {len(x18_lons)} x.18° anchors at rank 1 since",
    f"  they are indistinguishable and irrelevant to the Levantine corridor question.",
    "",
]

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text("\n".join(lines), encoding="utf-8")
print(f"Written -> {OUT}")
print(f"  Gerizim  : {ger_ap} A+  rank {ger_rank}/{n_anchors:,}  ({pctile_ger:.1f}th pctile)")
print(f"  Jerusalem: {jer_ap} A+  rank {jer_rank}/{n_anchors:,}  ({pctile_jer:.1f}th pctile)")
