"""
emit_constants.py
=================
Emits LaTeX macros for pure constants that are fully determined by
config.json and lib/beru.py — no corpus loading required.

GROUP 0 — Geographic and metrological constants:
  \\GerizimLon, \\GerizimLonDMS, \\GerizimLatDMS
  \\JerusalemLon
  \\GerizimJerusalemSep, \\GerizimJerusalemKm
  \\NwhcTotal, \\NclusterTotalPreCoord, \\NclusterTotal

USAGE
-----
  python3 analysis/global/emit_constants.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import json
import collections
from data.unesco_corpus import load_corpus
from lib.beru import (
    BERU as BERU_DEG, HARMONIC_STEP,
    TIER_APP, TIER_APLUS, TIER_A_MAX,
    TIER_APP_KM, TIER_APLUS_KM, TIER_A_KM,
    P_NULL_APP, P_NULL_AP, P_NULL_A,
    TIER_CMINUS, TIER_CMINUS2,
)

_cfg = json.loads((ROOT / "config.json").read_text())

# ── Metrological constants from config ────────────────────────────────────────
N_HARMONICS      = int(round(360.0 / (BERU_DEG * HARMONIC_STEP)))   # 120
NULL_RATE_APP    = P_NULL_APP    # = 2 * TIER_APP   / HARMONIC_STEP
NULL_RATE_AP     = P_NULL_AP     # = 2 * TIER_APLUS / HARMONIC_STEP
NULL_RATE_A      = P_NULL_A      # = 2 * TIER_A_MAX / HARMONIC_STEP
SIM_SEED         = _cfg["simulation"]["random_seed"]         # 42

GER_LON     = _cfg["anchors"]["gerizim"]["longitude"]
GER_LAT_DMS = _cfg["anchors"]["gerizim"]["lat_dms"]       # e.g. "N32~12~44"
GER_LON_DMS = _cfg["anchors"]["gerizim"]["lon_dms"]       # e.g. "E35~16~08"
JER_LON     = _cfg["anchors"]["jerusalem"]["longitude"]
KM_PER_DEG  = _cfg["units"]["km_per_degree"]              # 111.0 (equatorial convention)

# Gerizim–Jerusalem separation (equatorial, no cos-lat, per study convention)
sep_deg = round(abs(GER_LON - JER_LON), 3)
sep_km  = round(sep_deg * KM_PER_DEG, 1)

# Total UNESCO corpus size (all site types, including natural)
from data.unesco_corpus import cultural_sites_with_coords
corpus    = load_corpus()
n_whc     = len(corpus)
n_cult_mixed         = sum(1 for s in corpus if s.is_cultural_or_mixed)   # after category filter, before coord filter
_cm_coords = cultural_sites_with_coords(corpus)
n_cult_mixed_coords  = len(_cm_coords)             # final analysis corpus

# Regional breakdown of the analysis corpus (Cultural/Mixed with coordinates)
_region_counts = collections.Counter(s.regions for s in _cm_coords)
_n_corpus = n_cult_mixed_coords
# Primary UNESCO regions
_n_europe = _region_counts.get("Europe and North America", 0)
_n_asia   = _region_counts.get("Asia and the Pacific", 0)
_n_latam  = _region_counts.get("Latin America and the Caribbean", 0)
_n_arab   = _region_counts.get("Arab States", 0)
_n_africa = _region_counts.get("Africa", 0)

print("=" * 72)
print("  EMIT CONSTANTS — GROUP 0")
print("  Source: config.json + data/unesco_corpus.py")
print("=" * 72)

print()
print("  % LaTeX macros (GROUP 0):")
print(f"  \\newcommand{{\\GerizimLon}}{{{GER_LON}}}  % anchor longitude, decimal degrees E")
print(f"  \\newcommand{{\\GerizimLonDMS}}{{{GER_LON_DMS}}}  % anchor longitude, DMS (from config)")
print(f"  \\newcommand{{\\GerizimLatDMS}}{{{GER_LAT_DMS}}}  % anchor latitude, DMS (from config)")
print(f"  \\newcommand{{\\JerusalemLon}}{{{JER_LON:.3f}}}  % Jerusalem longitude, decimal degrees E")
print(f"  \\newcommand{{\\GerizimJerusalemSep}}{{{sep_deg}}}  % Gerizim–Jerusalem separation (degrees, equatorial)")
print(f"  \\newcommand{{\\GerizimJerusalemKm}}{{{sep_km}}}  % Gerizim–Jerusalem separation (km, equatorial)")
print(f"  \\newcommand{{\\NwhcTotal}}{{{n_whc}}}  % total UNESCO WHC inscribed sites (all types)")
print(f"  \\newcommand{{\\NclusterTotalPreCoord}}{{{n_cult_mixed}}}  % Cultural+Mixed before coordinate filter")
print(f"  \\newcommand{{\\NclusterTotal}}{{{n_cult_mixed_coords}}}  % Cultural+Mixed with valid coordinates (analysis corpus)")
print(f"  \\newcommand{{\\NexcNatural}}{{{n_whc - n_cult_mixed}}}  % Natural-only sites excluded from analysis corpus")
print(f"  \\newcommand{{\\NexcNoCoord}}{{{n_cult_mixed - n_cult_mixed_coords}}}  % Cultural/Mixed sites excluded for missing coordinates")
print(f"  \\newcommand{{\\NfineSweepSpacings}}{{7}}  % number of spacings in the fine unit sweep (Table tab:unitsweep_fine)")
print(f"  \\newcommand{{\\NharmonicsCoarse}}{{{N_HARMONICS}}}  % number of 0.1-beru harmonics in 360 degrees")
print(f"  \\newcommand{{\\simSeed}}{{{SIM_SEED}}}  % random seed used in all permutation/simulation tests")
print(f"  \\newcommand{{\\NullRateApp}}{{{NULL_RATE_APP*100:.1f}\\%}}  % geometric null rate for Tier-A++")
print(f"  \\newcommand{{\\NullRateAp}}{{{NULL_RATE_AP*100:.0f}\\%}}  % geometric null rate for Tier-A+")
print(f"  \\newcommand{{\\NullRateA}}{{{NULL_RATE_A*100:.0f}\\%}}  % geometric null rate for Tier-A")
# ── Tier threshold macros ────────────────────────────────────────────────────
print(f"  \\newcommand{{\\ApThreshBeru}}{{{TIER_APLUS}}}  % A+ threshold (beru)")
print(f"  \\newcommand{{\\ApThreshKm}}{{{TIER_APLUS_KM:.1f}}}  % A+ threshold (km, equatorial)")
print(f"  \\newcommand{{\\AppThreshBeru}}{{{TIER_APP}}}  % A++ threshold (beru)")
print(f"  \\newcommand{{\\AppThreshKm}}{{{TIER_APP_KM:.2f}}}  % A++ threshold (km, equatorial)")
print(f"  \\newcommand{{\\ATierThreshBeru}}{{{TIER_A_MAX}}}  % A-tier threshold (beru)")
print(f"  \\newcommand{{\\ATierThreshKm}}{{{TIER_A_KM:.0f}}}  % A-tier threshold (km, equatorial)")
print(f"  \\newcommand{{\\CMinusThreshBeru}}{{{TIER_CMINUS}}}  % C− threshold (beru, = A+ mirror)")
print(f"  \\newcommand{{\\CMinusThreshKm}}{{{TIER_CMINUS * BERU_DEG * 111.0:.1f}}}  % C− threshold (km)")
print(f"  \\newcommand{{\\NullRateApPct}}{{{NULL_RATE_AP*100:.4g}}}  % A+ null rate as plain number (no %)")
print(f"  \\newcommand{{\\NullRateApFormula}}{{$2 \\times \\ApThreshBeru{{}} / 0.1$}}  % null rate formula")
# ── Static convenience macros ─────────────────────────────────────────────────
print(f"  \\newcommand{{\\threeStar}}{{***}}  % *** significance shorthand for p<0.001 literal comparisons")
print(f"  \\newcommand{{\\twoStar}}{{**}}    % ** significance shorthand")
print(f"  \\newcommand{{\\oneStar}}{{*}}     % * significance shorthand")
print(f"  \\newcommand{{\\nsLabel}}{{ns}}    % ns label — use for editorial ns where no single p-value macro exists")

# ── Regional fractions (Cultural/Mixed analysis corpus) ───────────────────────
print(f"  \\newcommand{{\\CorpusPctEurope}}{{{round(100*_n_europe/_n_corpus)}}}  % pct Europe & North America in analysis corpus")
print(f"  \\newcommand{{\\CorpusPctAsiaPac}}{{{round(100*_n_asia/_n_corpus)}}}  % pct Asia and the Pacific in analysis corpus")
print(f"  \\newcommand{{\\CorpusPctLatAm}}{{{round(100*_n_latam/_n_corpus)}}}  % pct Latin America & Caribbean in analysis corpus")
print(f"  \\newcommand{{\\CorpusPctArabStates}}{{{round(100*_n_arab/_n_corpus)}}}  % pct Arab States in analysis corpus")
print(f"  \\newcommand{{\\CorpusPctAfrica}}{{{round(100*_n_africa/_n_corpus)}}}  % pct Africa in analysis corpus")

print()
print("  ✓ All GROUP 0 constants emitted.")

# ── ResultsStore ─────────────────────────────────────────────────────────────
from lib.results_store import ResultsStore
ResultsStore().write_many({
    "ApThreshBeru":    TIER_APLUS,
    "ApThreshKm":      round(TIER_APLUS_KM, 1),
    "AppThreshBeru":   TIER_APP,
    "AppThreshKm":     round(TIER_APP_KM, 2),
    "ATierThreshBeru": TIER_A_MAX,
    "ATierThreshKm":   round(TIER_A_KM, 0),
    "CMinusThreshBeru": TIER_CMINUS,
    "NullRateApPct":   round(NULL_RATE_AP * 100, 4),
})
