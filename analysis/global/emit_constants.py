"""
emit_constants.py
=================
Emits LaTeX macros for pure constants that are fully determined by
config.json and lib/beru.py — no corpus loading required.

GROUP 0 — Geographic and metrological constants:
  \\GerizimLon, \\GerizimLonDMS, \\GerizimLatDMS
  \\JerusalemLon
  \\GerizimJerusalemSep, \\GerizimJerusalemKm
  \\NwhcTotal

USAGE
-----
  python3 analysis/global/emit_constants.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import json
from data.unesco_corpus import load_corpus

_cfg = json.loads((ROOT / "config.json").read_text())

GER_LON     = _cfg["anchors"]["gerizim"]["longitude"]
GER_LAT_DMS = _cfg["anchors"]["gerizim"]["lat_dms"]       # e.g. "N32~12~44"
GER_LON_DMS = _cfg["anchors"]["gerizim"]["lon_dms"]       # e.g. "E35~16~08"
JER_LON     = _cfg["anchors"]["jerusalem"]["longitude"]
KM_PER_DEG  = _cfg["units"]["km_per_degree"]              # 111.0 (equatorial convention)

# Gerizim–Jerusalem separation (equatorial, no cos-lat, per study convention)
sep_deg = round(abs(GER_LON - JER_LON), 3)
sep_km  = round(sep_deg * KM_PER_DEG, 1)

# Total UNESCO corpus size (all site types, including natural)
corpus    = load_corpus()
n_whc     = len(corpus)

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

print()
print("  ✓ All GROUP 0 constants emitted.")
