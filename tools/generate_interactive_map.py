#!/usr/bin/env python3
"""
generate_interactive_map.py — Build the interactive Leaflet HTML tier map.

Reads:
    keywords.json              — religion sets, mound evolution, founding keywords
    config.json (via lib/beru) — GERIZIM anchor, BERU unit, tier thresholds
    data/unesco_corpus         — UNESCO cultural/mixed sites with coordinates
    data/store/unesco/wikidata_stupas_q180987.csv — Wikidata Q180987 stupas

Writes:
    supplementary/interactive_map/silkroad_interactive.html  — self-contained browser map

Re-run whenever the corpus, keywords.json, or config.json changes:
    cd gerizim-paper-a
    python3 tools/generate_interactive_map.py
"""

import csv
import json
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (GERIZIM, deviation as beru_dev, tier_label as tier,
                       HARMONIC_STEP, BERU, KM_PER_DEGREE,
                       TIER_APP_KM, TIER_APLUS_KM, TIER_A_KM,
                       TIER_C_MAX, TIER_CMINUS, TIER_CMINUS2)
from lib.dome_filter import is_dome_site_raw, FORM_KEYWORD_RES

_KM_PER_BERU  = BERU * KM_PER_DEGREE
HARM_STEP_DEG = HARMONIC_STEP * BERU          # degrees per harmonic step (3.0°)
_KM_APP  = round(TIER_APP_KM)                 # ≈5 km
_KM_AP   = round(TIER_APLUS_KM)              # ≈11 km
_KM_A    = round(TIER_A_KM)                  # ≈21 km
_KM_C    = round(TIER_C_MAX   * _KM_PER_BERU)
_KM_CM   = round(TIER_CMINUS  * _KM_PER_BERU)
_KM_CMM  = round(TIER_CMINUS2 * _KM_PER_BERU)

OUTDIR = PROJECT_ROOT / "supplementary" / "interactive_map"
OUTDIR.mkdir(exist_ok=True)

WIKIDATA_CSV = PROJECT_ROOT / "data" / "store" / "unesco" / "wikidata_stupas_q180987.csv"
KEYWORDS_PATH = PROJECT_ROOT / "keywords.json"

SHOW_TIERS = {"A++", "A+", "A", "B", "C", "C-", "C--"}

SILK_ROAD_DIR = PROJECT_ROOT / "data" / "store" / "silk_road"


def _load_route_polylines(filename):
    """Return {route_name: [[lat, lon], ...]} from a curated polyline CSV.

    Rows carry route,order,lon,lat; grouped by route and sorted by order.
    Returns Leaflet-order [lat, lon] pairs.
    """
    import csv as _csv
    path = SILK_ROAD_DIR / filename
    by_route = {}
    if not path.exists():
        return by_route
    with open(path, newline="", encoding="utf-8") as fh:
        for row in _csv.DictReader(r for r in fh if not r.startswith("#")):
            try:
                r  = row["route"].strip()
                o  = int(row["order"])
                ln = float(row["lon"])
                lt = float(row["lat"])
            except (ValueError, KeyError):
                continue
            by_route.setdefault(r, []).append((o, lt, ln))
    return {k: [[lt, ln] for _, lt, ln in sorted(v)] for k, v in by_route.items()}


def _load_owtrad_edges():
    """Return list of [[lat1,lon1],[lat2,lon2]] edge pairs from OWTRAD route network."""
    import csv as _csv
    path = SILK_ROAD_DIR / "owtrad_routes.csv"
    edges = []
    if not path.exists():
        return edges
    with open(path, newline="", encoding="utf-8") as fh:
        for row in _csv.DictReader(r for r in fh if not r.startswith("#")):
            try:
                edges.append([
                    [float(row["lat1"]), float(row["lon1"])],
                    [float(row["lat2"]), float(row["lon2"])],
                ])
            except (ValueError, KeyError):
                pass
    return edges


def _load_silk_road_routes():
    """Return (overland_trunk, owtrad_edges, maritime_routes)."""
    overland = _load_route_polylines("overland_trunk.csv")
    maritime = _load_route_polylines("maritime_routes.csv")
    overland_list = [overland[k] for k in ("north", "south") if overland.get(k)]
    maritime_list = [maritime[k] for k in ("med", "sea") if maritime.get(k)]
    owtrad_edges  = _load_owtrad_edges()
    return overland_list, owtrad_edges, maritime_list


# ══════════════════════════════════════════════════════════════════════════════
#  Keyword extractors — driven entirely by keywords.json
# ══════════════════════════════════════════════════════════════════════════════

def _load_keyword_config():
    return json.loads(KEYWORDS_PATH.read_text())


def _build_founding_patterns(kw_cfg):
    sections = [
        "founding_capital", "sacred_origin", "founding_monument",
        "founding_axis", "ancient_landscape",
    ]
    patterns = {}
    for sec in sections:
        cfg = kw_cfg.get(sec, {})
        for kw in cfg.get("unambiguous", []) + cfg.get("ambiguous", []):
            if kw not in patterns:
                patterns[kw] = re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    return patterns


def _find_religions(text, religion_sets):
    return [
        rel for rel, kws in religion_sets.items()
        if any(re.search(r"\b" + re.escape(k) + r"\b", text, re.I) for k in kws)
    ]


def _find_dome_kws(text):
    return [kw for kw, rx in FORM_KEYWORD_RES.items() if rx.search(text)]


_MOUND_CHAIN_KEYS = ["mound_positive_context", "mound_unambiguous", "stupa", "dome"]
_MOUND_CHAIN_LABELS = ["mound", "tumulus/barrow", "stupa/dagoba", "dome/tholos"]


def _find_mound_stages(text, mound_evo):
    stages = []
    for label, key in zip(_MOUND_CHAIN_LABELS, _MOUND_CHAIN_KEYS):
        kws = mound_evo.get(key, [])
        if any(re.search(r"\b" + re.escape(k) + r"\b", text, re.I) for k in kws):
            stages.append(label)
    return stages


def _find_founding(text, founding_patterns, cap=4):
    return [kw for kw, pat in founding_patterns.items() if pat.search(text)][:cap]


# ══════════════════════════════════════════════════════════════════════════════
#  Data extraction
# ══════════════════════════════════════════════════════════════════════════════

def build_site_data():
    kw_cfg = _load_keyword_config()
    religion_sets = {k: v for k, v in kw_cfg.get("religion_sets", {}).items()
                     if isinstance(v, list)}
    mound_evo = kw_cfg.get("mound_evolution", {})
    founding_patterns = _build_founding_patterns(kw_cfg)

    print("Loading UNESCO corpus…")
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    print(f"  {len(cultural)} cultural/mixed sites with coordinates")

    sites = []

    for s in cultural:
        lat = getattr(s, "latitude", None)
        lon = getattr(s, "longitude", None)
        if lat is None or lon is None:
            continue
        d = beru_dev(lon)
        t = tier(d)
        if t not in SHOW_TIERS:
            continue
        full = s.full_text or ""
        sites.append({
            "name":         s.site,
            "lat":          lat,
            "lon":          lon,
            "tier":         t,
            "dev":          round(d, 5),
            "dome":         is_dome_site_raw(s),
            "corpus":       "UNESCO",
            "year":         getattr(s, "year", None),
            "religions":    _find_religions(full, religion_sets),
            "dome_kws":     _find_dome_kws(full),
            "mound_stages": _find_mound_stages(full, mound_evo),
            "founding":     _find_founding(full, founding_patterns),
        })

    print(f"  {len(sites)} UNESCO sites in shown tiers")

    print("Loading Wikidata Q180987 stupas…")
    wiki_count = 0
    with open(WIKIDATA_CSV, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(row for row in fh if not row.startswith("#"))
        for rec in reader:
            try:
                lat = float(rec["lat"])
                lon = float(rec["lon"])
            except (ValueError, KeyError):
                continue
            d = beru_dev(lon)
            t = tier(d)
            if t not in SHOW_TIERS:
                continue
            sites.append({
                "name":         rec.get("name", ""),
                "lat":          lat,
                "lon":          lon,
                "tier":         t,
                "dev":          round(d, 5),
                "dome":         True,
                "corpus":       "Wikidata Q180987",
                "year":         None,
                "religions":    ["Buddhism"],
                "dome_kws":     ["stupa"],
                "mound_stages": ["stupa/dagoba"],
                "founding":     [],
            })
            wiki_count += 1

    print(f"  {wiki_count} Wikidata sites in shown tiers")

    by_tier = {}
    for s in sites:
        by_tier[s["tier"]] = by_tier.get(s["tier"], 0) + 1
    print(f"  Tier breakdown: {dict(sorted(by_tier.items()))}")
    print(f"  Total: {len(sites)} sites")

    return sites


# ══════════════════════════════════════════════════════════════════════════════
#  HTML generation
# ══════════════════════════════════════════════════════════════════════════════

def build_html(sites, overland_routes=None, owtrad_edges=None, maritime_routes=None,
               midpoints=None):
    sites_json    = json.dumps(sites, ensure_ascii=False)
    overland_routes = overland_routes or []
    owtrad_edges    = owtrad_edges    or []
    maritime_routes = maritime_routes or []
    midpoints       = midpoints       or []

    html = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>Harmonic Tier Map — Gerizim Analysis</title>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: 'Georgia', serif; background: #111; }
#map { width: 100vw; height: 100vh; }

#panel {
  position: absolute; top: 12px; left: 58px; z-index: 1000;
  background: rgba(255,255,255,0.96); border-radius: 8px;
  padding: 10px 14px; font-size: 12px; max-width: 240px;
  box-shadow: 0 2px 12px rgba(0,0,0,0.35);
}
#panel-header {
  display: flex; align-items: center; justify-content: space-between; gap: 8px;
  margin-bottom: 5px;
}
#panel h2 { font-size: 13px; color: #222; margin: 0; }
#panel-toggle {
  background: none; border: none; cursor: pointer;
  font-size: 14px; color: #888; padding: 0; line-height: 1; flex-shrink: 0;
}
#panel-toggle:hover { color: #333; }
#panel-body { }
#panel.minimized #panel-body { display: none; }
#panel p  { color: #555; line-height: 1.4; font-size: 10.5px; margin-bottom: 8px; }

.tier-row, .rt-row {
  display: flex; align-items: center; gap: 6px; margin: 3px 0;
  cursor: pointer; user-select: none;
}
.tier-row input, .rt-row input { cursor: pointer; flex-shrink: 0; }
.tier-swatch { width: 13px; height: 13px; border-radius: 50%; flex-shrink: 0; }
.tier-label  { font-size: 11px; color: #333; flex: 1; }
.tier-count  { font-size: 10px; color: #888; }

.section-divider { margin-top: 8px; padding-top: 8px; border-top: 1px solid #ddd; }
.b-divider        { margin-top: 6px; padding-top: 6px; border-top: 1px solid #eee; }

.rt-line-solid  { width: 18px; height: 3px; background: #8e44ad; flex-shrink: 0; }
.rt-line-dashed { width: 18px; height: 0; border-top: 2.5px dashed #16a085; flex-shrink: 0; }
.rt-line-grid   { width: 18px; height: 0; border-top: 1.5px dashed #aaa; flex-shrink: 0; }
.rt-label { font-size: 11px; color: #333; }

#legend {
  position: absolute; bottom: 20px; right: 12px; z-index: 1000;
  background: rgba(255,255,255,0.95); border-radius: 8px;
  padding: 8px 12px; font-size: 11px;
  box-shadow: 0 2px 12px rgba(0,0,0,0.25);
  min-width: 170px;
}
#legend-header {
  display: flex; align-items: center; justify-content: space-between; gap: 8px;
}
#legend h3 { font-size: 11px; color: #333; margin: 0; }
#legend-header { display: flex; align-items: center; justify-content: space-between; gap: 8px; }
#legend-toggle {
  background: none; border: none; cursor: pointer;
  font-size: 14px; color: #888; padding: 0; line-height: 1; flex-shrink: 0;
}
#legend-toggle:hover { color: #333; }
#legend-body { margin-top: 6px; }
#legend.minimized #legend-body { display: none; }
.leg-row { display: flex; align-items: center; gap: 7px; margin: 4px 0; color: #333; }
.leg-sym { width: 18px; height: 18px; display: flex; align-items: center; justify-content: center; flex-shrink: 0; }

.leaflet-tooltip {
  background: rgba(15,15,30,0.93); color: #f0f0f0;
  border: 1px solid rgba(255,255,255,0.1); border-radius: 7px;
  font-size: 12px; line-height: 1.5; padding: 9px 12px;
  box-shadow: 0 3px 14px rgba(0,0,0,0.5);
  max-width: 280px; width: max-content; pointer-events: none;
  white-space: normal !important; word-break: break-word; overflow-wrap: break-word;
}
.leaflet-tooltip::before { display: none; }
.tt-name { font-weight: bold; font-size: 13px; margin-bottom: 4px; white-space: normal; word-break: break-word; }
.tt-tier { font-size: 11.5px; margin-bottom: 3px; }
.tt-row  { font-size: 10.5px; color: #ccc; margin: 2px 0; display: flex; gap: 5px; }
.tt-label { color: #888; white-space: nowrap; }
.tt-val   { color: #e8e8e8; }
.tt-chips { display: flex; flex-wrap: wrap; gap: 3px; margin-top: 4px; }
.chip { display: inline-block; padding: 1px 6px; border-radius: 10px; font-size: 10px; white-space: nowrap; }
.chip-rel   { background: rgba(52,152,219,0.3);  color: #90caf9; border: 1px solid rgba(52,152,219,0.4); }
.chip-dome  { background: rgba(230,126,34,0.25); color: #ffcc80; border: 1px solid rgba(230,126,34,0.35); }
.chip-mound { background: rgba(39,174,96,0.25);  color: #a5d6a7; border: 1px solid rgba(39,174,96,0.35); }
.chip-found { background: rgba(155,89,182,0.25); color: #ce93d8; border: 1px solid rgba(155,89,182,0.35); }

/* ── Draw tool ── */
#draw-btn {
  margin-top: 8px; width: 100%; padding: 5px 8px;
  font-size: 11px; font-family: inherit;
  background: #f0f4ff; border: 1px solid #c5d0e8; border-radius: 5px;
  cursor: pointer; color: #2c4a8e; text-align: center; transition: background 0.15s;
}
#draw-btn:hover { background: #dce6ff; }
#draw-btn.active { background: #fde8e8; border-color: #e74c3c; color: #c0392b; }

/* ── Result overlay ── */
#draw-result {
  display: none; position: absolute; top: 50%; left: 50%;
  transform: translate(-50%, -50%); z-index: 2000;
  background: rgba(15,15,35,0.97); color: #f0f0f0;
  border: 1px solid rgba(255,255,255,0.15); border-radius: 10px;
  box-shadow: 0 8px 32px rgba(0,0,0,0.6);
  min-width: 430px; max-width: 530px; font-size: 12px;
}
.dr-header {
  display: flex; align-items: center; justify-content: space-between;
  padding: 10px 14px 8px; border-bottom: 1px solid rgba(255,255,255,0.1);
}
.dr-title { font-size: 13px; font-weight: bold; color: #fff; }
.dr-close { background: none; border: none; color: #aaa; font-size: 16px; cursor: pointer; padding: 0 2px; line-height: 1; }
.dr-close:hover { color: #fff; }
.dr-body { padding: 10px 14px 12px; }
.dr-section { font-size: 10.5px; color: #aaa; margin: 8px 0 5px; font-style: italic; }
.dr-section:first-child { margin-top: 0; }
.dr-table { width: 100%; border-collapse: collapse; font-size: 11px; }
.dr-table th { text-align: left; color: #999; font-weight: normal; font-size: 10px; padding: 2px 6px 4px; border-bottom: 1px solid rgba(255,255,255,0.1); }
.dr-table td { padding: 3px 6px; vertical-align: middle; }
.dr-table tbody tr:hover { background: rgba(255,255,255,0.04); }
.dr-table .gap-row td { padding: 2px 0; border-top: 1px solid rgba(255,255,255,0.07); }
.tier-badge { padding: 1px 6px; border-radius: 8px; color: #fff; font-size: 10px; font-weight: bold; white-space: nowrap; display: inline-block; }
.sig  { color: #f1c40f; font-weight: bold; font-size: 11px; }
.ns   { color: #555; font-size: 11px; }
.dr-note { margin-top: 10px; font-size: 10px; color: #555; font-style: italic; }
</style>
</head>
<body>
<div id="map"></div>

<div id="panel">
  <div id="panel-header">
    <h2>Harmonic Tier Map</h2>
    <button id="panel-toggle" title="Minimize / maximize">&#8722;</button>
  </div>
  <div id="panel-body">
  <p>UNESCO cultural + Wikidata Q180987 stupas graded by b&#275;ru deviation from Gerizim. Hover for details.</p>
  <div id="tier-controls">
    <div class="tier-row" data-tier="A++">
      <input type="checkbox" checked/>
      <span class="tier-swatch" style="background:#922b21"></span>
      <span class="tier-label">A<sup>++</sup> &#8804;__KM_APP__ km</span>
      <span class="tier-count" id="cnt-APP"></span>
    </div>
    <div class="tier-row" data-tier="A+">
      <input type="checkbox" checked/>
      <span class="tier-swatch" style="background:#c0392b"></span>
      <span class="tier-label">A<sup>+</sup> &#8804;__KM_AP__ km</span>
      <span class="tier-count" id="cnt-AP"></span>
    </div>
    <div class="tier-row" data-tier="A">
      <input type="checkbox" checked/>
      <span class="tier-swatch" style="background:#e67e22"></span>
      <span class="tier-label">A &#8804;__KM_A__ km</span>
      <span class="tier-count" id="cnt-A"></span>
    </div>
    <div class="tier-row" data-tier="C">
      <input type="checkbox" checked/>
      <span class="tier-swatch" style="background:#7fb3d3"></span>
      <span class="tier-label">C &#8804;__KM_C__ km (mid)</span>
      <span class="tier-count" id="cnt-C"></span>
    </div>
    <div class="tier-row" data-tier="C-">
      <input type="checkbox" checked/>
      <span class="tier-swatch" style="background:#2471a3"></span>
      <span class="tier-label">C<sup>&#8722;</sup> &#8804;__KM_CM__ km (mid)</span>
      <span class="tier-count" id="cnt-CM"></span>
    </div>
    <div class="tier-row" data-tier="C--">
      <input type="checkbox" checked/>
      <span class="tier-swatch" style="background:#1a3a6b"></span>
      <span class="tier-label">C<sup>&#8722;&#8722;</sup> &#8804;__KM_CMM__ km (mid)</span>
      <span class="tier-count" id="cnt-CMM"></span>
    </div>
    <div class="tier-row b-divider" data-tier="B">
      <input type="checkbox"/>
      <span class="tier-swatch" style="background:#aaaaaa"></span>
      <span class="tier-label">B (remainder)</span>
      <span class="tier-count" id="cnt-B"></span>
    </div>
  </div>
  <div class="section-divider">
    <div style="font-size:10px;color:#888;margin-bottom:4px">Corpus</div>
    <div class="rt-row" id="toggle-unesco">
      <input type="checkbox" checked/>
      <span class="rt-label">UNESCO cultural/mixed</span>
      <span class="tier-count" id="cnt-corpus-UNESCO"></span>
    </div>
    <div class="rt-row" id="toggle-wikidata">
      <input type="checkbox" checked/>
      <span class="rt-label">Wikidata Q180987</span>
      <span class="tier-count" id="cnt-corpus-Wikidata"></span>
    </div>
  </div>
  <div class="section-divider">
    <div class="rt-row" id="toggle-overland">
      <input type="checkbox" checked/>
      <div class="rt-line-solid"></div>
      <span class="rt-label">Overland Silk Road</span>
    </div>
    <div class="rt-row" id="toggle-maritime">
      <input type="checkbox" checked/>
      <div class="rt-line-dashed"></div>
      <span class="rt-label">Maritime Silk Road</span>
    </div>
  </div>
  <div class="section-divider">
    <div class="rt-row" id="toggle-grid">
      <input type="checkbox" checked/>
      <div class="rt-line-grid"></div>
      <span class="rt-label">0.1-b&#275;ru harmonic grid</span>
    </div>
  </div>
  <div class="section-divider">
    <div class="rt-row" id="toggle-midpoints">
      <input type="checkbox"/>
      <span style="width:13px;height:13px;border-radius:50%;background:#7d3c98;flex-shrink:0;display:inline-block"></span>
      <span class="rt-label">OWTRAD edge midpoints</span>
    </div>
  </div>
  <div class="section-divider">
    <button id="draw-btn">&#11042; Draw area</button>
  </div>
  </div><!-- end panel-body -->
</div>

<div id="draw-result"></div>

<div id="legend">
  <div id="legend-header">
    <h3>Legend</h3>
    <button id="legend-toggle" title="Minimize / maximize">&#8722;</button>
  </div>
  <div id="legend-body">
    <div class="leg-row">
      <div class="leg-sym"><svg width="14" height="14"><circle cx="7" cy="7" r="6" fill="#888" stroke="white" stroke-width="1.5"/></svg></div>
      &#9679; UNESCO dome/stupa
    </div>
    <div class="leg-row">
      <div class="leg-sym"><svg width="14" height="14"><rect x="1" y="1" width="12" height="12" fill="#888" stroke="white" stroke-width="1.5"/></svg></div>
      &#9632; UNESCO non-dome
    </div>
    <div class="leg-row">
      <div class="leg-sym"><svg width="14" height="14"><polygon points="7,1 13,13 1,13" fill="#888" stroke="white" stroke-width="1.5"/></svg></div>
      &#9650; Wikidata Q180987
    </div>
    <div style="margin-top:8px;padding-top:6px;border-top:1px solid #ddd;font-size:10px;color:#777;line-height:1.9">
      <span style="background:rgba(52,152,219,0.3);color:#3498db;padding:0 5px;border-radius:8px">&#9632;</span> Religion &nbsp;
      <span style="background:rgba(230,126,34,0.25);color:#e67e22;padding:0 5px;border-radius:8px">&#9632;</span> Dome type<br/>
      <span style="background:rgba(39,174,96,0.25);color:#27ae60;padding:0 5px;border-radius:8px">&#9632;</span> Mound evo &nbsp;
      <span style="background:rgba(155,89,182,0.25);color:#9b59b6;padding:0 5px;border-radius:8px">&#9632;</span> Founding
    </div>
  </div>
</div>

<script>
const SITES = """ + sites_json + """;
const GERIZIM = """ + str(GERIZIM) + """;
const HARM_STEP_DEG = """ + str(HARM_STEP_DEG) + """;
const OVERLAND_ROUTES = """ + json.dumps(overland_routes) + """;
const OWTRAD_EDGES    = """ + json.dumps(owtrad_edges) + """;
const MARITIME_ROUTES = """ + json.dumps(maritime_routes) + """;
const MIDPOINTS       = """ + json.dumps(midpoints) + """;
const OWTRAD_ZOOM_THRESHOLD = 5;

const TIER_COLOR  = { "A++":"#922b21","A+":"#c0392b","A":"#e67e22","C":"#7fb3d3","C-":"#2471a3","C--":"#1a3a6b","B":"#aaaaaa" };
const TIER_RADIUS = { "A++":9,"A+":7,"A":5,"C":5,"C-":7,"C--":9,"B":4 };

const WORLD_BOUNDS = L.latLngBounds(L.latLng(-85,-180), L.latLng(85,180));
const map = L.map("map", {
  center: [25, 50], zoom: 3, minZoom: 2, maxZoom: 18,
  maxBounds: WORLD_BOUNDS, maxBoundsViscosity: 1.0, worldCopyJump: false,
});
L.tileLayer("https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png", {
  attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OSM</a> &copy; <a href="https://carto.com/">CARTO</a>',
  subdomains: "abcd", maxZoom: 18, noWrap: true, bounds: WORLD_BOUNDS,
}).addTo(map);

const gridLayer = L.layerGroup();
const GRID_STEPS = Math.ceil(180 / HARM_STEP_DEG);
for (let k = -GRID_STEPS; k <= GRID_STEPS; k++) {
  const lon = GERIZIM + k * HARM_STEP_DEG;
  if (lon < -180 || lon > 180) continue;
  L.polyline([[-85,lon],[85,lon]], {
    color: k===0 ? "#444444" : "#888888",
    weight: k===0 ? 1.4 : 0.8,
    opacity: k===0 ? 0.85 : 0.65,
    dashArray: k===0 ? null : "4 7",
  }).addTo(gridLayer);
}
gridLayer.addTo(map);

// Overland: trunk (low zoom) swaps to dense OWTRAD edges at OWTRAD_ZOOM_THRESHOLD
const trunkLayer  = L.layerGroup(
  OVERLAND_ROUTES.map(pts => L.polyline(pts, { color:"#8e44ad", weight:2.2, opacity:0.65 }))
);
const owtradLayer = L.layerGroup(
  OWTRAD_EDGES.map(([a,b]) => L.polyline([a,b], { color:"#8e44ad", weight:1.0, opacity:0.45 }))
);
const maritimeLayer = L.layerGroup(
  MARITIME_ROUTES.map(pts => L.polyline(pts, { color:"#16a085", weight:2.2, opacity:0.65, dashArray:"8 5" }))
);

let overlandOn = true;
function _applyOverlandZoom() {
  if (!overlandOn) return;
  const detailed = map.getZoom() >= OWTRAD_ZOOM_THRESHOLD;
  if (detailed) {
    if (map.hasLayer(trunkLayer))  map.removeLayer(trunkLayer);
    if (!map.hasLayer(owtradLayer)) owtradLayer.addTo(map);
  } else {
    if (map.hasLayer(owtradLayer)) map.removeLayer(owtradLayer);
    if (!map.hasLayer(trunkLayer))  trunkLayer.addTo(map);
  }
}
_applyOverlandZoom();
map.on("zoomend", _applyOverlandZoom);
maritimeLayer.addTo(map);

// OWTRAD edge midpoints — purple fill, tier-coloured ring, hidden by default
const TIER_COLOR_MID = { "A++":"#922b21","A+":"#c0392b","A":"#e67e22","C":"#7fb3d3","C-":"#2471a3","C--":"#1a3a6b","B":"#aaaaaa" };
const midpointLayer = L.layerGroup();
MIDPOINTS.forEach(m => {
  const r = 5, size = r*2+6, c = r+3;
  const ring = TIER_COLOR_MID[m.tier] || "#aaaaaa";
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${size}" height="${size}"><circle cx="${c}" cy="${c}" r="${r}" fill="#7d3c98" stroke="${ring}" stroke-width="2"/></svg>`;
  const icon = L.divIcon({ html: svg, iconSize:[size,size], iconAnchor:[c,c], className:"" });
  L.marker([m.mid_lat, m.mid_lon], { icon })
   .bindTooltip(`<div class="tt-name">Midpoint: ${m.node1} → ${m.node2}</div>
     <div class="tt-tier">Tier <strong style="color:#fff">${m.tier}</strong> &middot; &delta; = ${m.dev}</div>
     <div class="tt-row"><span class="tt-label">Mid lon:</span><span class="tt-val">${m.mid_lon.toFixed(3)}°E</span></div>
     <div class="tt-row"><span class="tt-label">Dataset:</span><span class="tt-val">${m.dataset}</span></div>`,
     { sticky:true, opacity:1, direction:"top" })
   .addTo(midpointLayer);
});

function makeIcon(shape, color, r) {
  const size = r * 2 + 5, c = r + 2.5;
  let body;
  if (shape === "circle") {
    body = `<circle cx="${c}" cy="${c}" r="${r}" fill="${color}" stroke="white" stroke-width="1.5"/>`;
  } else if (shape === "square") {
    body = `<rect x="2.5" y="2.5" width="${r*2}" height="${r*2}" fill="${color}" stroke="white" stroke-width="1.5"/>`;
  } else {
    body = `<polygon points="${c},2 ${c+r},${size-2} ${c-r},${size-2}" fill="${color}" stroke="white" stroke-width="1.5"/>`;
  }
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${size}" height="${size}">${body}</svg>`;
  return L.divIcon({ html: svg, iconSize:[size,size], iconAnchor:[c,c], className:"" });
}

const MOUND_ORDER = ["mound","tumulus/barrow","stupa/dagoba","dome/tholos"];

function buildTooltip(s) {
  const name  = s.name.replace(/<[^>]*>/g,"").replace(/&amp;/g,"&").replace(/&lt;/g,"<").replace(/&gt;/g,">");
  const shape = s.corpus === "Wikidata Q180987" ? "stupa (Wikidata)"
              : s.dome ? "dome/stupa" : "non-dome cultural";
  const moundOrdered = MOUND_ORDER.filter(m => (s.mound_stages||[]).includes(m));
  const moundHtml  = moundOrdered.map(m=>`<span class="chip chip-mound">${m}</span>`).join(" \u2192 ");
  const relChips   = (s.religions||[]).map(r=>`<span class="chip chip-rel">${r}</span>`).join(" ");
  const domeChips  = (s.dome_kws||[]).map(d=>`<span class="chip chip-dome">${d}</span>`).join(" ");
  const foundChips = (s.founding||[]).map(f=>`<span class="chip chip-found">${f}</span>`).join(" ");
  return `
    <div class="tt-name">${name}</div>
    <div class="tt-tier">Tier <strong style="color:#ffffff">${s.tier}</strong>
      &nbsp;&middot;&nbsp; \u03b4 = ${s.dev} &nbsp;&middot;&nbsp; ${s.corpus}</div>
    <div class="tt-row"><span class="tt-label">Shape:</span><span class="tt-val">${shape}</span></div>
    <div class="tt-row"><span class="tt-label">Inscribed:</span><span class="tt-val">${s.year||"\u2014"}</span></div>
    ${relChips   ? `<div class="tt-chips">${relChips}</div>` : ""}
    ${domeChips  ? `<div class="tt-chips">${domeChips}</div>` : ""}
    ${moundHtml  ? `<div class="tt-chips">${moundHtml}</div>` : ""}
    ${foundChips ? `<div class="tt-chips">${foundChips}</div>` : ""}
  `;
}

const allMarkers = []; // {marker, tier, corpus}
const visibleTiers   = new Set(Object.keys(TIER_COLOR).filter(t => t !== "B"));
const visibleCorpora = new Set(["UNESCO", "Wikidata Q180987"]);

function applyFilters() {
  allMarkers.forEach(({marker, tier, corpus}) => {
    const show = visibleTiers.has(tier) && visibleCorpora.has(corpus);
    if (show && !map.hasLayer(marker)) marker.addTo(map);
    else if (!show && map.hasLayer(marker)) map.removeLayer(marker);
  });
}

const tierCounts   = {};
const corpusCounts = {};

SITES.forEach(s => {
  const color  = TIER_COLOR[s.tier] || "#aaaaaa";
  const r      = TIER_RADIUS[s.tier] || 4;
  const shape  = s.corpus === "Wikidata Q180987" ? "triangle"
               : s.dome ? "circle" : "square";
  const marker = L.marker([s.lat, s.lon], { icon: makeIcon(shape, color, r) });
  marker.bindTooltip(buildTooltip(s), { sticky: true, opacity: 1, direction: "top" });
  allMarkers.push({marker, tier: s.tier, corpus: s.corpus});
  tierCounts[s.tier]     = (tierCounts[s.tier]    ||0) + 1;
  corpusCounts[s.corpus] = (corpusCounts[s.corpus]||0) + 1;
});

applyFilters();

const cntIds = {"A++":"APP","A+":"AP","A":"A","C":"C","C-":"CM","C--":"CMM","B":"B"};
Object.entries(tierCounts).forEach(([t,n]) => {
  const el = document.getElementById("cnt-"+cntIds[t]);
  if (el) el.textContent = "("+n+")";
});
["UNESCO","Wikidata Q180987"].forEach(corp => {
  const id = "cnt-corpus-" + (corp === "Wikidata Q180987" ? "Wikidata" : corp);
  const el = document.getElementById(id);
  if (el) el.textContent = "(" + (corpusCounts[corp]||0) + ")";
});

document.querySelectorAll("#tier-controls .tier-row").forEach(row => {
  const cb = row.querySelector("input");
  const t  = row.dataset.tier;
  row.addEventListener("click", e => {
    if (e.target !== cb) cb.checked = !cb.checked;
    cb.checked ? visibleTiers.add(t) : visibleTiers.delete(t);
    applyFilters();
  });
});

function bindCorpusToggle(id, corpusKey) {
  document.getElementById(id).addEventListener("click", function(e) {
    const cb = this.querySelector("input");
    if (e.target !== cb) cb.checked = !cb.checked;
    cb.checked ? visibleCorpora.add(corpusKey) : visibleCorpora.delete(corpusKey);
    applyFilters();
  });
}
bindCorpusToggle("toggle-unesco",   "UNESCO");
bindCorpusToggle("toggle-wikidata", "Wikidata Q180987");

function bindToggle(id, layer) {
  document.getElementById(id).addEventListener("click", function(e) {
    const cb = this.querySelector("input");
    if (e.target !== cb) cb.checked = !cb.checked;
    cb.checked ? layer.addTo(map) : map.removeLayer(layer);
  });
}
document.getElementById("toggle-overland").addEventListener("click", function(e) {
  const cb = this.querySelector("input");
  if (e.target !== cb) cb.checked = !cb.checked;
  overlandOn = cb.checked;
  if (!overlandOn) {
    map.removeLayer(trunkLayer);
    map.removeLayer(owtradLayer);
  } else {
    _applyOverlandZoom();
  }
});
bindToggle("toggle-maritime",   maritimeLayer);
bindToggle("toggle-grid",       gridLayer);
bindToggle("toggle-midpoints",  midpointLayer);

// ── Legend minimize / maximize ────────────────────────────────────────────────
document.getElementById("legend-toggle").addEventListener("click", () => {
  const leg = document.getElementById("legend");
  const btn = document.getElementById("legend-toggle");
  const minimized = leg.classList.toggle("minimized");
  btn.textContent = minimized ? "+" : "\u2212";
  btn.title = minimized ? "Maximize" : "Minimize";
});

// ── Panel minimize / maximize ─────────────────────────────────────────────────
document.getElementById("panel-toggle").addEventListener("click", () => {
  const panel = document.getElementById("panel");
  const btn   = document.getElementById("panel-toggle");
  const minimized = panel.classList.toggle("minimized");
  btn.textContent = minimized ? "+" : "\u2212";
  btn.title = minimized ? "Maximize" : "Minimize";
});

// ── Statistics ────────────────────────────────────────────────────────────────
function lgamma(z) {
  // Lanczos approximation, numerically stable for z > 0
  const c = [0.99999999999980993,676.5203681218851,-1259.1392167224028,
    771.32342877765313,-176.61502916214059,12.507343278686905,
    -0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
  if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z);
  z -= 1;
  let x = c[0];
  for (let i = 1; i < 9; i++) x += c[i] / (z + i);
  const t = z + 7.5;
  return 0.5*Math.log(2*Math.PI) + (z+0.5)*Math.log(t) - t + Math.log(x);
}

function logBinom(n, k) {
  if (k < 0 || k > n || n < 0) return -Infinity;
  if (k === 0 || k === n) return 0;
  return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}

// One-tailed hypergeometric P(X >= k)
function hypergeomSF(N, K, n, k) {
  let p = 0;
  const lo = Math.max(k, Math.max(0, n-(N-K)));
  const hi = Math.min(n, K);
  for (let i = lo; i <= hi; i++) {
    const lp = logBinom(K,i) + logBinom(N-K,n-i) - logBinom(N,n);
    if (isFinite(lp)) p += Math.exp(lp);
  }
  return Math.min(1, p);
}

// Fisher's exact one-tailed (enrichment in selection)
function fisherExact(k, n, K, N) {
  // k = tier hits in selection, n = selection size, K = tier hits globally, N = corpus size
  return hypergeomSF(N, K, n, k);
}

// Multivariate hypergeometric joint: P(X_A >= kA AND X_C >= kC)
// Tiers are mutually exclusive — sites belong to exactly one tier
function jointHypergeom(N, KA, KC, n, kA, kC) {
  if (kA === 0 && kC === 0) return 1;
  const KO = N - KA - KC;
  let p = 0;
  const hiA = Math.min(n, KA);
  const hiC = Math.min(n, KC);
  for (let i = kA; i <= hiA; i++) {
    for (let j = kC; j <= hiC; j++) {
      const rem = n - i - j;
      if (rem < 0 || rem > KO) continue;
      const lp = logBinom(KA,i) + logBinom(KC,j) + logBinom(KO,rem) - logBinom(N,n);
      if (isFinite(lp)) p += Math.exp(lp);
    }
  }
  return Math.min(1, p);
}

function sigStars(p) {
  if (p < 0.001) return "***";
  if (p < 0.01)  return "**";
  if (p < 0.05)  return "*";
  return "ns";
}
function fmtP(p) {
  if (p === 1) return "1.0000";
  if (p < 0.0001) return "< 0.0001";
  return p.toFixed(4);
}

// ── Point-in-polygon (ray casting) ───────────────────────────────────────────
function pointInPoly(lat, lng, poly) {
  let inside = false;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
    const yi = poly[i].lat, xi = poly[i].lng;
    const yj = poly[j].lat, xj = poly[j].lng;
    if (((yi > lat) !== (yj > lat)) && (lng < (xj-xi)*(lat-yi)/(yj-yi)+xi))
      inside = !inside;
  }
  return inside;
}

// ── Analysis ──────────────────────────────────────────────────────────────────
const TIER_BADGE_COLOR = TIER_COLOR;
const A_TIERS = ["A++","A+","A"];
const C_TIERS = ["C","C-","C--"];

function analyzeSelection(poly) {
  const activeSites = SITES.filter(s => visibleTiers.has(s.tier) && visibleCorpora.has(s.corpus));
  const N = activeSites.length;
  const inSel = activeSites.filter(s => pointInPoly(s.lat, s.lon, poly));
  const n = inSel.length;
  if (n === 0) { alert("No sites found inside the drawn area."); return; }

  const selCnt = {}, globCnt = {};
  activeSites.forEach(s => { globCnt[s.tier] = (globCnt[s.tier]||0)+1; });
  inSel.forEach(s => { selCnt[s.tier]  = (selCnt[s.tier] ||0)+1; });

  // Fisher for each tier
  const fisherRes = {};
  [...A_TIERS, ...C_TIERS].forEach(t => {
    const k = selCnt[t]||0, K = globCnt[t]||0;
    fisherRes[t] = { k, K, p: fisherExact(k, n, K, N) };
  });

  // Joint hypergeometric pairs
  const selA  = A_TIERS.reduce((s,t)=>s+(selCnt[t]||0), 0);
  const selC  = C_TIERS.reduce((s,t)=>s+(selCnt[t]||0), 0);
  const globA = A_TIERS.reduce((s,t)=>s+(globCnt[t]||0), 0);
  const globC = C_TIERS.reduce((s,t)=>s+(globCnt[t]||0), 0);

  const jointPairs = [
    { label: "A\u207A\u207A &#x2229; C\u207B\u207B", kA: selCnt["A++"]||0, kC: selCnt["C--"]||0, KA: globCnt["A++"]||0, KC: globCnt["C--"]||0 },
    { label: "A\u207A &#x2229; C\u207B",              kA: selCnt["A+"] ||0, kC: selCnt["C-"] ||0, KA: globCnt["A+"] ||0, KC: globCnt["C-"] ||0 },
    { label: "A &#x2229; C",                          kA: selCnt["A"]  ||0, kC: selCnt["C"]  ||0, KA: globCnt["A"]  ||0, KC: globCnt["C"]  ||0 },
    { label: "All-A &#x2229; All-C",                  kA: selA, kC: selC, KA: globA, KC: globC },
  ];
  jointPairs.forEach(r => { r.p = jointHypergeom(N, r.KA, r.KC, n, r.kA, r.kC); });

  showResult(n, N, fisherRes, jointPairs);
}

function showResult(n, N, fisherRes, jointPairs) {
  const el = document.getElementById("draw-result");

  function starSpan(p) {
    const s = sigStars(p);
    return s === "ns" ? `<span class="ns">ns</span>` : `<span class="sig">${s}</span>`;
  }
  function fisherRow(t) {
    const {k, K, p} = fisherRes[t];
    return `<tr>
      <td><span class="tier-badge" style="background:${TIER_BADGE_COLOR[t]}">${t}</span></td>
      <td style="text-align:right">${k} / ${n}</td>
      <td style="text-align:right">${K} / ${N}</td>
      <td style="text-align:right">${fmtP(p)} ${starSpan(p)}</td>
    </tr>`;
  }
  function jointRow(r) {
    return `<tr>
      <td>${r.label}</td>
      <td style="text-align:right">${r.kA}</td>
      <td style="text-align:right">${r.kC}</td>
      <td style="text-align:right">${fmtP(r.p)} ${starSpan(r.p)}</td>
    </tr>`;
  }

  el.innerHTML = `
    <div class="dr-header">
      <span class="dr-title">&#11042; Area selection &mdash; ${n} of ${N} sites</span>
      <button class="dr-close" onclick="resetDraw()">&#10005;</button>
    </div>
    <div class="dr-body">
      <div class="dr-section">Fisher&rsquo;s exact &mdash; tier enrichment in selection vs. global corpus (one-tailed)</div>
      <table class="dr-table">
        <thead><tr><th>Tier</th><th style="text-align:right">In area</th><th style="text-align:right">Global</th><th style="text-align:right">p</th></tr></thead>
        <tbody>
          ${A_TIERS.map(fisherRow).join("")}
          <tr class="gap-row"><td colspan="4"></td></tr>
          ${C_TIERS.map(fisherRow).join("")}
        </tbody>
      </table>
      <div class="dr-section">Hypergeometric joint &mdash; P(A&#x209C; &ge; k<sub>A</sub> &cap; C&#x209C; &ge; k<sub>C</sub>) in same selection</div>
      <table class="dr-table">
        <thead><tr><th>Pair</th><th style="text-align:right">k<sub>A</sub></th><th style="text-align:right">k<sub>C</sub></th><th style="text-align:right">p (joint)</th></tr></thead>
        <tbody>${jointPairs.map(jointRow).join("")}</tbody>
      </table>
      <div class="dr-note">&#x2605;&#x2605;&#x2605; p&lt;0.001 &nbsp; &#x2605;&#x2605; p&lt;0.01 &nbsp; &#x2605; p&lt;0.05 &nbsp;&mdash;&nbsp; Click outside to reset.</div>
    </div>`;
  el.style.display = "block";
}

// ── Freehand lasso draw tool ──────────────────────────────────────────────────
let drawMode = false, lassoPts = [], lassoLine = null, lassoPolygon = null;
const drawBtn = document.getElementById("draw-btn");

function resetDraw() {
  lassoPts = [];
  if (lassoLine)    { map.removeLayer(lassoLine);    lassoLine    = null; }
  if (lassoPolygon) { map.removeLayer(lassoPolygon); lassoPolygon = null; }
  document.getElementById("draw-result").style.display = "none";
  drawMode = false;
  drawBtn.classList.remove("active");
  drawBtn.textContent = "\u2B22 Draw area";
  map.getContainer().style.cursor = "";
  map.dragging.enable();
}

drawBtn.addEventListener("click", () => {
  drawMode = !drawMode;
  drawBtn.classList.toggle("active", drawMode);
  drawBtn.textContent = drawMode ? "\u2715 Cancel draw" : "\u2B22 Draw area";
  map.getContainer().style.cursor = drawMode ? "crosshair" : "";
  if (!drawMode) resetDraw();
});

// Click outside result box resets everything
document.addEventListener("mousedown", e => {
  const el = document.getElementById("draw-result");
  if (el.style.display === "block" && !el.contains(e.target) && e.target !== drawBtn)
    resetDraw();
}, true);

// Lasso drawing — intercept at container level to suppress map drag
// Skip if click originated inside the panel or result overlay
map.getContainer().addEventListener("mousedown", e => {
  if (!drawMode) return;
  if (document.getElementById("panel").contains(e.target)) return;
  if (document.getElementById("draw-result").contains(e.target)) return;
  map.dragging.disable();
  const latlng = map.mouseEventToLatLng(e);
  lassoPts = [latlng];
  if (lassoPolygon) { map.removeLayer(lassoPolygon); lassoPolygon = null; }
  if (lassoLine)    { map.removeLayer(lassoLine);    lassoLine    = null; }
}, true);

map.getContainer().addEventListener("mousemove", e => {
  if (!drawMode || lassoPts.length === 0) return;
  if (document.getElementById("panel").contains(e.target)) return;
  const latlng = map.mouseEventToLatLng(e);
  lassoPts.push(latlng);
  if (lassoLine) map.removeLayer(lassoLine);
  lassoLine = L.polyline(lassoPts, { color:"#e74c3c", weight:2, dashArray:"4 3", opacity:0.9 }).addTo(map);
}, true);

map.getContainer().addEventListener("mouseup", e => {
  if (!drawMode || lassoPts.length < 6) return;
  map.dragging.enable();
  if (lassoLine) { map.removeLayer(lassoLine); lassoLine = null; }
  lassoPolygon = L.polygon(lassoPts, {
    color:"#e74c3c", weight:1.5, fillColor:"#e74c3c", fillOpacity:0.07
  }).addTo(map);
  analyzeSelection(lassoPts);
  drawMode = false;
  drawBtn.classList.remove("active");
  drawBtn.textContent = "\u2B22 Draw area";
  map.getContainer().style.cursor = "";
}, true);
</script>
</body>
</html>
"""
    html = (html
        .replace("__KM_APP__",  str(_KM_APP))
        .replace("__KM_AP__",   str(_KM_AP))
        .replace("__KM_A__",    str(_KM_A))
        .replace("__KM_C__",    str(_KM_C))
        .replace("__KM_CM__",   str(_KM_CM))
        .replace("__KM_CMM__",  str(_KM_CMM))
    )
    return html


# ══════════════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 60)
    print("  generate_interactive_map.py")
    print("=" * 60)
    print()

    sites = build_site_data()
    print()

    print("Loading Silk Road routes…")
    overland_routes, owtrad_edges, maritime_routes = _load_silk_road_routes()
    print(f"  {len(overland_routes)} overland trunk polyline(s), {len(owtrad_edges)} OWTRAD edges, {len(maritime_routes)} maritime polyline(s)")

    print("Building edge midpoints…")
    import csv as _csv
    midpoints = []
    with open(SILK_ROAD_DIR / "owtrad_routes.csv", newline="", encoding="utf-8") as fh:
        for row in _csv.DictReader(r for r in fh if not r.startswith("#")):
            try:
                lon1 = float(row["lon1"]); lat1 = float(row["lat1"])
                lon2 = float(row["lon2"]); lat2 = float(row["lat2"])
            except (ValueError, KeyError):
                continue
            mid_lon = (lon1 + lon2) / 2
            mid_lat = (lat1 + lat2) / 2
            from lib.beru import deviation as _dev, tier_label as _tier
            d = _dev(mid_lon)
            midpoints.append({
                "mid_lat": round(mid_lat, 5),
                "mid_lon": round(mid_lon, 5),
                "node1":   row.get("node1", ""),
                "node2":   row.get("node2", ""),
                "dataset": row.get("dataset", ""),
                "tier":    _tier(d),
                "dev":     round(d, 5),
            })
    print(f"  {len(midpoints)} edge midpoints")
    print()

    print("Building HTML…")
    html = build_html(sites, overland_routes, owtrad_edges, maritime_routes, midpoints)

    out = OUTDIR / "silkroad_interactive.html"
    out.write_text(html, encoding="utf-8")
    print(f"  Written: {out}  ({len(html):,} bytes)")
    print()
    print("Open in browser:")
    print(f"  open {out}")
