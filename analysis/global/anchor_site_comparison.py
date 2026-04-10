"""
anchor_site_comparison.py
=========================
Computes beru deviations for individually named sites from BOTH the
Gerizim and Jerusalem anchors.  Produces the LaTeX macro block for
GROUP 16 (Lumbini, Takht-e Soleyman, Khoja Ahmed Yasawi) and
GROUP 22 (World Peace Pagoda) in paper_a_primary_unesco.tex.

Every macro printed by this script maps 1-to-1 to a \\newcommand in
the manuscript.  If any value changes, re-run this script and update
the .tex file.

USAGE
-----
  python3 analysis/global/anchor_site_comparison.py

OUTPUT
------
  Prints LaTeX macro definitions to stdout.
  GROUP 16: Anchor-comparison macros (Lumbini, Takht, Khoja, NappTotal)
  GROUP 22: World Peace Pagoda macros (WPPDeltaKm, WPPDeltaBeru)
"""

from __future__ import annotations

import sys
from pathlib import Path
import json as _json
import numpy as np

# ── Repo root ─────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, TIER_APP, TIER_APLUS

# ── Constants from config.json ────────────────────────────────────────────────
_CFG      = _json.loads((ROOT / "config.json").read_text())
JERUSALEM = _CFG["anchors"]["jerusalem"]["longitude"]    # 35.2317 — from config
KM_PER_DEG = _CFG["units"]["km_per_degree"]             # 111.0

# World Peace Pagoda, Lumbini — Wikidata Q6540965 coordinate (not in UNESCO XML)
# The UNESCO XML records the Lumbini WHS centroid at 83.27611°E (used for
# \LumbiniDevKm etc above).  The Pagoda itself has a distinct Wikidata entry
# at 83.2743°E / 27.4833°N.  These are two different points ~0.9 km apart.
# Beru deviation is computed with the equatorial convention (no cos-lat),
# matching every other site in the study.
WPP_LON = 83.2743            # Wikidata Q6540965 — Nipponzan-Myohoji pagoda


def beru_dev(lon: float, anchor: float) -> tuple[float, float, str]:
    """Return (dev_beru, dev_km, tier_label) for a site relative to anchor."""
    arc  = abs(lon - anchor)
    bv   = arc / BERU
    near = round(bv * 10) / 10
    dev  = abs(bv - near)
    km   = dev * BERU * KM_PER_DEG   # equatorial, no cos-lat
    if dev <= TIER_APP:
        tier = "A++"
    elif dev <= TIER_APLUS:
        tier = "A+"
    else:
        tier = "A" if dev <= 0.010 else "B"
    return dev, km, tier


def find_site(sites: list[dict], fragment: str) -> dict | None:
    frag = fragment.lower()
    for s in sites:
        if frag in s["name"].lower():
            return s
    return None


def main():
    # Load UNESCO corpus
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    sites = [{"name": s.site, "lon": s.longitude, "lat": s.latitude}
             for s in cultural]
    lons = np.array([s["lon"] for s in sites])

    # ── Count A++ from Gerizim ────────────────────────────────────────────
    arcs = np.abs(lons - GERIZIM)
    arcs = np.minimum(arcs, 360 - arcs)
    bvs  = arcs / BERU
    near = np.round(bvs * 10) / 10
    devs = np.abs(bvs - near)
    n_app = int(np.sum(devs <= TIER_APP))

    print("=" * 72)
    print("  GROUP 16 — Anchor-Comparison Macros")
    print("  Script: analysis/global/anchor_site_comparison.py")
    print("=" * 72)

    # ── Lumbini ───────────────────────────────────────────────────────────
    lumbini = find_site(sites, "lumbini")
    if lumbini is None:
        sys.exit("ERROR: Lumbini not found in UNESCO XML")
    dev_g, km_g, tier_g = beru_dev(lumbini["lon"], GERIZIM)
    dev_j, km_j, tier_j = beru_dev(lumbini["lon"], JERUSALEM)
    m_g = km_g * 1000

    print(f"\n  Lumbini ({lumbini['name']})")
    print(f"    lon = {lumbini['lon']}°E")
    print(f"    vs Gerizim:   δ = {dev_g:.6f} beru = {km_g:.2f} km = {m_g:.0f} m  [{tier_g}]")
    print(f"    vs Jerusalem: δ = {dev_j:.6f} beru = {km_j:.2f} km  [{tier_j}]")
    print(f"\n  % LaTeX macros:")
    print(f"  \\newcommand{{\\LumbiniDevKm}}{{{km_g:.2f}}}          % Lumbini deviation from Gerizim (km)")
    print(f"  \\newcommand{{\\LumbiniDevM}}{{{int(round(m_g))}}}           % Lumbini deviation from Gerizim (m)")
    print(f"  \\newcommand{{\\LumbiniDevBeru}}{{{dev_g:.6f}}}      % Lumbini deviation from Gerizim (beru)")
    print(f"  \\newcommand{{\\LumbiniTier}}{{{tier_g}}}            % Lumbini tier relative to Gerizim")
    print(f"  \\newcommand{{\\LumbiniJerDevKm}}{{{km_j:.2f}}}        % Lumbini deviation from Jerusalem (km)")
    print(f"  \\newcommand{{\\LumbiniJerDevBeru}}{{{dev_j:.6f}}}    % Lumbini deviation from Jerusalem (beru)")
    print(f"  \\newcommand{{\\LumbiniJerTier}}{{{tier_j}}}          % Lumbini tier relative to Jerusalem")

    # ── Takht-e Soleyman ─────────────────────────────────────────────────
    takht = find_site(sites, "takht")
    if takht is None:
        sys.exit("ERROR: Takht-e Soleyman not found in UNESCO XML")
    dev_g, km_g, tier_g = beru_dev(takht["lon"], GERIZIM)
    dev_j, km_j, tier_j = beru_dev(takht["lon"], JERUSALEM)
    m_j = km_j * 1000

    print(f"\n  Takht-e Soleyman ({takht['name']})")
    print(f"    lon = {takht['lon']}°E")
    print(f"    vs Gerizim:   δ = {dev_g:.6f} beru = {km_g:.2f} km  [{tier_g}]")
    print(f"    vs Jerusalem: δ = {dev_j:.6f} beru = {km_j:.2f} km = {m_j:.0f} m  [{tier_j}]")
    print(f"\n  % LaTeX macros:")
    print(f"  \\newcommand{{\\TakhtGerizimDevKm}}{{{km_g:.2f}}}       % Takht-e Soleyman deviation from Gerizim (km)")
    print(f"  \\newcommand{{\\TakhtGerizimDevBeru}}{{{dev_g:.6f}}}   % Takht-e Soleyman deviation from Gerizim (beru)")
    print(f"  \\newcommand{{\\TakhtGerizimTier}}{{{tier_g}}}         % Takht-e Soleyman tier relative to Gerizim")
    print(f"  \\newcommand{{\\TakhtJerDevKm}}{{{km_j:.2f}}}          % Takht-e Soleyman deviation from Jerusalem (km)")
    print(f"  \\newcommand{{\\TakhtJerDevM}}{{{int(round(m_j))}}}         % Takht-e Soleyman deviation from Jerusalem (m)")
    print(f"  \\newcommand{{\\TakhtJerDevBeru}}{{{dev_j:.6f}}}      % Takht-e Soleyman deviation from Jerusalem (beru)")
    print(f"  \\newcommand{{\\TakhtJerTier}}{{{tier_j}}}            % Takht-e Soleyman tier relative to Jerusalem")

    # ── Khoja Ahmed Yasawi ───────────────────────────────────────────────
    khoja = find_site(sites, "khoja")
    if khoja is None:
        khoja = find_site(sites, "yasawi")
    if khoja is None:
        sys.exit("ERROR: Khoja Ahmed Yasawi not found in UNESCO XML")
    dev_g, km_g, tier_g = beru_dev(khoja["lon"], GERIZIM)
    dev_j, km_j, tier_j = beru_dev(khoja["lon"], JERUSALEM)

    print(f"\n  Khoja Ahmed Yasawi ({khoja['name']})")
    print(f"    lon = {khoja['lon']}°E")
    print(f"    vs Gerizim:   δ = {dev_g:.6f} beru = {km_g:.2f} km  [{tier_g}]")
    print(f"    vs Jerusalem: δ = {dev_j:.6f} beru = {km_j:.2f} km  [{tier_j}]")
    print(f"\n  % LaTeX macros:")
    print(f"  \\newcommand{{\\KhojaGerizimDevKm}}{{{km_g:.2f}}}       % Khoja Ahmed Yasawi deviation from Gerizim (km)")
    print(f"  \\newcommand{{\\KhojaGerizimDevBeru}}{{{dev_g:.6f}}}   % Khoja Ahmed Yasawi deviation from Gerizim (beru)")
    print(f"  \\newcommand{{\\KhojaGerizimTier}}{{{tier_g}}}         % Khoja Ahmed Yasawi tier relative to Gerizim")
    print(f"  \\newcommand{{\\KhojaJerDevKm}}{{{km_j:.2f}}}          % Khoja Ahmed Yasawi deviation from Jerusalem (km)")
    print(f"  \\newcommand{{\\KhojaJerDevBeru}}{{{dev_j:.6f}}}      % Khoja Ahmed Yasawi deviation from Jerusalem (beru)")
    print(f"  \\newcommand{{\\KhojaJerTier}}{{{tier_j}}}            % Khoja Ahmed Yasawi tier relative to Jerusalem")

    # ── A++ count ─────────────────────────────────────────────────────────
    print(f"\n  A++ sites from Gerizim: {n_app}")
    print(f"  \\newcommand{{\\NappTotal}}{{{n_app}}}             % total A++ sites from Gerizim anchor")
    print(f"  \\newcommand{{\\GerizimAppCount}}{{{n_app}}}        % Gerizim A++ count (alias)")

    # ═════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 72)
    print("  GROUP 22 — World Peace Pagoda (Lumbini)")
    print("  Script: analysis/global/anchor_site_comparison.py")
    print("=" * 72)

    dev_wpp, km_wpp, tier_wpp = beru_dev(WPP_LON, GERIZIM)
    print(f"\n  World Peace Pagoda")
    print(f"    lon = {WPP_LON}°E  (Wikidata / OSM coordinate)")
    print(f"    vs Gerizim: δ = {dev_wpp:.6f} beru = {km_wpp:.2f} km  [{tier_wpp}]")
    print(f"\n  % LaTeX macros:")
    print(f"  \\newcommand{{\\WPPDeltaKm}}{{{km_wpp:.2f}}}          % World Peace Pagoda deviation from Gerizim (km)")
    print(f"  \\newcommand{{\\WPPDeltaBeru}}{{{dev_wpp:.6f}}}      % World Peace Pagoda deviation from Gerizim (beru)")

    print("\n" + "=" * 72)
    print("  DONE — all GROUP 16 + GROUP 22 macros printed above.")
    print("=" * 72)


if __name__ == "__main__":
    main()
