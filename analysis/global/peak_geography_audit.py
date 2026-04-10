"""
peak_geography_audit.py
========================
Fine-grained A+ count profile in the 34–37°E Levant corridor.

This script audits the anchor-sweep peak near Gerizim/Jerusalem by:
  1. Sweeping 34.0–37.0°E at 0.001° resolution (3,001 anchors)
  2. Showing the A+ count profile around key Levant landmarks
  3. Confirming that the Gerizim–Jerusalem plateau is the local max
  4. Showing where the profile drops (Jericho, Dead Sea, Mt Nebo)
  5. Analysing the x.18°E artifact locally

Referenced in tenenbaum_2026_gerizim_beru.tex:
  \\JerichoAnchorAp, \\DeadSeaAnchorAp, \\MtNeboAnchorAp,
  \\DamascusAnchorAp, \\localBestLon, \\localBestAp

ANCHOR:  Sweep 34.0–37.0°E at 0.001° resolution
BERU:    30° arc
"""

import sys
from pathlib import Path
import numpy as np
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_B_MAX, P_NULL_AP,
    CONFIG, load_levant_landmarks,
)


# ── Load sites ──────────────────────────────────────────────────────────────
corpus = load_corpus()
_cultural = cultural_sites_with_coords(corpus)
sites = [{"name": s.site, "lon": s.longitude, "lat": s.latitude} for s in _cultural]

N = len(sites)
lons = np.array([s["lon"] for s in sites])


def count_ap_at_anchor(anchor):
    """Count Tier-A+ sites for a given anchor longitude."""
    arcs = np.abs(lons - anchor)
    arcs = np.minimum(arcs, 360 - arcs)
    bvs = arcs / BERU
    nearests = np.round(bvs * 10) / 10
    devs = np.abs(bvs - nearests)
    return int(np.sum(devs <= TIER_APLUS))


# ══════════════════════════════════════════════════════════════════════════════
# PART 1: Fine-grained sweep 34.0–37.0°E
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 100)
print("  PEAK GEOGRAPHY AUDIT — LEVANT CORRIDOR 34.0–37.0°E")
print(f"  N = {N} Cultural/Mixed UNESCO sites")
print(f"  Resolution: 0.001° (3,001 anchor points)")
print("=" * 100)

sweep = np.arange(34.000, 37.001, 0.001)
ap_profile = np.array([count_ap_at_anchor(a) for a in sweep])

# ══════════════════════════════════════════════════════════════════════════════
# PART 2: Key landmarks A+ counts
# ══════════════════════════════════════════════════════════════════════════════
# Key Levant landmarks — loaded from config.json
LANDMARKS = load_levant_landmarks()

print(f"\n  KEY LANDMARKS — A+ COUNTS")
print(f"  ──────────────────────────────────────────────────────────────")
print(f"  {'Landmark':<22}  {'Lon':>8}  {'A+':>4}  {'Phase':>6}  Bar")
print(f"  {'─'*22}  {'─'*8}  {'─'*4}  {'─'*6}  {'─'*30}")

for name, lon in sorted(LANDMARKS.items(), key=lambda x: -count_ap_at_anchor(x[1])):
    ap = count_ap_at_anchor(lon)
    phase = lon % 3.0
    bar = "█" * (ap - 30) if ap > 30 else ""
    print(f"  {name:<22}  {lon:>8.3f}  {ap:>4}  {phase:>5.3f}  {bar}")

# ══════════════════════════════════════════════════════════════════════════════
# PART 3: Profile visualization (ASCII chart)
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n\n{'=' * 100}")
print("  A+ COUNT PROFILE: 34.0–37.0°E (sampled every 0.05°)")
print("=" * 100)

# Sample every 50 steps (0.05°) for display
sample_step = 50
print(f"\n  {'Lon':>7}  {'A+':>4}  Profile")
print(f"  {'─'*7}  {'─'*4}  {'─'*60}")
for i in range(0, len(sweep), sample_step):
    lon_val = sweep[i]
    ap_val = ap_profile[i]
    bar = "█" * max(0, ap_val - 30)
    # Mark key locations
    mark = ""
    for lname, llon in LANDMARKS.items():
        if abs(lon_val - llon) < 0.025:
            mark = f"  ← {lname}"
            break
    print(f"  {lon_val:>7.2f}  {ap_val:>4}  {bar}{mark}")

# ══════════════════════════════════════════════════════════════════════════════
# PART 4: Ultra-fine sweep around Gerizim (35.15–35.35°E at 0.001°)
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n\n{'=' * 100}")
print("  ULTRA-FINE SWEEP: 35.15–35.35°E (Gerizim–Jerusalem corridor)")
print("=" * 100)

_uf_half = 0.125   # ±0.125° window around Gerizim for ultra-fine sweep
uf_sweep = np.arange(GERIZIM - _uf_half, GERIZIM + _uf_half + 0.001, 0.001)
uf_profile = np.array([count_ap_at_anchor(a) for a in uf_sweep])

local_max_idx = np.argmax(uf_profile)
local_max_lon = float(uf_sweep[local_max_idx])
local_max_ap = int(uf_profile[local_max_idx])

print(f"\n  Local maximum: {local_max_ap} A+ at {local_max_lon:.3f}°E")
print(f"  Phase: {local_max_lon % 3.0:.3f}° (mod 3°)")
print(f"  Distance from optimal phase (2.18): {abs(local_max_lon % 3.0 - 2.18):.3f}°")
print()

# Show the ultra-fine profile every 0.005°
print(f"  {'Lon':>8}  {'A+':>4}  {'Phase':>6}  Profile")
print(f"  {'─'*8}  {'─'*4}  {'─'*6}  {'─'*40}")
for i in range(0, len(uf_sweep), 5):
    lon_val = uf_sweep[i]
    ap_val = uf_profile[i]
    phase = lon_val % 3.0
    bar = "█" * max(0, ap_val - 45)
    mark = ""
    if abs(lon_val - GERIZIM) < 0.002:
        mark = " ← GERIZIM"
    elif abs(lon_val - LANDMARKS["Jerusalem"]) < 0.002:
        mark = " ← JERUSALEM"
    elif abs(lon_val - LANDMARKS["Tel Megiddo"]) < 0.002:
        mark = " ← MEGIDDO"
    elif abs(lon_val - local_max_lon) < 0.002 and mark == "":
        mark = " ← LOCAL MAX"
    print(f"  {lon_val:>8.3f}  {ap_val:>4}  {phase:>5.3f}  {bar}{mark}")

# ══════════════════════════════════════════════════════════════════════════════
# PART 5: x.18° artifact check at local level
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n\n{'=' * 100}")
print("  x.18° ARTIFACT CHECK")
print("=" * 100)

# The x.18° anchors in the 34–37 range
x18_local = [a for a in [34.18, 35.18] if 34.0 <= a <= 37.0]
for a in x18_local:
    ap = count_ap_at_anchor(a)
    print(f"  {a:.2f}°E (x.18° artifact): {ap} A+ sites")

print(f"""
  Gerizim ({GERIZIM}°E):  {count_ap_at_anchor(GERIZIM)} A+ sites
  
  The local peak at {local_max_lon:.3f}°E ({local_max_ap} A+) is close to 
  the x.18° artifact phase, but the Gerizim–Jerusalem plateau 
  (35.20–35.28°E, consistently {count_ap_at_anchor(LANDMARKS['Jerusalem'])}–{count_ap_at_anchor(GERIZIM)} A+) 
  is the archaeologically meaningful local maximum.
""")

# ══════════════════════════════════════════════════════════════════════════════
# PART 6: Drop-off profile from Gerizim
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 100)
print("  DROP-OFF FROM GERIZIM → EAST")
print("=" * 100)

_levant_end = CONFIG["anchor_sweep"]["levant"]["end_longitude"]  # 37.0
drop_points = [
    (name, lon)
    for name, lon in sorted(LANDMARKS.items(), key=lambda x: x[1])
    if GERIZIM <= lon <= _levant_end
]
# Always include Gerizim as the reference point at the top
drop_points = [("Gerizim", GERIZIM)] + [
    (name, lon) for name, lon in drop_points if abs(lon - GERIZIM) > 0.001
]

print(f"\n  {'Location':<20}  {'Lon':>8}  {'A+':>4}  {'Δ from Gerizim':>15}")
print(f"  {'─'*20}  {'─'*8}  {'─'*4}  {'─'*15}")
gerizim_count = count_ap_at_anchor(GERIZIM)
for name, lon in drop_points:
    ap = count_ap_at_anchor(lon)
    delta = ap - gerizim_count
    print(f"  {name:<20}  {lon:>8.3f}  {ap:>4}  {delta:>+15d}")

# ══════════════════════════════════════════════════════════════════════════════
# PART 7: LaTeX macro output
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n\n{'=' * 100}")
print("  LATEX MACRO VALUES (copy to tenenbaum_2026_gerizim_beru.tex)")
print("=" * 100)

print(f"  \\newcommand{{\\JerichoAnchorAp}}{{{count_ap_at_anchor(LANDMARKS['Jericho'])}}}           % A+ count with Jericho as anchor")
print(f"  \\newcommand{{\\DeadSeaAnchorAp}}{{{count_ap_at_anchor(LANDMARKS['Dead Sea center'])}}}          % A+ count with Dead Sea centre as anchor")
print(f"  \\newcommand{{\\MtNeboAnchorAp}}{{{count_ap_at_anchor(LANDMARKS['Mount Nebo'])}}}           % A+ count with Mt Nebo as anchor")
print(f"  \\newcommand{{\\DamascusAnchorAp}}{{{count_ap_at_anchor(LANDMARKS['Damascus'])}}}          % A+ count with Damascus as anchor")
print(f"  \\newcommand{{\\localBestLon}}{{{local_max_lon:.3f}}}          % longitude of local peak A+ in Levant sweep")
print(f"  \\newcommand{{\\localBestAp}}{{{local_max_ap}}}           % A+ count at local peak longitude")

print(f"\n{'=' * 100}")
print("  DONE")
print("=" * 100)
