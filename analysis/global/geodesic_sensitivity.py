"""
geodesic_sensitivity.py
=======================
Computes the latitude-sensitivity check described in §3.1 of the
manuscript: how many A+ sites (angular criterion) are within 6.7 km
physically, and how many would be added/removed by switching to a
geodesic physical-km threshold.

Emits LaTeX macros:
  \\GeodesicApCurrent     — A+ count under current angular criterion  (= NclusterAp)
  \\GeodesicApPhysical    — A+ count under strict 6.7 km geodesic threshold
  \\GeodesicDropOut       — sites removed  (angular A+ but > 6.7 km geodesic)
  \\GeodesicGainIn        — sites gained   (angular non-A+ but ≤ 6.7 km geodesic)
  \\GeodesicThresholdKm   — physical threshold used (6.7)
  \\GeodesicKmPerDegAtSixty — east-west km per degree at 60°N (rounded to 1 dp)

All distances use an east-west haversine approximation:
  d = |Δλ| × 111.32 × cos(lat)
This matches the physical interpretation stated in the manuscript.

USAGE
-----
  python3 analysis/global/geodesic_sensitivity.py

GROUP: Geodesic / latitude sensitivity (Methods §3.1)
"""

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, HARMONIC_STEP,
    deviation, tier_label,
)

# ── Constants ─────────────────────────────────────────────────────────────────
KM_PER_DEG_EQUATORIAL = 111.32          # km per degree of longitude at equator
PHYS_THRESHOLD_KM     = 6.7            # geodesic A+ threshold to compare against
REF_LAT_FOR_DISPLAY   = 60.0           # latitude used in the manuscript example


def nearest_harmonic_lon(site_lon: float,
                         anchor: float = GERIZIM,
                         beru: float = BERU,
                         step: float = HARMONIC_STEP) -> float:
    """Longitude of the nearest 0.1-beru harmonic to site_lon."""
    arc       = site_lon - anchor           # signed
    beru_val  = arc / beru
    nearest   = round(beru_val / step) * step
    return anchor + nearest * beru


def east_west_km(site_lon: float, site_lat: float,
                 anchor: float = GERIZIM) -> float:
    """
    Approximate east-west physical distance (km) from site_lon to
    the nearest 0.1-beru harmonic meridian, evaluated at site_lat.
    """
    h_lon = nearest_harmonic_lon(site_lon)
    return abs(site_lon - h_lon) * KM_PER_DEG_EQUATORIAL * math.cos(math.radians(site_lat))


# ── Load corpus ───────────────────────────────────────────────────────────────
sites = cultural_sites_with_coords(load_corpus())

# ── Classify under both criteria ─────────────────────────────────────────────
angular_ap   = []   # currently A+ (angular)
angular_nonap = []  # currently non-A+

for s in sites:
    dev  = deviation(s.longitude)
    tier = tier_label(dev)
    phys = east_west_km(s.longitude, s.latitude)
    rec  = dict(site=s, tier=tier, dev=dev, phys_km=phys)
    if tier in ("A++", "A+"):
        angular_ap.append(rec)
    else:
        angular_nonap.append(rec)

drop_out = [r for r in angular_ap    if r["phys_km"] > PHYS_THRESHOLD_KM]
gain_in  = [r for r in angular_nonap if r["phys_km"] <= PHYS_THRESHOLD_KM]

ap_current  = len(angular_ap)
ap_physical = ap_current - len(drop_out) + len(gain_in)

km_per_deg_at_sixty = round(
    KM_PER_DEG_EQUATORIAL * math.cos(math.radians(REF_LAT_FOR_DISPLAY)), 1
)

# ── Console report ────────────────────────────────────────────────────────────
SEP = "=" * 72
print(SEP)
print("  GEODESIC SENSITIVITY CHECK — §3.1 latitude-sensitivity note")
print(f"  Physical threshold: {PHYS_THRESHOLD_KM} km (east-west at site latitude)")
print(SEP)
print()
print(f"  Angular A+ (current):         {ap_current}")
print(f"  Sites that DROP OUT (phys > {PHYS_THRESHOLD_KM} km):  {len(drop_out)}")
for r in sorted(drop_out, key=lambda x: x['phys_km'], reverse=True):
    s = r['site']
    print(f"    {s.site[:50]:<50}  lat={s.latitude:6.1f}  phys={r['phys_km']:.2f} km")
print()
print(f"  Sites that GAIN IN (phys ≤ {PHYS_THRESHOLD_KM} km):   {len(gain_in)}")
for r in sorted(gain_in, key=lambda x: x['phys_km']):
    s = r['site']
    ang_km = r['dev'] * BERU * 111.0
    print(f"    {s.site[:50]:<50}  lat={s.latitude:6.1f}  "
          f"ang={ang_km:.2f} km  phys={r['phys_km']:.2f} km")
print()
print(f"  A+ under geodesic threshold:  {ap_physical}  "
      f"({'stronger' if ap_physical > ap_current else 'weaker'} than angular)")
print()
print(f"  km per degree at {REF_LAT_FOR_DISPLAY}°N: {km_per_deg_at_sixty}")
print()

# ── LaTeX macros ──────────────────────────────────────────────────────────────
print("  % LaTeX macros (geodesic sensitivity):")
print(f"  \\newcommand{{\\GeodesicApCurrent}}{{{ap_current}}}"
      f"  % A+ count under angular criterion")
print(f"  \\newcommand{{\\GeodesicApPhysical}}{{{ap_physical}}}"
      f"  % A+ count under {PHYS_THRESHOLD_KM}-km geodesic threshold")
print(f"  \\newcommand{{\\GeodesicDropOut}}{{{len(drop_out)}}}"
      f"  % angular A+ sites that exceed {PHYS_THRESHOLD_KM} km physically")
print(f"  \\newcommand{{\\GeodesicGainIn}}{{{len(gain_in)}}}"
      f"  % non-A+ sites within {PHYS_THRESHOLD_KM} km physically")
print(f"  \\newcommand{{\\GeodesicThresholdKm}}{{{PHYS_THRESHOLD_KM}}}"
      f"  % physical km threshold used in geodesic check")
print(f"  \\newcommand{{\\GeodesicKmPerDegAtSixty}}{{{km_per_deg_at_sixty}}}"
      f"  % east-west km per degree at {REF_LAT_FOR_DISPLAY}°N")
print()
print("  ✓ Geodesic sensitivity macros emitted.")
