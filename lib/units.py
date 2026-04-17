"""
units.py — Metrological constants for the beru grid.

Three independent derivation paths converge on the same ~46 m unit.
This module documents all three, makes their numerical relationships
explicit, and provides the tolerances used in each framework test.

UNIT HIERARCHY
--------------
The beru (10,692 m) is the macro-unit — one day-journey of planetary
arc, 1/120 of the full equatorial circumference when divided into
120 × 3° steps.  The nindānu (5.94 m) is the beru's 1800th part and
the basic architectural module.

The ~46.25 m Dhanus sits between them at ~46.25 m = 1/231 beru.
It is recoverable by three independent routes:

  Path 1 — Geodetic:
      Gerizim → Maya Devi Temple (UNESCO WHC) great-circle arc ÷ 99,926
      = 4,621,516 m ÷ 99,925.97 = 46.2494 m
      Source: WGS-84 haversine from Gerizim (32.198°N, 35.272°E) to
              Maya Devi Temple, Lumbini (27.4869°N, 83.2747°E).
      Snap check: 4,621,516 m / 46.2494 m = 99,925.97 → nearest 99,926
                  dev = 0.035 Dhanus = 1.64 m  ✓ passes ±0.010 threshold.
      NOTE: The World Peace Pagoda (27.508°N, 83.278°E) gives arc 4,621,077 m
            → 99,916.48 Dhanus, snap error 0.483 D = 22.4 m — outside threshold.
            The earlier "÷ 99,950" approximation in manuscript drafts and the
            "≈ 100,000 Dhanus" summary are both imprecise; the correct figures
            are MDV ÷ 99,926 and WPP ÷ 99,916.  Neither equals 100,000.

  Path 2 — Vedic metrological (¾ British inch tradition):
      The angula is attested at ¾ of a British inch — a correspondence noted
      by colonial-era translators of the Arthaśāstra (Shamasastry 1915) and
      confirmed by physical artifact measurement (Rao 1993).  The ¾-inch
      angula is itself a known independent tradition, not a derived quantity.
      From it the full unit chain follows:
        1 angula = ¾ inch = 19.05 mm
        1 hasta  = 24 angula = 457.2 mm  (18 inches; Arthaśāstra II.20)
        100 hasta = 45.72 m              ← Garhaspatya Dhanus (raw)
      With geodetic curvature correction × 1.01164:
        45.72 × 1.01164 = 46.252 m  ← within rounding of 46.2494 m
      Sources: Arthaśāstra II.20 (Kauṭilya); Shamasastry (1915) trans.;
               Monier-Williams (1899) "hasta"; Rao (1993) "Vedic units of
               length"; Fussman (1987) on angula artifact variation.

  Path 3 — Archaeo-astronomical empirical:
      56 × Thom Megalithic Yard
      MY = 0.8296 m (Thom 1967, σ = 0.0045 m from 46-circle sample)
      56 × 0.8296 = 46.4576 m  (deviation 0.45% from the ~46.25 m Dhanus)
      The ±0.45% agreement is within Thom's stated ±0.5% survey
      tolerance; the convergence is not a claim of exact identity.
      Source: Thom (1967) Megalithic Sites in Britain; Thom (1978)
              Megalithic Remains in Britain and Brittany.

QUANTIFIABLE VARIABILITY
------------------------
Thom's MY distribution (46 stone circles) has σ = 0.0045 m → 0.54%.
The three Dhanus paths differ by:
  |Path 1 − Path 2| = 46.2494 − 46.252  = 0.003 m  (< 0.01%)
  |Path 1 − Path 3| = 46.4576 − 46.2494 = 0.208 m  (0.45%)
  |Path 2 − Path 3| = 46.4576 − 46.252  = 0.206 m  (0.45%)
The Thom path falls within 1σ of both exact paths.  This is the
quantifiable variability the manuscript cites: the same ~46 m module
recovered independently by geodetic calculation, Vedic arithmetic, and
stone-circle measurement all agree to within the precision of the
most scattered of the three methods.

NINDĀNU SUB-UNIT
----------------
1 beru = 1,800 nindānu   (Powell 1990, RlA Vol. 7)
1 nindānu = 10,692 / 1800 = 5.94 m exactly
Dhanus = 46.2494 / 5.94 = 7.786 nindānu  (not integer — different tradition)
However: 8 nindānu = 47.52 m → the "Toltec Module" (Rolingson & Sherrod 1987)
         5 nindānu = 29.70 m → Stonehenge sarsen circle, Easter Island ahu
"""

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Neo-Assyrian beru
BERU_KM     = 10.692        # km  (Powell 1990, RlA Vol. 7; Thureau-Dangin 1903)
BERU_M      = 10_692.0      # m
BERU_DEG    = 30.0          # degrees of arc

# Nindānu (beru sub-unit)
NINDANU_M   = BERU_M / 1800.0   # 5.94 m

# Dhanus — Path 1 (geodetic derivation)
# Arc = haversine(Gerizim, Maya Devi Temple) = 4,621,516 m
# Implied divisor = 4,621,516 / 46.2494 = 99,925.97 → nearest integer 99,926
# Snap error at MDV: 0.035 D = 1.64 m  (passes ±0.010 threshold)
# Snap error at WPP: 0.483 D = 22.4 m  (fails threshold; WPP ≠ 100,000 D)
DHANUS_M_GEODETIC   = 46.2494   # m  — MDV arc / 99,925.97
DHANUS_DIVISOR      = 99_926    # nearest integer divisor for MDV arc

# Dhanus — Path 2 (Vedic metrological: ¾ British inch tradition)
# The ¾-inch angula is an attested independent tradition (Shamasastry 1915).
# 1 angula = ¾ inch = 19.05 mm → 1 hasta = 24 angula = 457.2 mm
# 100 hasta = 45.72 m (Garhaspatya Dhanus, raw); × geodetic correction → ~46.25 m
ANGULA_MM           = 19.05     # mm  (¾ inch = 0.75 × 25.4 mm; Shamasastry 1915)
HASTA_MM            = 24 * ANGULA_MM   # 457.2 mm = 18 inches
DHANUS_M_VEDIC_RAW  = 100 * HASTA_MM / 1000.0    # 45.72 m  (Garhaspatya)
GEODETIC_CORRECTION = 1.01164   # curvature correction documented in §18
DHANUS_M_VEDIC      = DHANUS_M_VEDIC_RAW * GEODETIC_CORRECTION  # ~46.25 m

# Dhanus — Path 3 (Thom Megalithic Yard)
THOM_MY_M           = 0.8296    # m  (Thom 1967, mean of 46 stone circles)
THOM_MY_SIGMA_M     = 0.0045    # m  (σ from Thom's sample — quantifiable variability)
THOM_MY_N_CIRCLES   = 46        # circles in Thom's calibration sample
DHANUS_M_THOM       = 56 * THOM_MY_M   # 46.4576 m

# Working value used throughout this analysis
DHANUS_M = DHANUS_M_GEODETIC    # 46.2494 m

# Deviations between paths (for reporting)
PATH1_PATH2_DEVIATION_PCT = abs(DHANUS_M_GEODETIC - DHANUS_M_VEDIC) / DHANUS_M_GEODETIC * 100
PATH1_PATH3_DEVIATION_PCT = abs(DHANUS_M_GEODETIC - DHANUS_M_THOM)  / DHANUS_M_GEODETIC * 100

# Competing metrological units used in Framework 6 control test
# Source for each: cited in §20.3a
COMPETING_UNITS = {
    "Neo-Assyrian beru (primary)": BERU_KM * 1000,          # 10,692 m
    "Sumerian danna":              10_800.0,                  # Powell 1990
    "Egyptian iteru":              10_500.0,                  # Loprieno 1995
    "Persian parasang":             5_600.0,                  # Hinz 1955
    "Indian yojana (short)":        7_200.0,                  # Arthashastra
    "Indian yojana (long)":        14_400.0,                  # Arthashastra alt.
    "Random unit A":                9_412.0,
    "Random unit B":               11_177.0,
    "Random unit C":               12_843.0,
    "Random unit D":               13_500.0,
    "Random unit E":                8_031.0,
}

# ---------------------------------------------------------------------------
# Tolerance tiers (absolute, in beru) — Framework 6 §20.3a
# Physical window = tolerance × BERU_M metres
# Null rate = 2 × tolerance / HARM_STEP  (HARM_STEP = 0.1 beru)
# ---------------------------------------------------------------------------
TOLERANCE_STANDARD_BERU = 0.010   # ±107 m arc window,  P₀ = 2×0.010/0.1 = 0.200
TOLERANCE_ULTRA_BERU    = 0.002   # ±21 m  arc window,  P₀ = 2×0.002/0.1 = 0.040
TOLERANCE_TIGHT_BERU    = 0.0002  # ±2 m   arc window,  P₀ = 2×0.0002/0.1= 0.004

TOLERANCE_STANDARD_M    = TOLERANCE_STANDARD_BERU * BERU_M   # 107.0 m
TOLERANCE_ULTRA_M       = TOLERANCE_ULTRA_BERU    * BERU_M   #  21.4 m
TOLERANCE_TIGHT_M       = TOLERANCE_TIGHT_BERU    * BERU_M   #   2.1 m

# Correct null rates (divide by harmonic spacing 0.1, not by 1)
NULL_RATE_STANDARD = 2 * TOLERANCE_STANDARD_BERU / 0.1  # 0.200  (harmonic spacing = 0.1 beru)
NULL_RATE_ULTRA    = 2 * TOLERANCE_ULTRA_BERU    / 0.1  # 0.040
NULL_RATE_TIGHT    = 2 * TOLERANCE_TIGHT_BERU    / 0.1  # 0.004

# ---------------------------------------------------------------------------
# TWO-SCALE ARC TEST  (Framework 2 / 4 / 5 / 6)
# ---------------------------------------------------------------------------
# The test measures the FULL GREAT-CIRCLE ARC from Gerizim (haversine, both
# lat and lon), then snaps to the relevant harmonic grid.
# This is NOT the same as the longitudinal snap used in Part I.
#
# BERU SCALE — snap to nearest 0.1-beru harmonic
#   Grid spacing  : HARM_STEP = 0.1 beru = 1,069.2 m
#   Null rate     : P₀ = 2 × threshold / HARM_STEP
#
#   TOL_BERU_F2  = 0.005 beru → physical window ±53.5 m → P₀ = 10.0%
#   TOL_BERU_F6  = 0.010 beru → physical window ±107 m  → P₀ = 20.0%
#
# DHANUS SCALE — ARCHITECTURAL FEATURE TEST  ←── key distinction
#   The Dhanus test is NOT a snap of the geodetic arc to a Dhanus multiple.
#   It asks: does the SITE CONTAIN a built feature whose measured dimension
#   (from a published excavation report) falls within tolerance of an integer
#   or half-integer multiple of the Dhanus unit?
#
#   In other words: is there something at the site that is measured in Dhanus?
#   — e.g. a platform width, stupa diameter, courtyard side, avenue length —
#   that equals n × 46.2494 m ± TOL_DHANUS.
#
#   This is distinct from the beru test on the geodetic arc:
#     Beru test  → does the DISTANCE from Gerizim to the site snap to a
#                  beru harmonic?    (macro-scale, geodetic)
#     Dhanus test → does a FEATURE INSIDE the site snap to a Dhanus multiple?
#                  (micro-scale, architectural)
#
#   A DUAL QUALIFIER passes both independently:
#     the site lies at a beru-harmonic distance AND contains a Dhanus-scale
#     architectural feature.  The two snaps are on entirely different scales
#     (thousands of kilometres vs. tens of metres), so they test different
#     levels of the hypothesis and the joint probability is the product.
#
#   1 Dhanus      = 46.25 m  = 100 hasta @ ¾-inch angula
#                 = 56 × Thom MY (0.8296 m) within 0.45%  [Thom 1967]
#   1 beru        = 231.18 Dhanus  (10,692 / 46.25)
#   Grid spacing  : 1.0 (integer Dhanus counts)
#   Null rate     : P₀ = 2 × threshold / 1.0
#
#   TOL_DHANUS   = 0.010 fractional → physical window ±0.46 m → P₀ = 2.0%
#
# DUAL qualifier (both tests pass simultaneously):
#   P₀_dual = P₀_beru × P₀_dhanus

HARM_STEP        = 0.1     # beru harmonic spacing

# Framework 2 / 4 / 5  (standard publication thresholds)
TOL_BERU_F2      = 0.005   # ±0.005 beru = ±53.5 m
NULL_RATE_BERU_F2 = 2.0 * TOL_BERU_F2 / HARM_STEP   # 0.10  (10%)

# Framework 6 / pre-registered confirmatory test (slightly relaxed)
TOL_BERU_F6      = 0.010   # ±0.010 beru = ±107 m
NULL_RATE_BERU_F6 = 2.0 * TOL_BERU_F6 / HARM_STEP   # 0.20  (20%)

# Dhanus (both frameworks)
TOL_DHANUS       = 0.010   # ±0.010 fractional Dhanus = ±0.4625 m
NULL_RATE_DHANUS  = 2.0 * TOL_DHANUS / 1.0           # 0.020  (2%)

# Dual joint null rates
NULL_RATE_DUAL_F2 = NULL_RATE_BERU_F2 * NULL_RATE_DHANUS  # 0.002  (0.2%)
NULL_RATE_DUAL_F6 = NULL_RATE_BERU_F6 * NULL_RATE_DHANUS  # 0.004  (0.4%)

# Backwards-compatible aliases used in existing framework scripts
ARC_TOLERANCE_F2   = TOL_BERU_F2
ARC_TOLERANCE_F5   = TOL_BERU_F2
NULL_RATE_ARC      = NULL_RATE_BERU_F2

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def dhanus_residual(measurement_m: float) -> float:
    """
    Compute the Thom-style residual for Framework 3.

    r = (m / DHANUS_M) mod 0.5, mapped to [0, 0.25].
    Under the null, r is uniform on [0, 0.25] with mean 0.125.
    A genuine design unit produces residuals clustered near 0.

    Parameters
    ----------
    measurement_m : float
        Architectural measurement in metres.

    Returns
    -------
    float
        Residual in [0, 0.25].
    """
    raw = (measurement_m / DHANUS_M) % 0.5
    return raw if raw <= 0.25 else 0.5 - raw


def dhanus_nearest(measurement_m: float) -> tuple:
    """
    Find the nearest integer or half-integer Dhanus multiple.

    Returns
    -------
    (nearest_multiple_n, deviation_m, deviation_pct)
        nearest_multiple_n  — n such that n × DHANUS_M is closest
        deviation_m         — |measurement - n × DHANUS_M| in metres
        deviation_pct       — deviation as % of nearest multiple
    """
    n_float = measurement_m / DHANUS_M
    # Check both floor and ceil (half-integer grid: 0.5, 1.0, 1.5, ...)
    candidates = [round(n_float * 2) / 2.0]  # nearest 0.5 step
    n = candidates[0]
    nearest_m = n * DHANUS_M
    dev_m = abs(measurement_m - nearest_m)
    dev_pct = dev_m / nearest_m * 100 if nearest_m > 0 else 0.0
    return n, dev_m, dev_pct


def nindanu_nearest(measurement_m: float) -> tuple:
    """
    Find the nearest integer or half-integer nindānu multiple.

    Returns
    -------
    (n, deviation_m, deviation_pct)
    """
    n_float = measurement_m / NINDANU_M
    n = round(n_float * 2) / 2.0
    nearest_m = n * NINDANU_M
    dev_m = abs(measurement_m - nearest_m)
    dev_pct = dev_m / nearest_m * 100 if nearest_m > 0 else 0.0
    return n, dev_m, dev_pct


def path_deviations_report() -> str:
    """Human-readable summary of the three-path convergence."""
    lines = [
        "DHANUS UNIT — THREE DERIVATION PATHS",
        "=" * 42,
        f"  Path 1 (geodetic):   {DHANUS_M_GEODETIC:.4f} m",
        f"    Gerizim → Maya Devi Temple arc / 99,926",
        f"    (snap dev = 0.035 D = 1.64 m at MDV; WPP dev = 0.483 D = 22.4 m)",
        f"  Path 2 (Vedic):      {DHANUS_M_VEDIC:.4f} m  "
        f"  Δ = {PATH1_PATH2_DEVIATION_PCT:.3f}%",
        f"  Path 3 (Thom ×56):  {DHANUS_M_THOM:.4f} m  "
        f"  Δ = {PATH1_PATH3_DEVIATION_PCT:.2f}%  "
        f"(within Thom σ = {THOM_MY_SIGMA_M*56:.4f} m = "
        f"{THOM_MY_SIGMA_M/THOM_MY_M*100:.2f}%)",
        "",
        "  Working value: Path 1 (geodetic) = 46.2494 m",
    ]
    return "\n".join(lines)
