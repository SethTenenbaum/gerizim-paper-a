"""
beru.py — Beru-unit calculations.

All beru math lives here. Constants are loaded from config.json
so that no script hardcodes values independently.

A "beru" is 30° of arc, attested in the Babylonian astronomical compendium
MUL.APIN (c. 1000 BCE). The harmonic grid divides the beru into 10 equal
steps of 0.1 beru = 3° each. A site's "deviation" is its distance from
the nearest such harmonic line.
"""

import json
from pathlib import Path

# ---------------------------------------------------------------------------
# Load constants from the shared config file.
# ---------------------------------------------------------------------------
_CONFIG_PATH = Path(__file__).parent.parent / "config.json"
with open(_CONFIG_PATH) as f:
    CONFIG = json.load(f)

# Anchor longitude
GERIZIM = CONFIG["anchors"]["gerizim"]["longitude"]

# Unit parameters
BERU = CONFIG["units"]["beru"]["degrees"]            # 30.0
HARMONIC_STEP = CONFIG["units"]["harmonic_step"]      # 0.1
KM_PER_DEGREE = CONFIG["units"]["km_per_degree"]      # 111.0

# Tier thresholds (beru deviation)
TIER_APP = CONFIG["tiers"]["A++"]["max_deviation_beru"]   # 0.0002
TIER_APLUS = CONFIG["tiers"]["A+"]["max_deviation_beru"]  # 0.002
TIER_A_MAX = CONFIG["tiers"]["A"]["max_deviation_beru"]   # 0.010
TIER_B_MAX = CONFIG["tiers"]["B"]["max_deviation_beru"]   # 0.050

# Geometric null rates
P_NULL_APP = CONFIG["null_rates"]["tier_app"]    # 0.004
P_NULL_AP  = CONFIG["null_rates"]["tier_aplus"]  # 0.04
P_NULL_A   = CONFIG["null_rates"]["tier_a"]      # 0.20


# ---------------------------------------------------------------------------
# Core calculations
# ---------------------------------------------------------------------------

def deviation(lon: float, anchor: float = GERIZIM, beru: float = BERU,
              wrap: bool = False) -> float:
    """
    How far is this longitude from the nearest 0.1-beru harmonic?

    Takes a site longitude and an anchor longitude. Converts the arc
    between them into beru units, finds the nearest 0.1-beru line,
    and returns the absolute distance to that line (in beru).

    The result is always in the range [0, 0.05].

    Parameters
    ----------
    lon : float
        Site longitude in degrees.
    anchor : float
        Anchor longitude in degrees (default: Mount Gerizim).
    beru : float
        Size of one beru in degrees (default: 30.0).
    wrap : bool
        If True, take the shorter arc on a 360° circle.
        Use this for global anchor sweeps where arcs can exceed 180°.

    Example
    -------
    >>> deviation(65.272, anchor=35.272)  # exactly 1.0 beru
    0.0
    """
    arc = abs(lon - anchor)
    if wrap:
        arc = min(arc, 360.0 - arc)
    beru_val = arc / beru
    nearest_harmonic = round(beru_val / HARMONIC_STEP) * HARMONIC_STEP
    return abs(beru_val - nearest_harmonic)


def full_calculation(lon: float, anchor: float = GERIZIM, beru: float = BERU,
                     wrap: bool = False) -> dict:
    """
    Full beru breakdown for a single site longitude.

    Returns a dict with:
        arc_deg:  angular distance from anchor (degrees)
        beru_val: distance in beru units
        nearest:  nearest 0.1-beru harmonic value
        dev:      deviation from nearest harmonic (beru)
        dev_km:   deviation in approximate km
        tier:     tier label (A++, A+, A, B, or C)
    """
    arc = abs(lon - anchor)
    if wrap:
        arc = min(arc, 360.0 - arc)
    beru_val = arc / beru
    nearest = round(beru_val / HARMONIC_STEP) * HARMONIC_STEP
    dev = abs(beru_val - nearest)
    return {
        "arc_deg": arc,
        "beru_val": beru_val,
        "nearest": nearest,
        "dev": dev,
        "dev_km": dev * beru * KM_PER_DEGREE,
        "tier": tier_label(dev),
    }


def deviation_at_spacing(lon: float, spacing: float,
                         anchor: float = GERIZIM, beru: float = BERU) -> float:
    """
    Deviation from the nearest harmonic at an arbitrary spacing.

    Used by the unit-sweep analysis (Tables 6–7) to test whether the
    0.1-beru spacing is uniquely significant.

    Parameters
    ----------
    lon : float
        Site longitude in degrees.
    spacing : float
        Harmonic spacing in beru (e.g. 0.10, 0.05, 0.20).
    anchor : float
        Anchor longitude in degrees.
    beru : float
        Size of one beru in degrees.

    Returns
    -------
    float
        Deviation in beru from the nearest harmonic at the given spacing.
    """
    arc = abs(lon - anchor)
    beru_val = arc / beru
    nearest = round(beru_val / spacing) * spacing
    return abs(beru_val - nearest)


# ---------------------------------------------------------------------------
# Tier classification
# ---------------------------------------------------------------------------

def tier_label(dev: float) -> str:
    """
    Classify a beru deviation into a tier.

    Tiers (from closest to farthest from a harmonic):
        A++  deviation ≤ 0.0002 beru  (about 0.67 km)
        A+   deviation ≤ 0.002  beru  (about 6.7 km)
        A    deviation ≤ 0.010  beru  (about 33 km)
        B    deviation ≤ 0.050  beru  (about 167 km)
        C    deviation > 0.050  beru
    """
    if dev <= TIER_APP:
        return "A++"
    if dev <= TIER_APLUS:
        return "A+"
    if dev <= TIER_A_MAX:
        return "A"
    if dev <= TIER_B_MAX:
        return "B"
    return "C"


def is_aplus(tier: str) -> bool:
    """Is the site in the A+ or A++ tier?"""
    return tier in ("A++", "A+")


def is_a_or_better(tier: str) -> bool:
    """Is the site in tier A, A+, or A++?"""
    return tier in ("A++", "A+", "A")


def dev_to_km(dev: float) -> float:
    """Convert a beru deviation to approximate kilometres."""
    return dev * BERU * KM_PER_DEGREE


# ---------------------------------------------------------------------------
# Keyword access
# ---------------------------------------------------------------------------

def load_keywords(key: str) -> list:
    """
    Return a keyword list from the config ``keywords`` section.

    Handles both old-style flat lists and new-style dicts with
    "unambiguous"/"ambiguous" sub-lists. In the dict case, returns
    the combined list (unambiguous + ambiguous) for backward
    compatibility with scripts that do simple substring matching.

    For context-aware classification, use ``lib.founding_filter``
    instead.

    Parameters
    ----------
    key : str
        One of the named keyword groups in config.json:
        ``"founding_capital"``, ``"sacred_origin"``,
        ``"founding_monument"``, ``"founding_axis"``,
        ``"ancient_landscape"``.

    Returns
    -------
    list of str
        All keyword strings (unambiguous + ambiguous) for the group.

    Raises
    ------
    KeyError
        If *key* is not present in the ``keywords`` section.
    """
    val = CONFIG["keywords"][key]
    if isinstance(val, list):
        return val
    # New-style dict with unambiguous/ambiguous sub-lists
    if isinstance(val, dict):
        return val.get("unambiguous", []) + val.get("ambiguous", [])
    return val


def load_religion_sets() -> list:
    """
    Return the religion keyword sets as a list of (name, [kw, ...]) tuples,
    in the order they appear in config.json.
    """
    sets = CONFIG["keywords"]["religion_sets"]
    return list(sets.items())


def load_notable_anchors() -> dict:
    """
    Return the notable historical/archaeological anchors from config.json.

    Returns
    -------
    dict
        {label: longitude_float}  — ordered as in config.json.
        The ``"note"`` key is stripped automatically.
    """
    raw = CONFIG["notable_anchors"].copy()
    raw.pop("note", None)
    return raw


def load_levant_landmarks() -> dict:
    """
    Return the Levant landmark longitudes from config.json.

    Returns
    -------
    dict
        {name: longitude_float}
    """
    raw = CONFIG["levant_landmarks"].copy()
    raw.pop("note", None)
    return raw
