"""
landmarks.py — Named longitudes from config.json.

Provides a single dict of all notable sites used for comparison,
so no analysis script hardcodes coordinate values.
"""

from lib.beru import CONFIG, GERIZIM

# ---------------------------------------------------------------------------
# Landmark dictionary:  name  →  longitude  (degrees E)
# ---------------------------------------------------------------------------
_raw = CONFIG.get("landmarks", {})
LANDMARKS = {
    key: entry["longitude"]
    for key, entry in _raw.items()
    if isinstance(entry, dict) and "longitude" in entry
}

LANDMARK_LABELS = {
    key: entry.get("label", key)
    for key, entry in _raw.items()
    if isinstance(entry, dict) and "label" in entry
}


def get(name: str) -> float:
    """Return the longitude for a named landmark, or raise KeyError."""
    return LANDMARKS[name]


# Convenience constants used by many scripts
JERUSALEM = LANDMARKS.get("jerusalem")
MEGIDDO = LANDMARKS.get("megiddo")
GERIZIM_TEMPLE = GERIZIM  # Now same as primary anchor (35.272°E, UNESCO XML)
MECCA = LANDMARKS.get("mecca")
LUMBINI = LANDMARKS.get("lumbini")


# ---------------------------------------------------------------------------
# Notable-anchor comparison table  (for anchor_uniqueness_audit, etc.)
# ---------------------------------------------------------------------------

NOTABLE_ANCHORS = {
    f"{LANDMARK_LABELS.get(k, k)} ({v:.3f}°E)": v
    for k, v in LANDMARKS.items()
    if k in (
        "jerusalem", "mecca", "giza", "babylon", "persepolis",
        "greenwich", "ujjain", "baghdad", "damascus", "rome",
        "athens", "alexandria",
    )
}
# Gerizim always appears first
NOTABLE_ANCHORS = {
    f"Gerizim ({GERIZIM}°E)": GERIZIM,
    **NOTABLE_ANCHORS,
}


# ---------------------------------------------------------------------------
# Levant landmark dict (for peak_geography_audit, etc.)
# ---------------------------------------------------------------------------

LEVANT_LANDMARKS = {
    LANDMARK_LABELS.get(k, k): v
    for k, v in LANDMARKS.items()
    if k in (
        "jerusalem", "jericho", "dead_sea", "mt_nebo",
        "damascus", "tel_megiddo", "bethlehem", "nablus",
        "mt_ebal", "bethel", "tyre",
    )
}
LEVANT_LANDMARKS["Gerizim"] = GERIZIM
