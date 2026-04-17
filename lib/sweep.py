"""
sweep.py — Anchor sweep analysis.

Tests whether Mount Gerizim is a uniquely strong anchor by repeating
the A+ count at thousands of trial anchor longitudes and computing
Gerizim's percentile rank.
"""

import numpy as np
from typing import List, Dict, Optional

from lib.beru import CONFIG, BERU, TIER_APLUS, TIER_A_MAX, GERIZIM, deviation


def count_hits(anchor: float, longitudes: List[float],
               threshold: float, wrap: bool = False) -> int:
    """
    Count how many longitudes fall within `threshold` beru of a
    0.1-beru harmonic, measured from a given anchor.

    Parameters
    ----------
    anchor : float
        Anchor longitude in degrees.
    longitudes : list of float
        Site longitudes.
    threshold : float
        Maximum deviation in beru (e.g. TIER_APLUS = 0.002).
    wrap : bool
        If True, use great-circle wrapping (for global sweeps).
    """
    count = 0
    for lon in longitudes:
        dev = deviation(lon, anchor=anchor, wrap=wrap)
        if dev <= threshold:
            count += 1
    return count


def run_sweep(longitudes: List[float],
              start: Optional[float] = None,
              end: Optional[float] = None,
              step: Optional[float] = None,
              sweep_name: str = "levant") -> Dict:
    """
    Sweep a range of anchor longitudes and count A+ and A hits at each.

    Uses sweep parameters from config.json by default. Set sweep_name
    to "global" for the 0–360° sweep, or "levant" for the 34–37° sweep.

    Parameters
    ----------
    longitudes : list of float
        Site longitudes.
    start, end, step : float, optional
        Override sweep range and resolution.
    sweep_name : str
        Key in config["anchor_sweep"] to use for defaults.

    Returns
    -------
    dict with sweep_anchors, counts_aplus, counts_a arrays.
    """
    cfg = CONFIG["anchor_sweep"].get(sweep_name, CONFIG["anchor_sweep"]["levant"])
    start = start if start is not None else cfg["start_longitude"]
    end = end if end is not None else cfg["end_longitude"]
    step = step if step is not None else cfg["step"]
    wrap = (end - start) > 180  # auto-detect global sweep

    anchors = np.arange(start, end + step / 2, step)
    counts_aplus = np.array([
        count_hits(a, longitudes, TIER_APLUS, wrap=wrap) for a in anchors
    ])
    counts_a = np.array([
        count_hits(a, longitudes, TIER_A_MAX, wrap=wrap) for a in anchors
    ])

    return {
        "sweep_anchors": anchors,
        "counts_aplus": counts_aplus,
        "counts_a": counts_a,
    }


def percentile_rank(sweep_counts: np.ndarray, value: int) -> float:
    """
    What fraction of sweep anchors have a count ≤ this value?

    A high percentile means this anchor is unusually strong.

    Parameters
    ----------
    sweep_counts : ndarray
        Array of hit counts from a sweep.
    value : int
        Hit count to rank.

    Returns
    -------
    float
        Percentile (0–100).
    """
    return float(np.mean(sweep_counts <= value)) * 100


def summarize_anchor(longitudes: List[float], anchor: float,
                     sweep_result: Dict, wrap: bool = False) -> Dict:
    """
    Compute hit counts and percentile ranks for one anchor.

    Parameters
    ----------
    longitudes : list of float
        Site longitudes.
    anchor : float
        Anchor longitude to evaluate.
    sweep_result : dict
        Output of run_sweep().
    wrap : bool
        Use great-circle wrapping.

    Returns
    -------
    dict with count_aplus, count_a, pctile_aplus, pctile_a.
    """
    n_aplus = count_hits(anchor, longitudes, TIER_APLUS, wrap=wrap)
    n_a = count_hits(anchor, longitudes, TIER_A_MAX, wrap=wrap)
    pctile_aplus = percentile_rank(sweep_result["counts_aplus"], n_aplus)
    pctile_a = percentile_rank(sweep_result["counts_a"], n_a)
    return {
        "count_aplus": n_aplus,
        "count_a": n_a,
        "pctile_aplus": pctile_aplus,
        "pctile_a": pctile_a,
    }
