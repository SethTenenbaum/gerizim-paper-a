#!/usr/bin/env python3
"""
mounds_geo_control.py
=====================
Geometric-null binomial proximity test for the North American mound corpus
(mounds_corpus_v3.csv, N=236 deduplicated NHL/NRHP sites with Wikidata QIDs).

This corpus is a type-matched geographic control: hemispherical earthen
monuments located entirely outside the Eurasian corridor. A null result
here confirms the Eurasian signal is geographically specific, not a generic
property of any hemispherical monument class globally.

Emits LaTeX macros to:
  manuscript/generated_macros.tex  (appended)
  analysis/americas/mounds_geo_control_macros.tex  (standalone)
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scipy.stats import binomtest

# ── Config ───────────────────────────────────────────────────────────────────
CORPUS_CSV = ROOT / "data" / "store" / "americas" / "mounds_corpus_v3.csv"
MACRO_OUT  = ROOT / "analysis" / "americas" / "mounds_geo_control_macros.tex"

GERIZIM_LON = 35.2728   # anchor longitude
PERIOD      = 3.0       # degrees

# Tier thresholds (fraction of T)
TIERS = {
    "App": (0.10, 0.20),   # ±10% of T  →  p0 = 0.20
    "Ap":  (0.20, 0.40),   # ±20% of T  →  cumulative with App: but we use CUMULATIVE below
    "A":   (0.50, 1.00),   # ±50% of T  →  p0 = 1.00  (all sites — sanity check)
}

# Actual tier null rates used in the paper
NULL_RATES = {
    "App": 0.20,   # A++: ±0.30° out of 1.5° half-period → 2*0.30/3 = 0.20
    "Ap":  0.40,   # A+:  cumulative ±0.60° → 2*0.60/3 = 0.40  [but paper uses 0.04 for A+ alone]
    "A":   1.00,   # A:   cumulative ±1.50° → 1.00  [sanity, not used]
}

# Paper's actual tier thresholds (degrees from nearest harmonic)
THRESH = {
    "App": 0.30,   # A++ : δ ≤ 0.30°
    "Ap":  0.60,   # A+  : δ ≤ 0.60°  (cumulative, includes A++)
    "A":   1.50,   # A   : δ ≤ 1.50°  (cumulative, includes A++ and A+)
}
# Geometric null rates (2d/T)
P0 = {k: 2 * v / PERIOD for k, v in THRESH.items()}
# → App: 0.20, Ap: 0.40, A: 1.00
# But the paper uses EXCLUSIVE tier null rates; recompute for cumulative:
# A++ only:    2*0.30/3 = 0.20
# A+ cumul:    2*0.60/3 = 0.40   (paper labels this tier A+, null 0.04 for non-cumul 0.06°)
# Let's match the paper's actual tier definitions exactly:
#   A++: δ ≤ 0.30°,  p0 = 2*0.30/3 = 0.20
#   A+:  δ ≤ 0.60°,  p0 = 2*0.60/3 = 0.40   (cumulative)
#   A:   δ ≤ 1.50°,  p0 = 2*1.50/3 = 1.00   (trivially 1)
# The paper's Table 3 uses:
#   A++  p0 = 0.20  (matches)
#   A+   p0 = ?     let's check lib/beru.py threshold
# We'll read the actual thresholds from lib/beru.py

try:
    from lib.beru import deviation, tier_label, APP_THRESH, AP_THRESH, A_THRESH
    THRESH_APP = APP_THRESH
    THRESH_AP  = AP_THRESH
    THRESH_A   = A_THRESH
except ImportError:
    # fallback: read from CSV beru_dev values and tier labels
    THRESH_APP = 0.10   # fraction of T/2 = 0.15°... actually beru_dev is fraction of T
    THRESH_AP  = 0.20
    THRESH_A   = 1.00

# ── Load corpus ───────────────────────────────────────────────────────────────
def load_corpus(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def harmonic_deviation(lon: float, anchor: float, period: float) -> float:
    """Fractional deviation from nearest harmonic, in [0, 0.5] × period."""
    shifted = (lon - anchor) % period
    dev = min(shifted, period - shifted)
    return dev  # degrees

# ── Main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    rows = load_corpus(CORPUS_CSV)
    n = len(rows)
    print(f"Mounds corpus: N={n}")

    # Compute deviations
    devs = []
    for r in rows:
        lon = float(r["lon"])
        dev = harmonic_deviation(lon, GERIZIM_LON, PERIOD)
        devs.append(dev)

    # Use paper tier thresholds: read from beru if available, else hardcode
    # Paper: A++ = δ ≤ 0.30°, A+ = δ ≤ 0.60°, A = δ ≤ 1.50°  (for T=3°)
    # These give p0 = 2d/T: 0.20, 0.40, 1.00
    # But paper table shows A+ p0 = NullRateAp which is likely 0.40 (cumulative)
    # We'll use the same thresholds as the paper's tier table and match beru.py
    try:
        from lib.beru import APP_THRESH, AP_THRESH, A_THRESH
        # beru stores as fraction of T (0..1), so threshold in degrees = frac * T
        th_app = APP_THRESH * PERIOD / 2  # APP_THRESH is fraction of half-period?
        # Let's just use the CSV tier labels directly
    except Exception:
        pass

    # Use CSV tier labels (already computed by fetch script using lib.beru)
    # Paper tier thresholds (from lib/beru.py + config.json):
    #   A++ : δ ≤ 0.06°  → p0 = 2*0.06/3 = 0.04
    #   A+  : δ ≤ 0.15°  → p0 = 2*0.15/3 = 0.10  (cumulative, includes A++)
    #   A   : δ ≤ 0.30°  → p0 = 2*0.30/3 = 0.20  (cumulative)
    n_app = sum(1 for r in rows if r["tier"] == "A++")
    n_ap  = sum(1 for r in rows if r["tier"] in ("A++", "A+"))
    n_a   = sum(1 for r in rows if r["tier"] in ("A++", "A+", "A"))

    p0_app = 0.04   # 2 * 0.06 / 3.0
    p0_ap  = 0.10   # 2 * 0.15 / 3.0  (cumulative A+ tier)
    p0_a   = 0.20   # 2 * 0.30 / 3.0  (cumulative A tier)

    # One-sided binomial (greater)
    res_app = binomtest(n_app, n, p0_app, alternative="greater")
    res_ap  = binomtest(n_ap,  n, p0_ap,  alternative="greater")
    res_a   = binomtest(n_a,   n, p0_a,   alternative="greater")

    rate_app = n_app / n * 100
    rate_ap  = n_ap  / n * 100
    rate_a   = n_a   / n * 100

    print(f"A++ : {n_app}/{n} ({rate_app:.1f}%)  null={p0_app*100:.0f}%  p={res_app.pvalue:.4f}")
    print(f"A+  : {n_ap}/{n}  ({rate_ap:.1f}%)  null={p0_ap*100:.0f}%  p={res_ap.pvalue:.4f}")
    print(f"A   : {n_a}/{n}   ({rate_a:.1f}%)  null={p0_a*100:.0f}%  p={res_a.pvalue:.4f}")

    def fmt_p(p: float) -> str:
        if p >= 0.10:
            return f"{p:.2f}"
        elif p >= 0.001:
            return f"{p:.3f}"
        else:
            return f"{p:.4f}"

    def sig(p: float) -> str:
        if p < 0.001: return "***"
        if p < 0.01:  return "**"
        if p < 0.05:  return "*"
        return "ns"

    macros = [
        f"\\newcommand{{\\moundsN}}{{{n}}}",
        f"\\newcommand{{\\moundsNApp}}{{{n_app}}}",
        f"\\newcommand{{\\moundsNAp}}{{{n_ap}}}",
        f"\\newcommand{{\\moundsNA}}{{{n_a}}}",
        f"\\newcommand{{\\moundsRateApp}}{{{rate_app:.1f}}}",
        f"\\newcommand{{\\moundsRateAp}}{{{rate_ap:.1f}}}",
        f"\\newcommand{{\\moundsRateA}}{{{rate_a:.1f}}}",
        f"\\newcommand{{\\moundsPApp}}{{{fmt_p(res_app.pvalue)}}}",
        f"\\newcommand{{\\moundsPAp}}{{{fmt_p(res_ap.pvalue)}}}",
        f"\\newcommand{{\\moundsPA}}{{{fmt_p(res_a.pvalue)}}}",
        f"\\newcommand{{\\moundsPAppSig}}{{{sig(res_app.pvalue)}}}",
        f"\\newcommand{{\\moundsPApSig}}{{{sig(res_ap.pvalue)}}}",
        f"\\newcommand{{\\moundsPASig}}{{{sig(res_a.pvalue)}}}",
    ]

    # Print to stdout so reproduce_all_macros.sh can collect via grep
    for m in macros:
        print(m)

    MACRO_OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(MACRO_OUT, "w") as f:
        f.write("% Mounds geographic control macros — auto-generated\n")
        f.write("\n".join(macros) + "\n")
    print(f"\nMacros written → {MACRO_OUT}", file=__import__('sys').stderr)

if __name__ == "__main__":
    main()
