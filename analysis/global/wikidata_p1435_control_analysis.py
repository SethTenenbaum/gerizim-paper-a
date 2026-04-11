"""
wikidata_p1435_control_analysis.py
==================================
Computes all Wikidata P1435 global control macros (GROUP 23) for
paper_a_primary_unesco.tex.

Compares the UNESCO Cultural/Mixed corpus (N=1011) against the
Wikidata P1435 globally-balanced heritage control corpus (N~51,102)
for beru-harmonic enrichment.

TESTS PERFORMED
---------------
  (a) Binomial test: P1435 control vs 4% geometric null
  (b) 2-proportion z-test: UNESCO > P1435 control
  (c) Odds ratio with 95% CI (Woolf log-odds SE)
  (d) Permutation test (50,000 draws, seed=42)
  (e) Mantel-Haenszel stratified by UNESCO region
  (f) Leave-one-out: drop Europe, recompute OR

DATA
----
  UNESCO:  data/store/unesco/unesco.xml (loaded via data.unesco_corpus)
  Control: data/store/wikidata/p1435_global_control.csv
           (produced by data/scripts/fetch_p1435_global.py)

OUTPUT
------
  Prints LaTeX macro definitions to stdout.
  All macros map 1-to-1 to \\WDctrl* commands in the manuscript.

USAGE
-----
  python3 analysis/global/wikidata_p1435_control_analysis.py

CHANGELOG
---------
  1.1.0 (2026-04-11)
    - Added Arab States region breakdown macros (GROUP 23g):
      \\WDctrlArabN, \\WDctrlArabAp, \\WDctrlArabRate.
      These are descriptive only (P1435 control corpus, not matched
      against UNESCO Arab States); used in manuscript to ground the
      Arabian Peninsula observation in code-generated numbers.
  1.0.4 (2026-04-10)
    - BUG FIX: \\WDctrlBinomP was emitted as "0.00" when the binomial
      p-value rounds below 0.01 under :.2f formatting.  Replaced with
      fmt_p() helper that outputs "$< 0.001$" for p < 0.001, consistent
      with the corrected macro in paper_a_primary_unesco.tex and
      generated_macros.tex.
    - Console print (a) changed from :.2f to :.4g for the same reason.
  1.0.0 (2026-04-09)
    - Initial release.
"""

from __future__ import annotations

__version__ = "1.1.0"

import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.stats import binomtest

# ── Repo root ─────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import json as _json
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, TIER_APLUS, P_NULL_AP

# ── Constants from config.json ────────────────────────────────────────────────
_CFG           = _json.loads((ROOT / "config.json").read_text())
BERU_DEG       = BERU                          # 30.0° — from lib.beru
HARMONIC_STEP  = BERU * 0.1                   # 0.1 beru = 3.0°
HALF_STEP      = HARMONIC_STEP / 2
P_NULL         = P_NULL_AP                    # 0.04 — geometric null for Tier-A+
CTRL_CSV       = ROOT / "data" / "store" / "wikidata" / "p1435_global_control.csv"
N_PERM         = _CFG["simulation"]["n_permutations_wikidata"]  # 50,000
SEED           = _CFG["simulation"]["random_seed"]              # 42

# Region boxes matching fetch_p1435_global.py
REGION_BOXES = {
    "Europe":       {"lat": (35, 72),  "lon": (-25, 45)},
    "Asia-Pacific": {"lat": (-50, 55), "lon": (45, 180)},
    "Americas":     {"lat": (-60, 72), "lon": (-180, -30)},
    "Africa":       {"lat": (-40, 37), "lon": (-20, 55)},
    "Arab States":  {"lat": (12, 42),  "lon": (25, 60)},
}


def assign_region(lat: float, lon: float) -> str:
    for reg, box in REGION_BOXES.items():
        if box["lat"][0] <= lat <= box["lat"][1] and box["lon"][0] <= lon <= box["lon"][1]:
            return reg
    return "Other"


def beru_dev(lon: float) -> float:
    arc  = abs(lon - GERIZIM)
    arc  = min(arc, 360 - arc)
    bv   = arc / BERU
    near = round(bv * 10) / 10
    return abs(bv - near)


def load_control():
    """Load P1435 CSV file and compute beru deviations."""
    rows = []
    with open(CTRL_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["lat"])
                lon = float(row["lon"])
            except (ValueError, KeyError):
                continue
            dev = beru_dev(lon)
            region = row.get("region", assign_region(lat, lon))
            rows.append({
                "lat": lat, "lon": lon, "dev": dev,
                "aplus": dev <= TIER_APLUS, "region": region
            })
    return rows


def load_unesco():
    """Load UNESCO Cultural/Mixed sites and compute beru deviations."""
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    rows = []
    for s in cultural:
        dev = beru_dev(s.longitude)
        region = assign_region(s.latitude, s.longitude)
        rows.append({
            "lat": s.latitude, "lon": s.longitude, "dev": dev,
            "aplus": dev <= TIER_APLUS, "region": region
        })
    return rows


def mantel_haenszel(strata: list[tuple[int, int, int, int]]):
    """Compute Mantel-Haenszel common OR and CMH chi-squared."""
    num_or = 0.0
    den_or = 0.0
    num_chi = 0.0
    den_chi = 0.0
    for a, n1, c, n2 in strata:
        b = n1 - a
        d = n2 - c
        N = n1 + n2
        if N == 0:
            continue
        num_or += a * d / N
        den_or += b * c / N
        m = a + c
        E_a = n1 * m / N
        V_a = n1 * n2 * m * (N - m) / (N * N * (N - 1)) if N > 1 else 0
        num_chi += a - E_a
        den_chi += V_a
    common_or = num_or / den_or if den_or > 0 else float("inf")
    chi2 = num_chi ** 2 / den_chi if den_chi > 0 else 0
    p = 1 - stats.chi2.cdf(chi2, 1)
    return common_or, chi2, p


def main():
    if not CTRL_CSV.exists():
        sys.exit(f"\n  ERROR: {CTRL_CSV} not found.\n"
                 f"  Run: python3 data/scripts/fetch_p1435_global.py\n")

    ctrl = load_control()
    unesco = load_unesco()

    n_c = len(ctrl)
    k_c = sum(1 for r in ctrl if r["aplus"])
    r_c = k_c / n_c

    n_u = len(unesco)
    k_u = sum(1 for r in unesco if r["aplus"])
    r_u = k_u / n_u

    print("=" * 72)
    print("  GROUP 23 — Wikidata P1435 Global Heritage Control")
    print(f"  Script: analysis/global/wikidata_p1435_control_analysis.py  v{__version__}")
    print("=" * 72)

    print(f"\n  UNESCO:  N = {n_u:,}, A+ = {k_u}, rate = {r_u*100:.2f}%")
    print(f"  Control: N = {n_c:,}, A+ = {k_c}, rate = {r_c*100:.2f}%")

    # (a) Binomial: control vs 4% null
    p_binom_c = binomtest(k_c, n_c, P_NULL, alternative="two-sided").pvalue

    # (b) 2-proportion z-test
    p_hat = (k_u + k_c) / (n_u + n_c)
    se = math.sqrt(p_hat * (1 - p_hat) * (1/n_u + 1/n_c))
    z = (r_u - r_c) / se if se > 0 else 0
    p_z = 1 - stats.norm.cdf(z)

    # (c) Odds ratio with 95% CI
    a, b = k_u, n_u - k_u
    c, d = k_c, n_c - k_c
    or_global = (a * d) / (b * c) if b * c > 0 else float("inf")
    log_or = math.log(or_global)
    se_log_or = math.sqrt(1/a + 1/b + 1/c + 1/d)
    or_lo = math.exp(log_or - 1.96 * se_log_or)
    or_hi = math.exp(log_or + 1.96 * se_log_or)

    # (d) Permutation test — vectorised via NumPy (fast)
    rng_np    = np.random.default_rng(SEED)
    n_total   = n_u + n_c
    n_hits    = k_u + k_c
    obs_diff  = r_u - r_c
    # Each permutation: draw k positions (without replacement) for the
    # "UNESCO" slot and count hits; the rest go to "control".
    # Equivalent to shuffling a 0/1 array but ~100× faster than Python list.
    perm_u    = rng_np.hypergeometric(n_hits, n_total - n_hits, n_u, size=N_PERM)
    perm_c    = n_hits - perm_u
    r_perm_u  = perm_u / n_u
    r_perm_c  = perm_c / n_c
    p_perm    = float(np.mean(r_perm_u - r_perm_c >= obs_diff))

    # (e) Mantel-Haenszel by region
    regions = sorted(set(r["region"] for r in unesco) | set(r["region"] for r in ctrl))
    strata = []
    for reg in regions:
        if "Americas" in reg:
            continue  # Americas excluded from MH (different population)
        ku_r = sum(1 for r in unesco if r["region"] == reg and r["aplus"])
        nu_r = sum(1 for r in unesco if r["region"] == reg)
        kc_r = sum(1 for r in ctrl if r["region"] == reg and r["aplus"])
        nc_r = sum(1 for r in ctrl if r["region"] == reg)
        if nu_r > 0 and nc_r > 0:
            strata.append((ku_r, nu_r, kc_r, nc_r))

    mh_or, mh_chi2, mh_p = mantel_haenszel(strata)

    # (f) Leave-one-out: drop Europe
    ku_noeu = sum(1 for r in unesco if r["region"] != "Europe" and r["aplus"])
    nu_noeu = sum(1 for r in unesco if r["region"] != "Europe")
    kc_noeu = sum(1 for r in ctrl if r["region"] != "Europe" and r["aplus"])
    nc_noeu = sum(1 for r in ctrl if r["region"] != "Europe")
    if nc_noeu > 0 and nu_noeu > 0:
        or_noeu, _ = stats.fisher_exact(
            [[ku_noeu, nu_noeu - ku_noeu], [kc_noeu, nc_noeu - kc_noeu]],
            alternative="greater"
        )
    else:
        or_noeu = float("inf")

    # (g) Arab States region breakdown (P1435 control only — descriptive)
    nc_arab  = sum(1 for r in ctrl if r["region"] == "Arab States")
    kc_arab  = sum(1 for r in ctrl if r["region"] == "Arab States" and r["aplus"])
    rc_arab  = kc_arab / nc_arab if nc_arab > 0 else 0.0

    # ── Print results ─────────────────────────────────────────────────────
    print(f"\n  (a) Binomial (control vs 4% null):  p = {p_binom_c:.4g}")
    print(f"  (b) 2-prop z-test (UNESCO > ctrl):  z = {z:.3f}, p = {p_z:.4f}")
    print(f"  (c) OR = {or_global:.2f}  (95% CI {or_lo:.2f}--{or_hi:.2f})")
    print(f"  (d) Permutation p ({N_PERM:,} draws):  p = {p_perm:.3f}")
    print(f"  (e) Mantel-Haenszel: OR = {mh_or:.3f}, χ² = {mh_chi2:.3f}, p = {mh_p:.4f}")
    print(f"  (f) Drop-Europe OR = {or_noeu:.2f}")
    print(f"  (g) Arab States (P1435 ctrl): N = {nc_arab:,}, A+ = {kc_arab}, rate = {rc_arab*100:.2f}%")

    # ── LaTeX macros ──────────────────────────────────────────────────────
    def fmt_p(p: float, decimals: int = 3) -> str:
        """Format a p-value for LaTeX: use '$< 0.001$' when rounding would give 0.000."""
        threshold = 10 ** (-decimals)
        if p < threshold:
            zeros = "0" * decimals
            return f"$< 0.{zeros[:-1]}1$"
        return f"{p:.{decimals}f}"

    print(f"\n  % LaTeX macros (GROUP 23):")
    # Format N with comma
    n_c_fmt = f"{n_c:,}".replace(",", "{,}")
    print(f"  \\newcommand{{\\WDctrlN}}{{{n_c_fmt}}}           % Wikidata P1435 control corpus size")
    print(f"  \\newcommand{{\\WDctrlAp}}{{{k_c}}}            % A+ hits in Wikidata P1435 control")
    print(f"  \\newcommand{{\\WDctrlRate}}{{{r_c*100:.2f}}}          % A+ rate, Wikidata P1435 control (%)")
    print(f"  \\newcommand{{\\WDctrlBinomP}}{{{fmt_p(p_binom_c)}}}          % p-value, binomial test P1435 control")
    print(f"  \\newcommand{{\\WDctrlOR}}{{{or_global:.2f}}}           % global odds ratio, P1435 control")
    print(f"  \\newcommand{{\\WDctrlORlo}}{{{or_lo:.2f}}}           % OR 95% CI lower bound")
    print(f"  \\newcommand{{\\WDctrlORhi}}{{{or_hi:.2f}}}           % OR 95% CI upper bound")
    print(f"  \\newcommand{{\\WDctrlZp}}{{{p_z:.4f}}}          % p-value, Z-test OR significance")
    print(f"  \\newcommand{{\\WDctrlPermP}}{{{p_perm:.3f}}}          % permutation p, Wikidata P1435 control")
    print(f"  \\newcommand{{\\WDctrlMHOR}}{{{mh_or:.3f}}}           % Mantel-Haenszel OR, stratified")
    print(f"  \\newcommand{{\\WDctrlMHchi}}{{{mh_chi2:.3f}}}          % Mantel-Haenszel chi-sq")
    print(f"  \\newcommand{{\\WDctrlMHp}}{{{mh_p:.4f}}}           % p-value, Mantel-Haenszel test")
    print(f"  \\newcommand{{\\WDctrlDropEuropeOR}}{{{or_noeu:.2f}}}     % OR after dropping European sites")
    # Arab States region — descriptive only (P1435 control)
    nc_arab_fmt = f"{nc_arab:,}".replace(",", "{,}")
    print(f"  \\newcommand{{\\WDctrlArabN}}{{{nc_arab_fmt}}}           % P1435 control: Arab States site count")
    print(f"  \\newcommand{{\\WDctrlArabAp}}{{{kc_arab}}}            % P1435 control: Arab States A+ count")
    print(f"  \\newcommand{{\\WDctrlArabRate}}{{{rc_arab*100:.2f}}}          % P1435 control: Arab States A+ rate (%)")

    print(f"\n" + "=" * 72)
    print(f"  DONE — all GROUP 23 macros printed above.")
    print("=" * 72)


if __name__ == "__main__":
    main()
