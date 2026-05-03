"""
tier_threshold_audit.py
=======================
Sensitivity audit: sweep a range of null-rate percentages for the A+ tier
(and symmetrically for C-) and report binomial p-values, odds ratios, and
site counts for the full UNESCO Cultural+Mixed corpus.

The primary analysis uses A+ = 4 % (0.060°, ≈7 km).  This script sweeps
null rates from 1 % to 10 % in 0.5 % steps, recalculating the tier threshold
as  threshold_deg = null_rate × 1.5°  (half the 3° harmonic spacing).

Results are written to:
  results/tier_threshold_audit.csv
  results/tier_threshold_audit_macros.tex  (LaTeX table fragment)

Run from repo root:
    python3 analysis/global/tier_threshold_audit.py
"""

import sys
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, fisher_exact

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, HARMONIC_STEP

# ── Config ─────────────────────────────────────────────────────────────────
HARMONIC_STEP_DEG = HARMONIC_STEP * BERU   # 3.0°
HALF_STEP          = HARMONIC_STEP_DEG / 2  # 1.5°

# Sweep: null rates 1 % … 10 % in 0.5 % steps
NULL_RATES = [r / 1000 for r in range(10, 105, 5)]  # 0.010 … 0.100

# Primary null rate from config (4 %)
_cfg = json.loads((ROOT / "config.json").read_text())
PRIMARY_NULL_RATE = _cfg["tiers"]["A+"]["null_rate"]


# ── Helper ─────────────────────────────────────────────────────────────────
def harmonic_deviation_deg(lon: float, anchor: float = GERIZIM) -> float:
    """Deviation from nearest 3° harmonic line, in degrees."""
    arc = abs(lon - anchor) % HARMONIC_STEP_DEG
    if arc > HALF_STEP:
        arc = HARMONIC_STEP_DEG - arc
    return arc


# ── Load corpus ────────────────────────────────────────────────────────────
corpus = load_corpus()
sites  = cultural_sites_with_coords(corpus)
N      = len(sites)

devs = np.array([harmonic_deviation_deg(s.longitude) for s in sites])

# ── Sweep ──────────────────────────────────────────────────────────────────
rows = []
for p0 in NULL_RATES:
    thresh_deg = p0 * HALF_STEP          # half-step × null_rate
    n_in       = int((devs <= thresh_deg).sum())
    n_out      = N - n_in
    binom_p    = binomtest(n_in, N, p0, alternative="greater").pvalue
    # Fisher exact: 2×2 table  [in/out] × [observed/expected under null]
    expected_in  = round(N * p0)
    expected_out = N - expected_in
    _, fisher_p  = fisher_exact(
        [[n_in, n_out], [expected_in, expected_out]],
        alternative="greater"
    )
    odds = (n_in / N) / p0 if p0 > 0 else float("inf")
    flag = "← primary" if abs(p0 - PRIMARY_NULL_RATE) < 1e-9 else ""
    rows.append(dict(
        null_rate_pct = round(p0 * 100, 1),
        thresh_deg    = round(thresh_deg, 4),
        thresh_km     = round(thresh_deg * 111, 1),
        n_in          = n_in,
        N             = N,
        observed_rate = round(n_in / N, 4),
        odds_ratio    = round(odds, 2),
        binom_p       = f"{binom_p:.2e}",
        fisher_p      = f"{fisher_p:.2e}",
        flag          = flag,
    ))

# ── Write CSV ──────────────────────────────────────────────────────────────
out_csv = ROOT / "results" / "tier_threshold_audit.csv"
out_csv.parent.mkdir(exist_ok=True)
fieldnames = ["null_rate_pct","thresh_deg","thresh_km","n_in","N",
              "observed_rate","odds_ratio","binom_p","fisher_p","flag"]
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)
print(f"Wrote {out_csv}")

# ── Write LaTeX table fragment ─────────────────────────────────────────────
out_tex = ROOT / "results" / "tier_threshold_audit_macros.tex"
lines = [
    r"\begin{tabular}{rrrrrrrr}",
    r"\toprule",
    r"$P_0$ (\%) & $\delta_{\max}$ (\textdegree{}) & $\approx$km "
    r"& $n$ & Obs.\,rate & OR & Binom.\,$p$ & Fisher\,$p$ \\",
    r"\midrule",
]
for r in rows:
    bold = r["flag"] != ""
    prefix = r"\textbf{" if bold else ""
    suffix = "}" if bold else ""
    row_str = (
        f"{prefix}{r['null_rate_pct']:.1f}{suffix} & "
        f"{r['thresh_deg']:.4f} & "
        f"{r['thresh_km']:.1f} & "
        f"{r['n_in']} & "
        f"{float(r['observed_rate']):.3f} & "
        f"{r['odds_ratio']:.2f} & "
        r"\texttt{" + r['binom_p'] + r"} & "
        r"\texttt{" + r['fisher_p'] + r"} \\"
    )
    lines.append(row_str)
lines += [r"\bottomrule", r"\end{tabular}"]
out_tex.write_text("\n".join(lines) + "\n")
print(f"Wrote {out_tex}")

# ── Console summary ────────────────────────────────────────────────────────
print(f"\nN = {N}  |  Harmonic step = {HARMONIC_STEP_DEG}°  |  Primary P0 = {PRIMARY_NULL_RATE*100:.1f}%")
print(f"{'P0 (%)':>7}  {'thresh°':>8}  {'~km':>6}  {'n_in':>5}  {'obs%':>6}  {'OR':>5}  {'binom_p':>10}  {'fisher_p':>10}")
for r in rows:
    marker = " <-- primary" if r["flag"] else ""
    print(
        f"{r['null_rate_pct']:>7.1f}  "
        f"{r['thresh_deg']:>8.4f}  "
        f"{r['thresh_km']:>6.1f}  "
        f"{r['n_in']:>5}  "
        f"{float(r['observed_rate'])*100:>5.1f}%  "
        f"{r['odds_ratio']:>5.2f}  "
        f"{r['binom_p']:>10}  "
        f"{r['fisher_p']:>10}"
        f"{marker}"
    )
