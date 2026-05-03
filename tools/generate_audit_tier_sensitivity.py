"""
generate_audit_tier_sensitivity.py
====================================
Produce supplementary/audit/tier_sensitivity_audit.txt

Comprehensive tier-threshold sensitivity audit:
  - 8 subsets × 6 tier levels at primary null rates  (Part 1)
  - A+ null-rate sweep 1%–10% for all subsets        (Part 2)

Also writes machine-readable CSVs:
  supplementary/audit/tier_sensitivity_part1.csv
  supplementary/audit/tier_sensitivity_part2.csv

Run from repo root:
    python3 tools/generate_audit_tier_sensitivity.py
"""

import csv
import re
import sys
import json
from pathlib import Path
from datetime import datetime, timezone
from scipy.stats import binomtest, fisher_exact

sys.path.insert(0, str(Path(__file__).parent.parent))

from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, HARMONIC_STEP,
    P_NULL_APP, P_NULL_AP, P_NULL_A,
    deviation, tier_label,
    is_aplusplus, is_aplus, is_a_or_better,
    is_c_or_better, is_cminus_or_better,
)
from lib.dome_filter import is_dome_site
from lib.founding_filter import classify_site
from lib.beru import load_religion_sets

OUT_TXT  = Path(__file__).parent.parent / "supplementary" / "audit" / "tier_sensitivity_audit.txt"
OUT_CSV1 = Path(__file__).parent.parent / "supplementary" / "audit" / "tier_sensitivity_part1.csv"
OUT_CSV2 = Path(__file__).parent.parent / "supplementary" / "audit" / "tier_sensitivity_part2.csv"

_cfg = json.loads((Path(__file__).parent.parent / "config.json").read_text())

HARMONIC_STEP_DEG = HARMONIC_STEP * BERU   # 3.0°
HALF_STEP          = HARMONIC_STEP_DEG / 2  # 1.5°

NULL_RATES = [r / 1000 for r in range(10, 105, 5)]  # 0.010 … 0.100

PRIMARY_NR = {
    "A++": _cfg["tiers"]["A++"]["null_rate"],
    "A+":  _cfg["tiers"]["A+"]["null_rate"],
    "A":   _cfg["tiers"]["A"]["null_rate"],
    "C":   _cfg["tiers"]["C"]["null_rate"],
    "C-":  _cfg["tiers"]["C-"]["null_rate"],
    "C--": _cfg["tiers"]["C--"]["null_rate"],
}
TIER_NAMES = ["A++", "A+", "A", "C", "C-", "C--"]
NR_KEYS    = ["App", "Ap", "A", "C", "Cm", "Cmm"]
NR_VALS    = [PRIMARY_NR["A++"], PRIMARY_NR["A+"], PRIMARY_NR["A"],
              PRIMARY_NR["C"],   PRIMARY_NR["C-"], PRIMARY_NR["C--"]]

SEP  = "─" * 110
SEP2 = "═" * 110

# ── Helpers ──────────────────────────────────────────────────────────────────
def dev_deg(lon):
    arc = abs(lon - GERIZIM) % HARMONIC_STEP_DEG
    if arc > HALF_STEP:
        arc = HARMONIC_STEP_DEG - arc
    return arc

def mid_deg(lon):
    return HALF_STEP - dev_deg(lon)

def binom_p(k, n, p0):
    if n == 0: return 1.0
    return binomtest(k, n, p0, alternative="greater").pvalue

def fisher_p(k, n, p0):
    if n == 0: return 1.0
    ek = round(n * p0); enk = n - ek
    if ek == 0 or enk == 0: return 1.0
    _, p = fisher_exact([[k, n-k], [ek, enk]], alternative="greater")
    return p

def OR(k, n, p0):
    if n == 0 or p0 == 0: return float("nan")
    return (k / n) / p0

def sig(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"

def count_tier(devs_mids, p0_app, p0_ap, p0_a):
    t_app = p0_app * HALF_STEP
    t_ap  = p0_ap  * HALF_STEP
    t_a   = p0_a   * HALF_STEP
    counts = {t: 0 for t in TIER_NAMES}
    for d, m in devs_mids:
        if d <= t_app: counts["A++"] += 1
        if d <= t_ap:  counts["A+"]  += 1
        if d <= t_a:   counts["A"]   += 1
        if m <= t_a:   counts["C"]   += 1
        if m <= t_ap:  counts["C-"]  += 1
        if m <= t_app: counts["C--"] += 1
    return counts

# ── Load corpus ───────────────────────────────────────────────────────────────
print("Loading corpus …")
corpus   = load_corpus()
cultural = cultural_sites_with_coords(corpus)

RELIGION_SETS = load_religion_sets()
_relig_kws = [kw for _, kws in RELIGION_SETS for kw in kws]
_STUPA_RE  = re.compile(r"\bstupas?\b", re.IGNORECASE)
_MOUND_RES = [re.compile(p, re.IGNORECASE) for p in
              [r"\btumulus\b", r"\btumuli\b", r"\bbarrow\b",
               r"\bbarrows\b", r"\bkofun\b"]]

print("Classifying sites …")
_sites = []
for s in cultural:
    dd = dev_deg(s.longitude)
    md_ = mid_deg(s.longitude)
    ft  = s.full_text or ""
    is_dome_  = is_dome_site(s)
    is_mound_ = is_dome_ or any(r.search(ft) for r in _MOUND_RES)
    cats = classify_site(s)
    _sites.append({
        "name":       s.site,
        "lon":        s.longitude,
        "dev":        dd,
        "mid":        md_,
        "is_relig":   any(kw in ft for kw in _relig_kws),
        "is_dome":    is_dome_,
        "is_mound":   is_mound_,
        "is_stupa_kw": bool(_STUPA_RE.search(ft)),
        "is_founding": len(cats) > 0,
    })

N_FULL = len(_sites)

# Wikidata Q180987
_wiki_csv = Path(__file__).parent.parent / "data" / "store" / "unesco" / "wikidata_stupas_q180987.csv"
_wiki_lons = []
if _wiki_csv.exists():
    with open(_wiki_csv) as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    for row in csv.DictReader(lines):
        try:
            _wiki_lons.append(float(row["lon"]))
        except (ValueError, KeyError):
            pass
N_WIKI = len(_wiki_lons)

def _dm(lon):
    return (dev_deg(lon), mid_deg(lon))

SUBSETS = [
    ("full",       f"Full corpus (N={N_FULL})",
     [(s["dev"], s["mid"]) for s in _sites]),
    ("religion",   f"Religion-founding (N={sum(s['is_relig'] for s in _sites)})",
     [(s["dev"], s["mid"]) for s in _sites if s["is_relig"]]),
    ("dome",       f"Dome/tholos/stupa (N={sum(s['is_dome'] for s in _sites)})",
     [(s["dev"], s["mid"]) for s in _sites if s["is_dome"]]),
    ("mound",      f"Dome + mound (N={sum(s['is_mound'] for s in _sites)})",
     [(s["dev"], s["mid"]) for s in _sites if s["is_mound"]]),
    ("founding",   f"Founding sites (N={sum(s['is_founding'] for s in _sites)})",
     [(s["dev"], s["mid"]) for s in _sites if s["is_founding"]]),
    ("stupa_kw",   f"UNESCO stupa keyword (N={sum(s['is_stupa_kw'] for s in _sites)})",
     [(s["dev"], s["mid"]) for s in _sites if s["is_stupa_kw"]]),
    ("stupa_q180", f"Wikidata Q180987 stupa (N={N_WIKI})",
     [_dm(lon) for lon in _wiki_lons]),
    ("stupa_geo",  f"Wikidata stupa Java 107-115°E (N={sum(1 for lon in _wiki_lons if 107<=lon<=115)})",
     [_dm(lon) for lon in _wiki_lons if 107 <= lon <= 115]),
]

# ── Build output ───────────────────────────────────────────────────────────────
lines_out = []
w = lines_out.append

ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
w(SEP2)
w("TIER-THRESHOLD SENSITIVITY AUDIT")
w(f"Generated: {ts}")
w(f"Anchor: Gerizim {GERIZIM}°E  |  Harmonic step: {HARMONIC_STEP_DEG}°  |  Half-step: {HALF_STEP}°")
w(f"UNESCO Cultural+Mixed N={N_FULL}  |  Wikidata Q180987 N={N_WIKI}")
w(SEP2)

# ── PART 1 ─────────────────────────────────────────────────────────────────────
w("")
w("PART 1 — PRIMARY NULL RATES: all subsets × all tier levels")
w("          Primary null rates: A++/C--=2%  A+/C-=4%  A/C=8%")
w(SEP)
col_hdr = (
    f"{'Subset':<32}  {'N':>5}  "
    + "  ".join(
        f"{'tier':>4} {'p0':>4} {'n':>4} {'obs%':>5} {'OR':>5} {'p_binom':>9} {'sig':>3} {'p_fish':>9} {'sig':>3}"
        for tier in TIER_NAMES
    )
)
# Simpler column header
w(f"{'Subset':<32}  {'N':>5}  "
  f"  {'── A++ ──':^50}  {'── A+ ──':^50}  {'── A ──':^50}"
  f"  {'── C ──':^50}  {'── C- ──':^50}  {'── C-- ──':^50}")
w(f"{'':32}  {'':5}  "
  + ("  " + f"{'p0':>4} {'n':>4} {'obs%':>5} {'OR':>5} {'binom_p':>9} {'sig':>3} {'fish_p':>9} {'sig':>3}") * 6)
w(SEP)

csv1_rows = []
csv1_fields = ["subset", "label", "N"]
for t, k in zip(TIER_NAMES, NR_KEYS):
    csv1_fields += [f"p0_{k}", f"n_{k}", f"obs_pct_{k}", f"OR_{k}", f"binom_{k}", f"sig_{k}", f"fish_{k}", f"fsig_{k}"]

for key, label, dm_list in SUBSETS:
    N = len(dm_list)
    if N == 0:
        continue
    counts = count_tier(dm_list, PRIMARY_NR["A++"], PRIMARY_NR["A+"], PRIMARY_NR["A"])
    row = {"subset": key, "label": label, "N": N}
    cells = []
    for tier, nr_key, p0 in zip(TIER_NAMES, NR_KEYS, NR_VALS):
        k = counts[tier]
        bp  = binom_p(k, N, p0)
        fp  = fisher_p(k, N, p0)
        or_ = OR(k, N, p0)
        obs = k / N * 100
        cells.append(f"{p0*100:>4.1f}% {k:>4} {obs:>5.1f}% {or_:>5.2f} {bp:>9.2e} {sig(bp):>3} {fp:>9.2e} {sig(fp):>3}")
        row[f"p0_{nr_key}"]   = p0
        row[f"n_{nr_key}"]    = k
        row[f"obs_pct_{nr_key}"] = round(obs, 2)
        row[f"OR_{nr_key}"]   = round(or_, 3)
        row[f"binom_{nr_key}"] = f"{bp:.3e}"
        row[f"sig_{nr_key}"]  = sig(bp)
        row[f"fish_{nr_key}"] = f"{fp:.3e}"
        row[f"fsig_{nr_key}"] = sig(fp)
    csv1_rows.append(row)
    w(f"{label:<32}  {N:>5}  " + "  ".join(cells))

w(SEP)

# ── PART 2 ─────────────────────────────────────────────────────────────────────
w("")
w("PART 2 — A+ NULL-RATE SWEEP (1%–10%), all subsets")
w("         Columns: A+ threshold hits  AND  C (anti-harmonic mirror) hits")
w("         OR > 1.0 and p < 0.05 on A+ = signal; OR ≈ 1.0 on C = specificity")
w(SEP)

csv2_rows = []
csv2_fields = ["subset", "label", "N", "p0_pct", "thresh_deg", "thresh_km",
               "n_Ap", "obs_pct_Ap", "OR_Ap", "binom_Ap", "sig_Ap", "fish_Ap",
               "n_C",  "obs_pct_C",  "OR_C",  "binom_C",  "sig_C"]

for key, label, dm_list in SUBSETS:
    N = len(dm_list)
    if N == 0:
        continue
    w("")
    w(f"  {label}")
    w(f"  {'P0 (%)':>7}  {'thresh°':>8}  {'~km':>6}  "
      f"{'A+ n':>5} {'obs%':>6} {'OR':>5} {'binom_p':>10} {'sig':>4}  "
      f"{'C  n':>5} {'obs%':>6} {'OR':>5} {'binom_p':>10} {'sig':>4}")
    w(f"  {'─'*7}  {'─'*8}  {'─'*6}  " + "─" * 60)
    for p0 in NULL_RATES:
        thresh = p0 * HALF_STEP
        k_ap = sum(1 for d, _ in dm_list if d <= thresh)
        k_c  = sum(1 for _, m in dm_list if m <= thresh)
        bp_ap = binom_p(k_ap, N, p0)
        bp_c  = binom_p(k_c,  N, p0)
        fp_ap = fisher_p(k_ap, N, p0)
        or_ap = OR(k_ap, N, p0)
        or_c  = OR(k_c,  N, p0)
        flag  = " <-- primary" if abs(p0 - PRIMARY_NR["A+"]) < 1e-9 else ""
        w(f"  {p0*100:>7.1f}  {thresh:>8.4f}  {thresh*111:>6.1f}  "
          f"{k_ap:>5} {k_ap/N*100:>5.1f}% {or_ap:>5.2f} {bp_ap:>10.2e} {sig(bp_ap):>4}  "
          f"{k_c:>5} {k_c/N*100:>5.1f}% {or_c:>5.2f} {bp_c:>10.2e} {sig(bp_c):>4}{flag}")
        csv2_rows.append({
            "subset": key, "label": label, "N": N,
            "p0_pct": round(p0 * 100, 1),
            "thresh_deg": round(thresh, 4),
            "thresh_km":  round(thresh * 111, 1),
            "n_Ap":       k_ap,
            "obs_pct_Ap": round(k_ap / N * 100, 2),
            "OR_Ap":      round(or_ap, 3),
            "binom_Ap":   f"{bp_ap:.3e}",
            "sig_Ap":     sig(bp_ap),
            "fish_Ap":    f"{fp_ap:.3e}",
            "n_C":        k_c,
            "obs_pct_C":  round(k_c / N * 100, 2),
            "OR_C":       round(or_c, 3),
            "binom_C":    f"{bp_c:.3e}",
            "sig_C":      sig(bp_c),
        })

w("")
w(SEP2)
w("END OF AUDIT")
w(SEP2)

# ── Write outputs ─────────────────────────────────────────────────────────────
OUT_TXT.write_text("\n".join(lines_out) + "\n")
print(f"Wrote {OUT_TXT}")

with open(OUT_CSV1, "w", newline="") as f:
    wtr = csv.DictWriter(f, fieldnames=csv1_fields)
    wtr.writeheader()
    wtr.writerows(csv1_rows)
print(f"Wrote {OUT_CSV1}")

with open(OUT_CSV2, "w", newline="") as f:
    wtr = csv.DictWriter(f, fieldnames=csv2_fields)
    wtr.writeheader()
    wtr.writerows(csv2_rows)
print(f"Wrote {OUT_CSV2}")

print("Done.")
