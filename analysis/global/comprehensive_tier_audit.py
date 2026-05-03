"""
comprehensive_tier_audit.py
============================
Full tier-threshold sensitivity audit across all named subsets and all
tier levels (A++, A+, A, C, C-, C--).

For each subset × null-rate combination the script computes:
  - site count n
  - hits at each tier level (A++, A+, A, C, C-, C--)
  - one-sided binomial p-value against the geometric null
  - odds ratio (observed rate / null rate)
  - Fisher exact p-value (observed vs expected under null)

SUBSETS TESTED
--------------
  full        Full UNESCO Cultural+Mixed corpus (N ≈ 1011)
  religion    Sites mentioning world-religion founding (Buddhism, Christianity,
              Judaism, Islam, Zoroastrianism, Hinduism)
  dome        Dome/tholos/stupa/spherical UNESCO sites (dome_filter)
  mound       Dome + tumulus/barrow/kofun sites (dome_filter mound extension)
  founding    Sites classified by founding_filter (any category F/S/M/X/L)
  stupa_kw    UNESCO sites matching stupa/stupas keyword only
  stupa_q180  Wikidata Q180987 stupa corpus (external, non-UNESCO)
  stupa_geo   Stupa sites in Java-region longitude window (107-115°E)

TIER SWEEP
----------
  Null rates 1 %-10 % in 0.5 % steps (18 steps).
  Primary null rates: A++ = 2%, A+ = 4%, A = 8% (from config.json).
  C-tier mirrors: C = 8%, C- = 4%, C-- = 2%.

OUTPUT
------
  results/comprehensive_tier_audit.csv   — full matrix
  results/comprehensive_tier_audit.txt   — human-readable summary

Run from repo root:
    python3 analysis/global/comprehensive_tier_audit.py
"""

import csv
import sys
import json
from pathlib import Path
from scipy.stats import binomtest, fisher_exact

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, HARMONIC_STEP,
    TIER_APP_DEG, TIER_APLUS_DEG, TIER_A_DEG,
    TIER_C_DEG, TIER_CMINUS_DEG, TIER_CMINUS2_DEG,
    P_NULL_APP, P_NULL_AP, P_NULL_A,
    deviation, tier_label,
    is_aplusplus, is_aplus, is_a_or_better,
    is_c_or_better, is_cminus_or_better,
)
from lib.dome_filter import is_dome_site, is_dome_site_raw, FORM_KEYWORDS
from lib.founding_filter import classify_site

# ── Config ──────────────────────────────────────────────────────────────────
_cfg = json.loads((ROOT / "config.json").read_text())
HARMONIC_STEP_DEG = HARMONIC_STEP * BERU   # 3.0°
HALF_STEP          = HARMONIC_STEP_DEG / 2  # 1.5°
MIDPOINT_DEV_DEG   = HALF_STEP              # 1.5° from harmonic = perfect midpoint

# Null-rate sweep: 1 % to 10 % in 0.5 % steps
NULL_RATES = [r / 1000 for r in range(10, 105, 5)]

# Primary fixed null rates from config
PRIMARY_NR = {
    "A++": _cfg["tiers"]["A++"]["null_rate"],
    "A+":  _cfg["tiers"]["A+"]["null_rate"],
    "A":   _cfg["tiers"]["A"]["null_rate"],
    "C":   _cfg["tiers"]["C"]["null_rate"],
    "C-":  _cfg["tiers"]["C-"]["null_rate"],
    "C--": _cfg["tiers"]["C--"]["null_rate"],
}

# ── Helpers ─────────────────────────────────────────────────────────────────
def dev_deg(lon: float) -> float:
    """Deviation from nearest 3° harmonic in DEGREES."""
    arc = abs(lon - GERIZIM) % HARMONIC_STEP_DEG
    if arc > HALF_STEP:
        arc = HARMONIC_STEP_DEG - arc
    return arc


def midpoint_dev_deg(lon: float) -> float:
    """Distance from the inter-harmonic midpoint in DEGREES."""
    d = dev_deg(lon)
    return HALF_STEP - d  # 0 = perfect midpoint, HALF_STEP = on harmonic


def binom_one_sided(k: int, n: int, p0: float) -> float:
    if n == 0:
        return 1.0
    return binomtest(k, n, p0, alternative="greater").pvalue


def fisher_one_sided(k: int, n: int, p0: float) -> float:
    if n == 0:
        return 1.0
    exp_k = round(n * p0)
    exp_nk = n - exp_k
    if exp_k == 0 or exp_nk == 0:
        return 1.0
    _, p = fisher_exact([[k, n - k], [exp_k, exp_nk]], alternative="greater")
    return p


def odds_ratio(k: int, n: int, p0: float) -> float:
    if n == 0 or p0 == 0:
        return float("nan")
    return (k / n) / p0


def sig(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"


# ── Load UNESCO corpus ───────────────────────────────────────────────────────
print("Loading UNESCO corpus …")
corpus = load_corpus()
cultural = cultural_sites_with_coords(corpus)

# Pre-compute deviations for all cultural sites
_sites_data = []
for s in cultural:
    dd = dev_deg(s.longitude)
    md = midpoint_dev_deg(s.longitude)
    _sites_data.append({
        "obj":     s,
        "lon":     s.longitude,
        "dev_deg": dd,
        "mid_deg": md,
        "tier":    tier_label(deviation(s.longitude)),
        "text":    s.full_text or "",
    })

N_FULL = len(_sites_data)
print(f"  Full corpus: {N_FULL} sites")

# ── Religion keywords ────────────────────────────────────────────────────────
from lib.beru import load_religion_sets
RELIGION_SETS = load_religion_sets()
_relig_kws_all = [kw for _, kws in RELIGION_SETS for kw in kws]

def _is_relig(sd):
    return any(kw in sd["text"] for kw in _relig_kws_all)

# ── Dome / mound keyword sets (dome_filter) ──────────────────────────────────
import re
_STUPA_RES = re.compile(r"\bstupas?\b", re.IGNORECASE)

def _has_dome_keyword(sd):
    return is_dome_site(sd["obj"])

def _has_mound_keyword(sd):
    if _has_dome_keyword(sd):
        return True
    # tumulus / barrow / kofun — always unambiguous in monumental context
    for pat_str in [r"\btumulus\b", r"\btumuli\b", r"\bbarrow\b",
                    r"\bbarrows\b", r"\bkofun\b"]:
        if re.search(pat_str, sd["text"], re.IGNORECASE):
            return True
    return False

def _has_stupa_kw(sd):
    return bool(_STUPA_RES.search(sd["text"]))

def _has_founding(sd):
    cats = classify_site(sd["obj"])
    return len(cats) > 0

# Precompute flags (expensive — ~1-2 s)
print("Classifying sites (dome/mound/founding) …")
for sd in _sites_data:
    sd["is_relig"]    = _is_relig(sd)
    sd["is_dome"]     = _has_dome_keyword(sd)
    sd["is_mound"]    = _has_mound_keyword(sd)
    sd["is_stupa_kw"] = _has_stupa_kw(sd)
    sd["is_founding"] = _has_founding(sd)

# ── Wikidata Q180987 stupa corpus ────────────────────────────────────────────
import csv as _csv_mod
WIKIDATA_CSV = ROOT / "data" / "store" / "unesco" / "wikidata_stupas_q180987.csv"
_wiki_lons = []
if WIKIDATA_CSV.exists():
    with open(WIKIDATA_CSV) as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    reader = _csv_mod.DictReader(lines)
    for row in reader:
        try:
            _wiki_lons.append(float(row["lon"]))
        except (ValueError, KeyError):
            pass
print(f"  Wikidata Q180987 stupa corpus: {len(_wiki_lons)} sites")

_wiki_data = [
    {"lon": lon, "dev_deg": dev_deg(lon), "mid_deg": midpoint_dev_deg(lon),
     "in_geo": 107.0 <= lon <= 115.0}
    for lon in _wiki_lons
]

# ── Subset definitions ────────────────────────────────────────────────────────
# Each entry: (key, label, list_of_devs_and_mids)
#   devs_mids = list of (dev_deg, mid_deg) tuples
UNESCO_SUBSETS = [
    ("full",       "Full corpus",           [(s["dev_deg"], s["mid_deg"]) for s in _sites_data]),
    ("religion",   "Religion founding",     [(s["dev_deg"], s["mid_deg"]) for s in _sites_data if s["is_relig"]]),
    ("dome",       "Dome/tholos/stupa",     [(s["dev_deg"], s["mid_deg"]) for s in _sites_data if s["is_dome"]]),
    ("mound",      "Dome + mound",          [(s["dev_deg"], s["mid_deg"]) for s in _sites_data if s["is_mound"]]),
    ("founding",   "Founding sites",        [(s["dev_deg"], s["mid_deg"]) for s in _sites_data if s["is_founding"]]),
    ("stupa_kw",   "UNESCO stupa keyword",  [(s["dev_deg"], s["mid_deg"]) for s in _sites_data if s["is_stupa_kw"]]),
]
WIKI_SUBSETS = [
    ("stupa_q180", "Wikidata Q180987 stupa", [(s["dev_deg"], s["mid_deg"]) for s in _wiki_data]),
    ("stupa_geo",  "Wikidata stupa (Java region 107-115°E)", [(s["dev_deg"], s["mid_deg"]) for s in _wiki_data if s["in_geo"]]),
]
ALL_SUBSETS = UNESCO_SUBSETS + WIKI_SUBSETS

for key, label, devs in ALL_SUBSETS:
    print(f"  {key:12s}  N={len(devs):4d}  {label}")

# ── Tier classification at arbitrary threshold ────────────────────────────────
TIER_NAMES = ["A++", "A+", "A", "C", "C-", "C--"]

def classify_at_p0(dev_deg_val: float, mid_deg_val: float,
                   p0_app: float, p0_ap: float, p0_a: float) -> dict:
    """
    Classify a site against given null rates (which map to degree thresholds).

    Returns dict: tier → bool
    A-tier thresholds: thresh = p0 × 1.5°
    C-tier thresholds: midpoint threshold = p0 × 1.5°  (symmetric)
    """
    t_app = p0_app * HALF_STEP
    t_ap  = p0_ap  * HALF_STEP
    t_a   = p0_a   * HALF_STEP
    return {
        "A++": dev_deg_val <= t_app,
        "A+":  dev_deg_val <= t_ap,
        "A":   dev_deg_val <= t_a,
        "C":   mid_deg_val <= t_a,
        "C-":  mid_deg_val <= t_ap,
        "C--": mid_deg_val <= t_app,
    }


def count_at_p0(devs_mids, p0_app, p0_ap, p0_a):
    counts = {t: 0 for t in TIER_NAMES}
    for dev_d, mid_d in devs_mids:
        c = classify_at_p0(dev_d, mid_d, p0_app, p0_ap, p0_a)
        for t in TIER_NAMES:
            if c[t]:
                counts[t] += 1
    return counts


# ── PART 1: Fixed primary null rates, all subsets, all tier levels ───────────
print("\n" + "=" * 100)
print("PART 1 — PRIMARY NULL RATES: all subsets × all tier levels")
print("=" * 100)

part1_rows = []
header = ["subset", "label", "N",
          "n_App", "p0_App", "binom_App", "sig_App", "OR_App",
          "n_Ap",  "p0_Ap",  "binom_Ap",  "sig_Ap",  "OR_Ap",
          "n_A",   "p0_A",   "binom_A",   "sig_A",   "OR_A",
          "n_C",   "p0_C",   "binom_C",   "sig_C",   "OR_C",
          "n_Cm",  "p0_Cm",  "binom_Cm",  "sig_Cm",  "OR_Cm",
          "n_Cmm", "p0_Cmm", "binom_Cmm", "sig_Cmm", "OR_Cmm"]

NR_APP = PRIMARY_NR["A++"]
NR_AP  = PRIMARY_NR["A+"]
NR_A   = PRIMARY_NR["A"]
NR_C   = PRIMARY_NR["C"]
NR_CM  = PRIMARY_NR["C-"]
NR_CMM = PRIMARY_NR["C--"]

fmt_hdr = (
    f"{'Subset':<18}  {'N':>5}  "
    f"{'A++ n':>6} {'p0':>5} {'p':>10} {'sig':>4} {'OR':>5}  "
    f"{'A+  n':>6} {'p0':>5} {'p':>10} {'sig':>4} {'OR':>5}  "
    f"{'A   n':>6} {'p0':>5} {'p':>10} {'sig':>4} {'OR':>5}  "
    f"{'C   n':>6} {'p0':>5} {'p':>10} {'sig':>4} {'OR':>5}  "
    f"{'C-  n':>6} {'p0':>5} {'p':>10} {'sig':>4} {'OR':>5}  "
    f"{'C-- n':>6} {'p0':>5} {'p':>10} {'sig':>4} {'OR':>5}"
)
print(fmt_hdr)
print("-" * len(fmt_hdr))

for key, label, devs_mids in ALL_SUBSETS:
    N = len(devs_mids)
    if N == 0:
        continue
    counts = count_at_p0(devs_mids, NR_APP, NR_AP, NR_A)
    row = {"subset": key, "label": label, "N": N}
    for tier_key, nr_key, p0 in [
        ("A++", "App",  NR_APP),
        ("A+",  "Ap",   NR_AP),
        ("A",   "A",    NR_A),
        ("C",   "C",    NR_C),
        ("C-",  "Cm",   NR_CM),
        ("C--", "Cmm",  NR_CMM),
    ]:
        k  = counts[tier_key]
        bp = binom_one_sided(k, N, p0)
        or_ = odds_ratio(k, N, p0)
        row[f"n_{nr_key}"]     = k
        row[f"p0_{nr_key}"]    = p0
        row[f"binom_{nr_key}"] = bp
        row[f"sig_{nr_key}"]   = sig(bp)
        row[f"OR_{nr_key}"]    = or_
    part1_rows.append(row)

    def _fmt(nr_key, p0):
        k  = row[f"n_{nr_key}"]
        bp = row[f"binom_{nr_key}"]
        s  = row[f"sig_{nr_key}"]
        o  = row[f"OR_{nr_key}"]
        return f"{k:>6} {p0*100:>4.1f}% {bp:>10.2e} {s:>4} {o:>5.2f}"

    line = (
        f"{label:<18}  {N:>5}  "
        f"{_fmt('App', NR_APP)}  "
        f"{_fmt('Ap',  NR_AP)}  "
        f"{_fmt('A',   NR_A)}  "
        f"{_fmt('C',   NR_C)}  "
        f"{_fmt('Cm',  NR_CM)}  "
        f"{_fmt('Cmm', NR_CMM)}"
    )
    print(line)

# ── PART 2: A+ sweep across null rates 1%-10%, all subsets ───────────────────
print("\n" + "=" * 100)
print("PART 2 — A+ NULL-RATE SWEEP (1%-10%), all subsets")
print("=" * 100)

part2_rows = []
p2_header = ["subset", "label", "N", "p0_pct", "thresh_deg", "thresh_km",
             "n_Ap", "obs_rate", "OR_Ap", "binom_Ap", "sig_Ap", "fisher_Ap",
             "n_C",  "OR_C",  "binom_C",  "sig_C"]

for key, label, devs_mids in ALL_SUBSETS:
    N = len(devs_mids)
    if N == 0:
        continue
    print(f"\n  {label} (N={N})")
    print(f"  {'P0 (%)':>7}  {'thresh°':>8}  {'~km':>6}  "
          f"{'A+':>4} {'obs%':>6} {'OR':>5} {'binom':>10} {'sig':>4} "
          f"{'C':>4} {'C_OR':>5} {'C_binom':>10} {'C_sig':>4}")
    for p0 in NULL_RATES:
        thresh = p0 * HALF_STEP
        k_ap   = sum(1 for d, _ in devs_mids if d <= thresh)
        k_c    = sum(1 for _, m in devs_mids if m <= thresh)
        bp_ap  = binom_one_sided(k_ap, N, p0)
        bp_c   = binom_one_sided(k_c,  N, p0)
        or_ap  = odds_ratio(k_ap, N, p0)
        or_c   = odds_ratio(k_c,  N, p0)
        fp_ap  = fisher_one_sided(k_ap, N, p0)
        print(f"  {p0*100:>7.1f}  {thresh:>8.4f}  {thresh*111:>6.1f}  "
              f"{k_ap:>4} {k_ap/N*100:>5.1f}% {or_ap:>5.2f} {bp_ap:>10.2e} {sig(bp_ap):>4} "
              f"{k_c:>4} {or_c:>5.2f} {bp_c:>10.2e} {sig(bp_c):>4}")
        part2_rows.append({
            "subset": key, "label": label, "N": N,
            "p0_pct": round(p0 * 100, 1),
            "thresh_deg": round(thresh, 4),
            "thresh_km": round(thresh * 111, 1),
            "n_Ap": k_ap,
            "obs_rate": round(k_ap / N, 4),
            "OR_Ap": round(or_ap, 3),
            "binom_Ap": f"{bp_ap:.3e}",
            "sig_Ap": sig(bp_ap),
            "fisher_Ap": f"{fp_ap:.3e}",
            "n_C":  k_c,
            "OR_C": round(or_c, 3),
            "binom_C": f"{bp_c:.3e}",
            "sig_C": sig(bp_c),
        })

# ── Write CSV outputs ────────────────────────────────────────────────────────
out_dir = ROOT / "results"
out_dir.mkdir(exist_ok=True)

csv1 = out_dir / "comprehensive_tier_audit_part1.csv"
with open(csv1, "w", newline="") as f:
    w = _csv_mod.DictWriter(f, fieldnames=header)
    w.writeheader()
    w.writerows(part1_rows)
print(f"\nWrote {csv1}")

csv2 = out_dir / "comprehensive_tier_audit_part2.csv"
with open(csv2, "w", newline="") as f:
    w = _csv_mod.DictWriter(f, fieldnames=p2_header)
    w.writeheader()
    w.writerows(part2_rows)
print(f"Wrote {csv2}")

# ── Write human-readable summary ─────────────────────────────────────────────
out_txt = out_dir / "comprehensive_tier_audit.txt"
with open(out_txt, "w") as f:
    f.write("COMPREHENSIVE TIER-THRESHOLD AUDIT\n")
    f.write("All subsets × all tier levels × null-rate sweep\n")
    f.write(f"Corpus: UNESCO Cultural+Mixed N={N_FULL}; "
            f"Wikidata Q180987 N={len(_wiki_lons)}\n")
    f.write("=" * 100 + "\n\n")
    f.write("PART 1 — Primary null rates\n")
    f.write(fmt_hdr + "\n")
    f.write("-" * len(fmt_hdr) + "\n")
    for row in part1_rows:
        def _fmt2(nr_key, p0):
            k  = row[f"n_{nr_key}"]
            bp = row[f"binom_{nr_key}"]
            s  = row[f"sig_{nr_key}"]
            o  = row[f"OR_{nr_key}"]
            return f"{k:>6} {p0*100:>4.1f}% {bp:>10.2e} {s:>4} {o:>5.2f}"
        N = row["N"]
        label = row["label"]
        line = (
            f"{label:<18}  {N:>5}  "
            f"{_fmt2('App', NR_APP)}  {_fmt2('Ap', NR_AP)}  "
            f"{_fmt2('A', NR_A)}  {_fmt2('C', NR_C)}  "
            f"{_fmt2('Cm', NR_CM)}  {_fmt2('Cmm', NR_CMM)}"
        )
        f.write(line + "\n")
    f.write("\n\nPART 2 — A+ sweep (see CSV for full matrix)\n")
    f.write("Significant (p<0.05) cells across the sweep:\n")
    sig_cells = [(r["subset"], r["p0_pct"], r["OR_Ap"], r["binom_Ap"], r["sig_Ap"])
                 for r in part2_rows if r["sig_Ap"] in ("*", "**", "***", "~")]
    for subset, p0, or_, bp, s in sig_cells:
        f.write(f"  {subset:<14} P0={p0:>4.1f}%  OR={or_:>5.3f}  p={bp}  {s}\n")
print(f"Wrote {out_txt}")

print("\nDone.")
