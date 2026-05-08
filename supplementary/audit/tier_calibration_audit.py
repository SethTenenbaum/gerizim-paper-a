"""
tier_calibration_audit.py
=========================
Auditable reproduction of the per-corpus natural proximity scales cited in the
Tier Classification section of the manuscript:

    "We define three proximity thresholds ... calibrated to bracket the
     empirical natural proximity scales of the corpora under study:
     the OWTRAD Silk Road network (X%), the UNESCO dome subset (X%),
     the full UNESCO Cultural/Mixed corpus (X%), and the Wikidata stupa
     corpus (X%)."

The 'natural proximity scale' for a corpus is the geometric-null rate p0
at which that corpus's one-sided binomial enrichment is most statistically
concentrated, i.e. the window where -log(p_binom) peaks across a
logarithmically-spaced grid of 15 window sizes.

This script:
  1. Loads all four corpora independently.
  2. Computes the binomial p-curve across all 15 windows for each.
  3. Identifies the peak window (natural scale) for each corpus.
  4. Prints a full table of hit counts, expected counts, ratios, and
     p-values at every window for every corpus.
  5. Prints the per-corpus natural scales and the macro values they
     correspond to, so anyone can verify the manuscript figures.

Run from the project root:
    python3 supplementary/audit/tier_calibration_audit.py

Output is purely printed; no files are written.
"""

import csv
import sys
import json
import math
import numpy as np
from pathlib import Path
from scipy.stats import binom
from scipy.special import betainc

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT))

from data.unesco_corpus import load_corpus
from lib.dome_filter import FORM_KEYWORDS, FORM_KEYWORD_RES

# ── Configuration ─────────────────────────────────────────────────────────────
_CFG          = json.loads((_ROOT / "config.json").read_text())
GERIZIM_LON   = _CFG["anchors"]["gerizim"]["longitude"]
HARMONIC_STEP = 3.0

# Identical 15-window grid used in multiscale_combined_p.py
NULL_RATES = np.geomspace(0.0033, 0.40, 15)
WINDOWS    = NULL_RATES * HARMONIC_STEP / 2.0
K          = len(WINDOWS)

# ── Helpers ───────────────────────────────────────────────────────────────────
def harmonic_devs(lons):
    shifted = (lons - GERIZIM_LON) % HARMONIC_STEP
    return np.minimum(shifted, HARMONIC_STEP - shifted)

def binom_p_curve(lons, N):
    devs   = harmonic_devs(lons)
    counts = np.array([np.sum(devs <= w) for w in WINDOWS])
    p_arr  = np.clip(np.array([
        binom.sf(obs - 1, N, p0)
        for obs, p0 in zip(counts, NULL_RATES)
    ]), 1e-300, 1.0)
    return counts, p_arr

# ── Load corpora ──────────────────────────────────────────────────────────────
corpus   = load_corpus()
cultural = [s for s in corpus if s.category != "Natural" and s.has_coords]

lons_full = np.array([s.longitude for s in cultural])
N_full    = len(lons_full)

lons_dome = np.array([
    s.longitude for s in cultural
    if any(FORM_KEYWORD_RES[kw].search(s.full_text) for kw in FORM_KEYWORDS)
])
N_dome = len(lons_dome)

_WIKI_CSV = _ROOT / "data" / "store" / "unesco" / "wikidata_stupas_q180987.csv"
_wiki_lons = []
if _WIKI_CSV.exists():
    with open(_WIKI_CSV, newline="") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    for row in csv.DictReader(lines):
        try:
            _wiki_lons.append(float(row["lon"]))
        except (ValueError, KeyError):
            continue
lons_wiki = np.array(_wiki_lons)
N_wiki    = len(lons_wiki)

_OWTRAD_NODES = _ROOT / "data" / "store" / "silk_road" / "owtrad_nodes.csv"
_owtrad_lons, _seen = [], set()
if _OWTRAD_NODES.exists():
    with open(_OWTRAD_NODES, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            try:
                lon = float(row["lon"])
                lat = float(row["lat"])
            except (ValueError, KeyError):
                continue
            key = (round(lon, 4), round(lat, 4))
            if key not in _seen:
                _seen.add(key)
                _owtrad_lons.append(lon)
lons_owtrad = np.array(_owtrad_lons)
N_owtrad    = len(lons_owtrad)

# ── Compute curves ────────────────────────────────────────────────────────────
obs_dome,   p_dome   = binom_p_curve(lons_dome,   N_dome)
obs_full,   p_full   = binom_p_curve(lons_full,   N_full)
obs_wiki,   p_wiki   = binom_p_curve(lons_wiki,   N_wiki)
obs_owtrad, p_owtrad = binom_p_curve(lons_owtrad, N_owtrad)

dome_peak_idx   = int(np.argmax(-np.log(p_dome)))
full_peak_idx   = int(np.argmax(-np.log(p_full)))
wiki_peak_idx   = int(np.argmax(-np.log(p_wiki)))
owtrad_peak_idx = int(np.argmax(-np.log(p_owtrad)))

# ── Print ──────────────────────────────────────────────────────────────────────
SEP  = "=" * 90
SEP2 = "-" * 90
HDR  = (f"  {'Null%':>7}  {'Window°':>9}  {'~km@30N':>8}  "
        f"{'Obs':>5}  {'Exp':>6}  {'Ratio':>6}  {'p_binom':>10}  {'Peak':>5}")

def print_corpus_table(label, N, obs_arr, p_arr, peak_idx):
    print()
    print(SEP)
    print(f"  CORPUS: {label}   (N = {N})")
    print(SEP)
    print(HDR)
    print("  " + SEP2)
    for i in range(K):
        exp   = N * NULL_RATES[i]
        ratio = obs_arr[i] / exp if exp > 0 else float("nan")
        km    = WINDOWS[i] * 111.0 * math.cos(math.radians(30))
        marker = "  <-- PEAK" if i == peak_idx else ""
        print(f"  {NULL_RATES[i]*100:>7.2f}%  {WINDOWS[i]:>9.4f}  {km:>8.1f}  "
              f"{obs_arr[i]:>5}  {exp:>6.1f}  {ratio:>6.2f}x  "
              f"{p_arr[i]:>10.6f}{marker}")
    print()
    print(f"  Natural scale: null rate = {NULL_RATES[peak_idx]*100:.1f}%  "
          f"({WINDOWS[peak_idx]:.4f} deg, "
          f"~{WINDOWS[peak_idx]*111.0*math.cos(math.radians(30)):.0f} km at 30N)")
    print()

print()
print(SEP)
print("  TIER CALIBRATION AUDIT")
print("  Reproduces per-corpus natural proximity scales cited in manuscript")
print(f"  Anchor: {GERIZIM_LON}°E  |  Period: {HARMONIC_STEP}°  |  Windows: {K} log-spaced")
print(SEP)

print_corpus_table("OWTRAD Silk Road vertices",          N_owtrad, obs_owtrad, p_owtrad, owtrad_peak_idx)
print_corpus_table("UNESCO dome/stupa subset",           N_dome,   obs_dome,   p_dome,   dome_peak_idx)
print_corpus_table("UNESCO Cultural/Mixed full corpus",  N_full,   obs_full,   p_full,   full_peak_idx)
print_corpus_table("Wikidata Q180987 stupa corpus",      N_wiki,   obs_wiki,   p_wiki,   wiki_peak_idx)

print()
print(SEP)
print("  SUMMARY: PER-CORPUS NATURAL SCALES")
print("  (Values cited in manuscript Tier Classification section)")
print(SEP)
print()
rows = [
    ("OWTRAD Silk Road",    N_owtrad, owtrad_peak_idx),
    ("UNESCO dome subset",  N_dome,   dome_peak_idx),
    ("UNESCO full corpus",  N_full,   full_peak_idx),
    ("Wikidata stupas",     N_wiki,   wiki_peak_idx),
]
print(f"  {'Corpus':<25}  {'N':>6}  {'Null rate':>10}  {'Window (deg)':>13}  "
      f"{'~km@30N':>9}  {'Grid index':>11}")
print("  " + "-" * 82)
for label, N, idx in rows:
    km = WINDOWS[idx] * 111.0 * math.cos(math.radians(30))
    print(f"  {label:<25}  {N:>6}  "
          f"{NULL_RATES[idx]*100:>9.1f}%  {WINDOWS[idx]:>13.4f}  "
          f"{km:>9.1f}  [{idx:>2}] of {K-1}")
print()
print("  Note: 'Grid index' is the position in np.geomspace(0.0033, 0.40, 15).")
print("  The manuscript rounds these to: OWTRAD 2.6%, dome 3.6%, full 10.2%,")
print("  Wikidata 20.2% -- corresponding to indices [6, 7, 10, 12].")

print()
print(SEP)
print("  LATEX MACROS (emit by running multiscale_combined_p.py)")
print(SEP)
print()
print(f"  \\naturalNullPctOwtrad  = {NULL_RATES[owtrad_peak_idx]*100:.1f}")
print(f"  \\naturalNullPctDome    = {NULL_RATES[dome_peak_idx]*100:.1f}")
print(f"  \\naturalNullPctFull    = {NULL_RATES[full_peak_idx]*100:.1f}")
print(f"  \\naturalNullPctWiki    = {NULL_RATES[wiki_peak_idx]*100:.1f}")
print()
