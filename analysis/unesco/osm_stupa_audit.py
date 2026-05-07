"""
osm_stupa_audit.py
==================
Rayleigh phase-concentration and harmonic-enrichment audit for the
OpenStreetMap stupa corpus — a second independently sourced stupa dataset
for portability confirmation.

Corpus : data/store/unesco/osm_stupas.csv
Fetch  : data/scripts/fetch_osm_stupas.py (Overpass API, man_made/historic/
         building/ruins = stupa, deduplicated at 0.01 deg)

Tags queried: man_made=stupa (primary), historic=stupa, building=stupa,
              ruins=stupa.  man_made=stupa accounts for ~95% of features.

TESTS
-----
1. Rayleigh phase-concentration permutation test at T=3 deg
   (same as wikidata_q180987_stupa_audit.py Test 3, primary result)
2. A-tier enrichment  : one-sided exact binomial vs 20% null
3. A+ enrichment      : one-sided exact binomial vs geometric null from config
4. Cluster asymmetry  : mean site count at A-bearing vs non-A harmonics
5. Shared-concentration Fisher combination with Wikidata and UNESCO full

All results written to ResultsStore and printed as LaTeX macros.
"""

import sys
import csv
import json
import numpy as np
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT))

from lib.results_store import ResultsStore
from lib.beru import deviation as beru_deviation
from lib.beru import tier_label as tier_from_deviation

_CFG      = json.loads((_ROOT / "config.json").read_text())
ANCHOR    = _CFG["anchors"]["gerizim"]["longitude"]
PERIOD    = 3.0          # 3° longitude grid spacing (1 beru)
N_PERM    = _CFG["simulation"]["n_permutations"]
SEED      = _CFG["simulation"]["random_seed"]
P_NULL_AP = _CFG["tiers"]["A+"]["null_rate"]   # 0.10
P_NULL_A  = _CFG["tiers"]["A"]["null_rate"]    # 0.20

AP_THRESH_DEG  = _CFG["tiers"]["A+"]["max_deviation_deg"]   # 0.15
APP_THRESH_DEG = _CFG["tiers"]["A++"]["max_deviation_deg"]  # 0.06
A_THRESH_DEG   = _CFG["tiers"]["A"]["max_deviation_deg"]    # 0.30

rng = np.random.default_rng(SEED)

OSM_CSV = _ROOT / "data" / "store" / "unesco" / "osm_stupas.csv"
WIKI_CSV = _ROOT / "data" / "store" / "unesco" / "wikidata_stupas_q180987.csv"

# ── Load OSM corpus ───────────────────────────────────────────────────────────

def load_csv_lons(path):
    lons = []
    with open(path, newline="", encoding="utf-8") as fh:
        lines = [l for l in fh if not l.startswith("#")]
        for row in csv.DictReader(lines):
            try:
                lons.append(float(row["lon"]))
            except (KeyError, ValueError):
                pass
    return np.array(lons)

lons_osm  = load_csv_lons(OSM_CSV)
lons_wiki = load_csv_lons(WIKI_CSV)
N_osm  = len(lons_osm)
N_wiki = len(lons_wiki)

print(f"\nOSM corpus loaded:      N = {N_osm}")
print(f"Wikidata corpus loaded: N = {N_wiki}")

# ── Tier classification ───────────────────────────────────────────────────────

def classify(lons):
    devs  = np.array([beru_deviation(lon, ANCHOR, PERIOD) for lon in lons])
    tiers = [tier_from_deviation(d) for d in devs]
    app   = sum(1 for t in tiers if t == "A++")
    ap    = sum(1 for t in tiers if t in ("A++", "A+"))
    a     = sum(1 for t in tiers if t in ("A++", "A+", "A"))
    return devs, tiers, app, ap, a

devs_osm, tiers_osm, n_app, n_ap, n_a = classify(lons_osm)

print(f"\nOSM tier counts:")
print(f"  A++: {n_app}  A+: {n_ap}  A: {n_a}  N={N_osm}")

# ── Rayleigh phase-concentration permutation test ────────────────────────────

def rayleigh_R(lons, period):
    phases = (2 * np.pi * lons / period) % (2 * np.pi)
    return float(np.abs(np.mean(np.exp(1j * phases))))

R_obs_osm = rayleigh_R(lons_osm, PERIOD)

perm_Rs = np.empty(N_PERM)
for i in range(N_PERM):
    shift = rng.uniform(0, PERIOD)
    perm_Rs[i] = rayleigh_R(lons_osm + shift, PERIOD)

rayleigh_p_osm = float(np.mean(perm_Rs >= R_obs_osm))

print(f"\nTest 1 — Rayleigh phase concentration (T={PERIOD}°):")
print(f"  R_obs = {R_obs_osm:.4f},  perm-p = {rayleigh_p_osm:.4f}")

# ── Binomial enrichment tests ────────────────────────────────────────────────

from scipy.stats import binomtest as _binomtest

def binom_p(k, n, p0):
    return _binomtest(k, n, p0, alternative="greater").pvalue

P_NULL_APP = _CFG["tiers"]["A++"]["null_rate"]  # 0.04

p_a_tier   = binom_p(n_a,   N_osm, P_NULL_A)
p_ap_tier  = binom_p(n_ap,  N_osm, P_NULL_AP)
p_app_tier = binom_p(n_app, N_osm, P_NULL_APP)

rate_a  = 100 * n_a  / N_osm
rate_ap = 100 * n_ap / N_osm
rate_app = 100 * n_app / N_osm

print(f"\nTest 2 — A-tier enrichment:")
print(f"  A: {n_a}/{N_osm} = {rate_a:.1f}%  p={p_a_tier:.4f}")
print(f"  A+: {n_ap}/{N_osm} = {rate_ap:.1f}%  p={p_ap_tier:.4f}")
print(f"  A++: {n_app}/{N_osm} = {rate_app:.1f}%  p={p_app_tier:.4f}")

# ── Cluster asymmetry ─────────────────────────────────────────────────────────

harmonics = np.arange(0, 360, PERIOD)
sites_per_harmonic = []
a_harmonic = []
for h in harmonics:
    in_window = np.sum(np.abs(((lons_osm - h + PERIOD/2) % PERIOD) - PERIOD/2) <= A_THRESH_DEG)
    has_ap = any(abs(((lon - h + PERIOD/2) % PERIOD) - PERIOD/2) <= AP_THRESH_DEG
                 for lon in lons_osm)
    sites_per_harmonic.append(in_window)
    a_harmonic.append(has_ap)

sites_per_harmonic = np.array(sites_per_harmonic, dtype=float)
a_harmonic = np.array(a_harmonic)

mean_ap_node    = float(np.mean(sites_per_harmonic[a_harmonic]))  if a_harmonic.any() else 0.0
mean_nonap_node = float(np.mean(sites_per_harmonic[~a_harmonic])) if (~a_harmonic).any() else 0.0
cluster_ratio   = mean_ap_node / mean_nonap_node if mean_nonap_node > 0 else -1.0

print(f"\nTest 3 — Cluster asymmetry:")
print(f"  A+-node mean={mean_ap_node:.2f}, non-A-node mean={mean_nonap_node:.2f}, ratio={cluster_ratio:.2f}")

# ── Shared-concentration Fisher combination: OSM + Wikidata + UNESCO full ────
# Load pre-computed Rayleigh p-values for Wikidata and UNESCO from store

store = ResultsStore()

# Wikidata Rayleigh p from store (set by wikidata_q180987_stupa_audit.py)
p_wiki_rayleigh = store.read("circConcStupaP")   # 0.0018
# UNESCO full Rayleigh p from store
p_unesco_rayleigh = store.read("circConcUnescoP")  # 0.603

if p_wiki_rayleigh is None or p_unesco_rayleigh is None:
    print("\n  WARNING: Wikidata or UNESCO Rayleigh p not found in store.")
    print("  Run circular_two_sample_tests.py and wikidata_q180987_stupa_audit.py first.")
    p_wiki_rayleigh   = p_wiki_rayleigh   or 0.0018
    p_unesco_rayleigh = p_unesco_rayleigh or 0.603

from scipy.stats import chi2

def fisher_combine(pvals):
    """Fisher's method: -2 * sum(log(p)) ~ chi2(2k)."""
    stat = -2.0 * sum(np.log(max(p, 1e-15)) for p in pvals)
    df   = 2 * len(pvals)
    return stat, float(chi2.sf(stat, df))

# Two-corpus combination: OSM + Wikidata
chi2_two, p_two_fisher = fisher_combine([rayleigh_p_osm, p_wiki_rayleigh])
# Three-corpus combination: OSM + Wikidata + UNESCO full
chi2_three, p_three_fisher = fisher_combine([rayleigh_p_osm, p_wiki_rayleigh, p_unesco_rayleigh])

print(f"\nTest 4 — Shared-concentration Fisher combination:")
print(f"  OSM Rayleigh p        = {rayleigh_p_osm:.4f}")
print(f"  Wikidata Rayleigh p   = {p_wiki_rayleigh:.4f}")
print(f"  UNESCO full Rayleigh p = {p_unesco_rayleigh:.4f}")
print(f"  OSM + Wikidata:     chi2({2*2}) = {chi2_two:.3f},   p = {p_two_fisher:.4f}")
print(f"  OSM+Wiki+UNESCO:    chi2({2*3}) = {chi2_three:.3f}, p = {p_three_fisher:.4f}")

# ── Tag breakdown summary ─────────────────────────────────────────────────────

tag_counts = {}
with open(OSM_CSV, newline="", encoding="utf-8") as fh:
    lines = [l for l in fh if not l.startswith("#")]
    for row in csv.DictReader(lines):
        tag = row.get("matched_tag", "unknown")
        tag_counts[tag] = tag_counts.get(tag, 0) + 1

n_man_made = tag_counts.get("man_made=stupa", 0)
n_historic  = tag_counts.get("historic=stupa", 0)
n_building  = tag_counts.get("building=stupa", 0)
n_ruins     = tag_counts.get("ruins=stupa", 0)

# ── Write to store ────────────────────────────────────────────────────────────

store.write("osmStupaTotal",       N_osm)
store.write("osmStupaRayleighR",   round(R_obs_osm, 4))
store.write("osmStupaRayleighP",   round(rayleigh_p_osm, 4))
store.write("osmStupaATierCount",  int(n_a))
store.write("osmStupaATierRate",   round(rate_a, 1))
store.write("osmStupaATierBinomP", round(p_a_tier, 4))
store.write("osmStupaApCount",     int(n_ap))
store.write("osmStupaApRate",      round(rate_ap, 1))
store.write("osmStupaApBinomP",    round(p_ap_tier, 4))
store.write("osmStupaAppCount",    int(n_app))
store.write("osmStupaAppRate",     round(rate_app, 1))
store.write("osmStupaAppBinomP",   round(p_app_tier, 4))
store.write("osmStupaClusterRatio",     round(cluster_ratio, 2))
store.write("osmStupaClusterApMean",    round(mean_ap_node, 2))
store.write("osmStupaClusterNonApMean", round(mean_nonap_node, 2))
store.write("osmStupaFisherTwoCorpusChi",   round(chi2_two, 3))
store.write("osmStupaFisherTwoCorpusP",     round(p_two_fisher, 4))
store.write("osmStupaFisherThreeCorpusChi", round(chi2_three, 3))
store.write("osmStupaFisherThreeCorpusP",   round(p_three_fisher, 4))
store.write("osmStupaNManMade", int(n_man_made))
store.write("osmStupaNHistoric", int(n_historic))
store.write("osmStupaNBuilding", int(n_building))
store.write("osmStupaNRuins",    int(n_ruins))

# ── Print LaTeX macros ────────────────────────────────────────────────────────

SEP = "=" * 72

def sig(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return r"$^\dagger$"
    return "ns"

print(f"\n{SEP}")
print("  SUMMARY")
print(SEP)
print(f"  OSM stupa corpus (man_made/historic/building/ruins=stupa, dedup 0.01°)")
print(f"  N = {N_osm}  (man_made={n_man_made}, building={n_building}, historic={n_historic}, ruins={n_ruins})")
print(f"  Rayleigh R={R_obs_osm:.4f}, perm-p={rayleigh_p_osm:.4f} {sig(rayleigh_p_osm)}")
print(f"  A+: {n_ap}/{N_osm} ({rate_ap:.1f}%), binomial p={p_ap_tier:.4f} {sig(p_ap_tier)}")
print(f"  A++: {n_app}/{N_osm} ({rate_app:.1f}%), binomial p={p_app_tier:.4f} {sig(p_app_tier)}")
print(f"  Cluster ratio (A-node/non-A-node): {cluster_ratio:.2f}")
print(f"  Fisher OSM+Wikidata: chi2(4)={chi2_two:.3f}, p={p_two_fisher:.4f} {sig(p_two_fisher)}")
print(f"  Fisher OSM+Wiki+UNESCO: chi2(6)={chi2_three:.3f}, p={p_three_fisher:.4f} {sig(p_three_fisher)}")

print(f"\n{SEP}")
print("  LATEX MACROS")
print(SEP)

macros = [
    ("osmStupaTotal",              N_osm),
    ("osmStupaRayleighR",          f"{R_obs_osm:.4f}"),
    ("osmStupaRayleighP",          f"{rayleigh_p_osm:.4f}"),
    ("osmStupaATierCount",         n_a),
    ("osmStupaATierRate",          f"{rate_a:.1f}"),
    ("osmStupaATierBinomP",        f"{p_a_tier:.4f}"),
    ("osmStupaApCount",            n_ap),
    ("osmStupaApRate",             f"{rate_ap:.1f}"),
    ("osmStupaApBinomP",           f"{p_ap_tier:.4f}"),
    ("osmStupaAppCount",           n_app),
    ("osmStupaAppRate",            f"{rate_app:.1f}"),
    ("osmStupaAppBinomP",          f"{p_app_tier:.4f}"),
    ("osmStupaClusterRatio",       f"{cluster_ratio:.2f}"),
    ("osmStupaClusterApMean",      f"{mean_ap_node:.2f}"),
    ("osmStupaClusterNonApMean",   f"{mean_nonap_node:.2f}"),
    ("osmStupaFisherTwoCorpusChi", f"{chi2_two:.3f}"),
    ("osmStupaFisherTwoCorpusP",   f"{p_two_fisher:.4f}"),
    ("osmStupaFisherThreeCorpusChi", f"{chi2_three:.3f}"),
    ("osmStupaFisherThreeCorpusP", f"{p_three_fisher:.4f}"),
    ("osmStupaNManMade",           n_man_made),
    ("osmStupaNHistoric",          n_historic),
    ("osmStupaNBuilding",          n_building),
    ("osmStupaNRuins",             n_ruins),
]

for name, val in macros:
    print(f"  \\newcommand{{\\{name}}}{{{val}}}")
