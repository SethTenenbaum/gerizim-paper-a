"""
cross_corpus_geometric_test.py
==============================
Cross-corpus geometric null test for 3° harmonic alignment.

RATIONALE
---------
The Rayleigh test is inappropriate here.  It tests for *uniform*
concentration of all sites around a single mean direction — it fires on
any geographic clustering, not harmonic alignment in particular.

The correct null is **geometric**:
  Under H₀ a site's longitude is unrelated to the harmonic grid.
  Because harmonics repeat every 3°, a uniformly random longitude
  falls within ±d° of a harmonic with probability p = 2d / 3.

  Tier thresholds (from lib/beru.py / config.json):
    A++  ±0.06°  →  p_null = 0.040   (4 % of the 3° cell)
    A+   ±0.15°  →  p_null = 0.100  (10 %)
    A    ±0.30°  →  p_null = 0.200  (20 %)

Each corpus is tested *independently* with a one-sided binomial test.
If 3 independent corpora all show the same directional excess, the
joint probability is given by Fisher's combination of their p-values:

    χ²_F  =  −2 · Σ ln(p_i)   ~   χ²(2k),  k = number of corpora

Three corpora:
  1. UNESCO dome subset      (keyword-filtered cultural WHS)
  2. Wikidata Q180987 stupas (direct P31/P279 stupa instances)
  3. OSM stupas              (building=stupa / historic=stupa)

Results are written to the ResultsStore and printed as LaTeX macros.
"""

import sys, json, csv
import numpy as np
from pathlib import Path
from scipy.stats import binomtest, chi2

_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(_ROOT))

from data.unesco_corpus import load_corpus
from lib.dome_filter    import is_dome_site_raw as is_dome_site
from lib.beru           import (
    deviation as beru_deviation,
    P_NULL_APP, P_NULL_AP, P_NULL_A,
    TIER_APP_DEG, TIER_APLUS_DEG, TIER_A_DEG,
    GERIZIM,
)
from lib.results_store  import ResultsStore

_CFG  = json.loads((_ROOT / "config.json").read_text())
SEED  = _CFG["simulation"]["random_seed"]
T     = 3.0   # harmonic period (degrees)

SEP = "=" * 72

# ── helper: deviation from nearest harmonic (degrees) ─────────────────────────
def dev_deg(lon: float) -> float:
    """Distance in degrees from the nearest 3° harmonic anchored at Gerizim."""
    arc = (lon - GERIZIM) % T
    return min(arc, T - arc)

def classify(lon: float):
    d = dev_deg(lon)
    if d <= TIER_APP_DEG:   return "A++"
    if d <= TIER_APLUS_DEG: return "A+"
    if d <= TIER_A_DEG:     return "A"
    return "other"

# ── 1. Load the three corpora ─────────────────────────────────────────────────

# UNESCO dome subset
_corpus   = load_corpus()
_cultural = [s for s in _corpus if s.category != "Natural" and s.has_coords]
_dome     = [s for s in _cultural if is_dome_site(s)]
lons_dome = np.array([s.longitude for s in _dome])

# Wikidata Q180987 stupas
STUPA_CSV = _ROOT / "data" / "store" / "unesco" / "wikidata_stupas_q180987.csv"
_wiki_lons = []
with open(STUPA_CSV, newline="", encoding="utf-8") as fh:
    lines = [l for l in fh if not l.startswith("#")]
    for row in csv.DictReader(lines):
        try:
            _wiki_lons.append(float(row["lon"]))
        except (KeyError, ValueError):
            pass
lons_wiki = np.array(_wiki_lons)

# OSM stupas
OSM_CSV = _ROOT / "data" / "store" / "unesco" / "osm_stupas.csv"
_osm_lons = []
with open(OSM_CSV, newline="", encoding="utf-8") as fh:
    lines = [l for l in fh if not l.startswith("#")]
    for row in csv.DictReader(lines):
        try:
            _osm_lons.append(float(row["lon"]))
        except (KeyError, ValueError):
            pass
lons_osm = np.array(_osm_lons)

corpora = [
    ("UNESCO domes",     lons_dome),
    ("Wikidata stupas",  lons_wiki),
    ("OSM stupas",       lons_osm),
]

# ── 2. Per-corpus binomial tests ──────────────────────────────────────────────

def binom_row(lons: np.ndarray, threshold_deg: float, p_null: float):
    """Return (n, hits, p_binom, enrichment_fold)."""
    n    = len(lons)
    hits = int(np.sum([dev_deg(lon) <= threshold_deg for lon in lons]))
    p    = binomtest(hits, n, p_null, alternative="greater").pvalue
    enr  = (hits / n) / p_null if n else 0.0
    return n, hits, p, enr

def sig_stars(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"

def fmt_p_macro(p):
    """Return a LaTeX-safe p-value display string for generated macros."""
    if p < 0.001:
        return r"{<}0.001"
    return round(p, 3)

tiers = [
    ("A++", TIER_APP_DEG,   P_NULL_APP),
    ("A+",  TIER_APLUS_DEG, P_NULL_AP),
    ("A",   TIER_A_DEG,     P_NULL_A),
]

print()
print(SEP)
print("  CROSS-CORPUS GEOMETRIC NULL TEST  —  3° harmonic alignment")
print(f"  Anchor: Gerizim {GERIZIM}°E   |   Period T = {T}°")
print(SEP)

# Store per-corpus p-values for Fisher combination
results = {}   # results[(corpus_name, tier_label)] = (n, hits, p, enr)

for tier_label, thresh, p_null in tiers:
    print(f"\n  Tier {tier_label}  (±{thresh}°,  geometric null p = {p_null:.3f})")
    print(f"  {'Corpus':<22} {'N':>6}  {'Hits':>5}  {'Exp':>6}  {'Fold':>5}  {'p_binom':>9}  sig")
    print(f"  {'-'*68}")
    for name, lons in corpora:
        n, hits, p, enr = binom_row(lons, thresh, p_null)
        exp = n * p_null
        print(f"  {name:<22} {n:>6}  {hits:>5}  {exp:>6.1f}  {enr:>5.2f}x  {p:>9.4f}  {sig_stars(p)}")
        results[(name, tier_label)] = (n, hits, p, enr)

# ── 3. Fisher combination across corpora ─────────────────────────────────────

print()
print(SEP)
print("  FISHER COMBINATION — joint probability across 3 independent corpora")
print(SEP)

for tier_label, thresh, p_null in tiers:
    ps = [results[(name, tier_label)][2] for name, _ in corpora]
    chi2_f = -2.0 * sum(np.log(max(p, 1e-15)) for p in ps)
    p_f    = float(chi2.sf(chi2_f, df=2 * len(ps)))
    print(f"\n  Tier {tier_label}: p-values = {[round(p,4) for p in ps]}")
    print(f"    Fisher χ²({2*len(ps)}) = {chi2_f:.3f},  combined p = {p_f:.5f}  {sig_stars(p_f)}")
    results[("Fisher", tier_label)] = (chi2_f, p_f)

# ── 4. Individual site listing for A++ hits ───────────────────────────────────

print()
print(SEP)
print("  A++ HITS BY CORPUS  (±0.06°, ≤6.7 km from harmonic)")
print(SEP)
for name, lons in corpora:
    hits_idx = [i for i, lon in enumerate(lons) if dev_deg(lon) <= TIER_APP_DEG]
    print(f"\n  {name}  ({len(hits_idx)} A++ hits / {len(lons)} total)")
    for i in hits_idx:
        lon = lons[i]
        d   = dev_deg(lon)
        nearest = GERIZIM + round((lon - GERIZIM) / T) * T
        print(f"    lon={lon:9.4f}°   dev={d:.4f}°   nearest harmonic={nearest:.3f}°")

# ── 5. LaTeX macros ───────────────────────────────────────────────────────────

macros = {}

_corpus_keys = [
    ("UNESCO domes",    "Dome"),
    ("Wikidata stupas", "Wiki"),
    ("OSM stupas",      "Osm"),
]
_tier_keys = [
    ("A++", "App"),
    ("A+",  "Ap"),
    ("A",   "A"),
]

for (cname, ckey), (tlabel, tkey) in [
    (c, t) for c in _corpus_keys for t in _tier_keys
]:
    n, hits, p, enr = results[(cname, tlabel)]
    exp = n * dict(zip(["A++","A+","A"], [P_NULL_APP, P_NULL_AP, P_NULL_A]))[tlabel]
    macros[f"geo{ckey}{tkey}N"]    = n
    macros[f"geo{ckey}{tkey}Hits"] = hits
    macros[f"geo{ckey}{tkey}Exp"]  = round(exp, 1)
    macros[f"geo{ckey}{tkey}Enr"]  = round(enr, 2)
    macros[f"geo{ckey}{tkey}P"]    = round(p, 4)
    macros[f"geo{ckey}{tkey}Sig"]  = sig_stars(p)

for (tlabel, tkey) in _tier_keys:
    chi2_f, p_f = results[("Fisher", tlabel)]
    macros[f"geoFisher{tkey}Chi"]  = round(chi2_f, 3)
    macros[f"geoFisher{tkey}P"]    = fmt_p_macro(p_f)
    macros[f"geoFisher{tkey}Sig"]  = sig_stars(p_f)

# ── Non-circular two-corpus Fisher (Wiki + OSM only) ─────────────────────────
_nc_corpora = ["Wikidata stupas", "OSM stupas"]
for (tlabel, tkey) in _tier_keys:
    ps_nc = [results[(cname, tlabel)][2] for cname in _nc_corpora]
    chi2_nc = -2.0 * sum(np.log(max(p, 1e-15)) for p in ps_nc)
    p_nc    = float(chi2.sf(chi2_nc, df=2 * len(ps_nc)))
    macros[f"geoFisherNc{tkey}Chi"] = round(chi2_nc, 3)
    macros[f"geoFisherNc{tkey}P"]   = fmt_p_macro(p_nc)
    macros[f"geoFisherNc{tkey}Sig"] = sig_stars(p_nc)
    print(f"  Non-circ Fisher Tier {tlabel}: χ²(4)={chi2_nc:.3f}, p={p_nc:.5f} {sig_stars(p_nc)}")

print()
print(SEP)
print("  LATEX MACROS")
print(SEP)
for k, v in macros.items():
    print(f"  \\newcommand{{\\{k}}}{{{v}}}")

ResultsStore().write_many(macros)
print()
print("  Results written to data/store/results.json")
