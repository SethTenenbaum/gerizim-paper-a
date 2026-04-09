"""
unesco_buddhist_heritage_test.py
================================
Beru-grid enrichment test for UNESCO World Heritage sites with Buddhist heritage.

MOTIVATION
──────────
The UNESCO World Heritage List contains ~61 Cultural/Mixed properties whose
nomination text references Buddhist heritage (temples, monasteries, stupas,
grottoes, etc.).  This is a broader category than the stupa-specific test
(Test 2) — it includes all UNESCO sites that UNESCO itself classifies as
relevant to Buddhist cultural heritage, regardless of architectural form.

This test asks: do UNESCO Buddhist heritage sites show non-uniform
distribution with respect to beru harmonics from Gerizim?

Unlike Test 2 (author-compiled stupa database), this test uses the full
UNESCO WHC XML corpus, filtered only by keyword.  The selection criterion
is purely textual — "does the UNESCO short_description or site name
reference Buddhist heritage?" — with no reference to longitude or beru.

METHODOLOGY
───────────
1. Parse the full UNESCO WHC XML (Cultural + Mixed sites with coords).
2. Select sites whose 'site' name or 'short_description' contains any of
   a pre-defined set of Buddhist heritage keywords.
3. Compute beru deviation from Gerizim for each site.
4. Test for Tier-A+ and Tier-A enrichment (binomial, one-tailed).
5. Compare to the non-Buddhist UNESCO corpus (cluster asymmetry).

KEYWORDS (frozen before any beru calculation):
  buddhist, buddha, stupa, dagoba, vihara, sangha, dharma,
  monastery (in Buddhist-country context), temple (in Buddhist context),
  borobudur, lumbini, ajanta, ellora, mogao, yungang, longmen,
  bagan, pyu, anuradhapura, kandy, dambulla, sanchi, nalanda,
  bodh gaya, mahabodhi, horyu, nara, kyoto, nikko, hiraizumi,
  sansa, bulguksa, seokguram, gyeongju, baekje,
  sukhothai, ayutthaya, si thep, luang prabang, potala,
  chengde, wutai, dengfeng, taxila, takht-i-bahi, paharpur,
  sambor prei kuk, angkor, prambanan, yogyakarta

ANCHOR:  Mount Gerizim 35.274°E
BERU:    30° arc
"""

import re
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
from scipy.stats import mannwhitneyu, binomtest, fisher_exact

# ── Use common corpus library ──────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_A_MAX, TIER_B_MAX, P_NULL_AP, P_NULL_A,
    HARMONIC_STEP, CONFIG,
    deviation as _beru_dev, tier_label, is_aplus, is_a_or_better,
    load_keywords,
)
from lib.stats import significance_label as sig

HALF_STEP = HARMONIC_STEP / 2.0   # 0.05 beru

# buddhist_heritage in config is {"note":..., "patterns":[...]}
_bud_cfg = CONFIG["keywords"]["buddhist_heritage"]
BUDDHIST_KEYWORDS = (
    _bud_cfg["patterns"] if isinstance(_bud_cfg, dict) and "patterns" in _bud_cfg
    else load_keywords("buddhist_heritage")
)

COMPILED_KEYWORDS = [re.compile(kw, re.IGNORECASE) for kw in BUDDHIST_KEYWORDS]


# ── Anti-node deviation helper ───────────────────────────────────────────────
def anti_node_deviation(lon: float,
                        anchor: float = GERIZIM,
                        beru: float = BERU,
                        step: float = HARMONIC_STEP) -> float:
    """
    Distance (in beru) from the nearest anti-node (harmonic midpoint).
    Anti-nodes sit at (n + 0.5)*step for integer n.
    Returns a value in [0, step/2].
    """
    arc      = abs(lon - anchor)
    arc      = min(arc, 360.0 - arc)
    beru_val = arc / beru
    shifted  = beru_val - HALF_STEP
    nearest_anti = round(shifted / step) * step + HALF_STEP
    return abs(beru_val - nearest_anti)


# ── Helpers ─────────────────────────────────────────────────────────────────
def beru_deviation(lon: float):
    arc      = abs(lon - GERIZIM)
    beru_val = arc / BERU
    nearest  = round(beru_val * 10) / 10
    dev      = _beru_dev(lon)
    return arc, beru_val, nearest, dev


def matches_buddhist(name: str, desc: str, extended: str = "") -> list:
    """Return list of keyword patterns that match in name+desc+extended."""
    text = f"{name} {desc} {extended}"
    hits = []
    for i, pat in enumerate(COMPILED_KEYWORDS):
        if pat.search(text):
            hits.append(BUDDHIST_KEYWORDS[i])
    return hits


# ── Load UNESCO sites ───────────────────────────────────────────────────────
def load_unesco():
    corpus = load_corpus()

    all_sites = []
    buddhist_sites = []

    for site_obj in corpus:
        name = site_obj.site
        cat  = site_obj.category
        desc = site_obj.short_description
        extended = site_obj.extended_description

        if cat == "Natural":
            continue

        if not site_obj.has_coords:
            continue

        lat = site_obj.latitude
        lon = site_obj.longitude

        arc, beru_val, nearest, dev = beru_deviation(lon)
        anti_dev = anti_node_deviation(lon)
        site = {
            "name": name,
            "lat":  lat,
            "lon":  lon,
            "cat":  cat,
            "desc": desc[:200],
            "arc":  arc,
            "beru": beru_val,
            "nearest": nearest,
            "dev":     dev,
            "anti_dev": anti_dev,
            "dev_km":  dev * BERU * 111.0,
            "anti_dev_km": anti_dev * BERU * 111.0,
            "tier":    tier_label(dev),
        }
        all_sites.append(site)

        kw_hits = matches_buddhist(name, desc, extended)
        if kw_hits:
            site["keywords"] = kw_hits
            buddhist_sites.append(site)

    return all_sites, buddhist_sites


# ── Main ────────────────────────────────────────────────────────────────────
all_sites, buddhist_sites = load_unesco()

N_all = len(all_sites)
N_bud = len(buddhist_sites)

print("=" * 95)
print("  UNESCO BUDDHIST HERITAGE SITES — BERU-GRID TEST")
print(f"  Anchor: Gerizim {GERIZIM}°E  |  BERU={BERU}°")
print(f"  Source: UNESCO WHC XML, keyword-filtered for Buddhist heritage")
print("=" * 95)
print()

# ── List all Buddhist sites ──────────────────────────────────────────────────
n_Ap = sum(1 for s in buddhist_sites if is_aplus(s["tier"]))
n_A  = sum(1 for s in buddhist_sites if is_a_or_better(s["tier"]))
n_B  = sum(1 for s in buddhist_sites if s["tier"] == "B")
n_C  = sum(1 for s in buddhist_sites if s["tier"] == "C")

print(f"  Total Buddhist heritage sites:  N = {N_bud}")
print(f"  Tier-A+ (≤6.7 km from harmonic):  {n_Ap}  ({100*n_Ap/N_bud:.1f}%)  [null: 4%]")
print(f"  Tier-A  (≤33 km):                  {n_A}  ({100*n_A/N_bud:.1f}%)")
print(f"  Tier-B  (≤166.5 km):              {n_B}  ({100*n_B/N_bud:.1f}%)")
print(f"  Tier-C  (>166.5 km):              {n_C}  ({100*n_C/N_bud:.1f}%)")
print()

# ── Binomial tests ───────────────────────────────────────────────────────────
bt_Ap = binomtest(n_Ap, N_bud, P_NULL_AP, alternative="greater")
bt_A  = binomtest(n_A, N_bud, P_NULL_A, alternative="greater")

print(f"  Tier-A+ binomial (H₁: rate > {P_NULL_AP:.0%}):  {n_Ap}/{N_bud} = {100*n_Ap/N_bud:.1f}%")
print(f"    p = {bt_Ap.pvalue:.4f}  {sig(bt_Ap.pvalue)}")
print(f"    enrichment = {(n_Ap/N_bud)/P_NULL_AP:.2f}×")
print()
print(f"  Tier-A binomial (H₁: rate > {P_NULL_A:.0%}):  {n_A}/{N_bud} = {100*n_A/N_bud:.1f}%")
print(f"    p = {bt_A.pvalue:.4f}  {sig(bt_A.pvalue)}")
print(f"    enrichment = {(n_A/N_bud)/P_NULL_A:.2f}×")
print()

# ── Per-site table ──────────────────────────────────────────────────────────
print("=" * 95)
print("  ALL BUDDHIST HERITAGE SITES (sorted by deviation)")
print("=" * 95)
print()
print(f"  {'#':>3}  {'Site':<60s}  {'lon':>9}  {'beru':>6}  {'dev':>7}  {'km':>6}  {'Tier':<4}  {'Keywords'}")
print(f"  {'─'*3}  {'─'*60}  {'─'*9}  {'─'*6}  {'─'*7}  {'─'*6}  {'─'*4}  {'─'*30}")

for i, s in enumerate(sorted(buddhist_sites, key=lambda x: x["dev"]), 1):
    kw_str = ", ".join(s.get("keywords", []))[:40]
    flag = " ◀◀" if is_aplus(s["tier"]) else (" ◀" if s["tier"] == "A" else "")
    print(f"  {i:>3}  {s['name'][:60]:<60s}  {s['lon']:>9.4f}  {s['beru']:>6.4f}  "
          f"{s['dev']:>7.5f}  {s['dev_km']:>6.1f}  {s['tier']:<4s}{flag}  {kw_str}")

print()

# ── Cluster asymmetry: compare Buddhist A+ sites to Buddhist non-A+ ────────
print("=" * 95)
print("  CLUSTER ASYMMETRY WITHIN BUDDHIST SITES")
print("=" * 95)
print()

CLUSTER_RADIUS_DEG = 0.5
lons = np.array([s["lon"] for s in buddhist_sites])
n_bs = len(buddhist_sites)

for i in range(n_bs):
    dlon = np.abs(lons - lons[i])
    dlon = np.minimum(dlon, 360 - dlon)
    buddhist_sites[i]["cluster_size"] = int(np.sum(dlon <= CLUSTER_RADIUS_DEG)) - 1

ap_bud = [s for s in buddhist_sites if is_aplus(s["tier"])]
nonap_bud = [s for s in buddhist_sites if not is_aplus(s["tier"])]

if ap_bud and nonap_bud:
    mean_ap = np.mean([s["cluster_size"] for s in ap_bud])
    mean_nonap = np.mean([s["cluster_size"] for s in nonap_bud])
    print(f"  A+ Buddhist sites ({len(ap_bud)}): mean cluster size = {mean_ap:.1f}")
    print(f"  Non-A+ Buddhist sites ({len(nonap_bud)}): mean cluster size = {mean_nonap:.1f}")

    if len(ap_bud) >= 2:
        U, mw_p = mannwhitneyu(
            [s["cluster_size"] for s in ap_bud],
            [s["cluster_size"] for s in nonap_bud],
            alternative="greater"
        )
        print(f"  Mann-Whitney U (A+ > non-A+): p = {mw_p:.4f}  {sig(mw_p)}")
    print()

# ── By geographic band ──────────────────────────────────────────────────────
print("=" * 95)
print("  BY GEOGRAPHIC BAND (Mauryan lineage vs. East Asian)")
print("=" * 95)
print()

# Define bands
MAURYAN_LON_RANGE = (60, 120)   # Bactria → Indonesia
EAST_ASIAN_RANGE  = (100, 180)  # China → Japan (overlaps slightly)

# More precise: South/SE Asian band vs East Asian band
def classify_region(s):
    lon = s["lon"]
    lat = s["lat"]
    name = s["name"].lower()
    # Japan
    if 125 < lon < 146 and 30 < lat < 46:
        return "East Asian"
    # Korea
    if 124 < lon < 132 and 33 < lat < 43:
        return "East Asian"
    # China proper (excluding western Silk Road)
    if 100 < lon < 125 and 20 < lat < 45:
        # Silk Road sites in western China
        if lon < 105 and lat > 35:
            return "Silk Road/Central Asian"
        if "mogao" in name or "yungang" in name or "longmen" in name:
            return "East Asian"
        if "silk road" in name or "xanadu" in name:
            return "Silk Road/Central Asian"
        return "East Asian"
    # South Asia
    if 65 < lon < 100 and 0 < lat < 40:
        return "South/SE Asian (Mauryan lineage)"
    # SE Asia (Myanmar, Thailand, Cambodia, Indonesia)
    if 93 < lon < 115 and -10 < lat < 25:
        return "South/SE Asian (Mauryan lineage)"
    # Central Asia
    if 55 < lon < 80 and 30 < lat < 50:
        return "Silk Road/Central Asian"
    return "Other"

for s in buddhist_sites:
    s["region"] = classify_region(s)

regions = defaultdict(list)
for s in buddhist_sites:
    regions[s["region"]].append(s)

for region in ["South/SE Asian (Mauryan lineage)", "East Asian",
               "Silk Road/Central Asian", "Other"]:
    rsites = regions.get(region, [])
    if not rsites:
        continue
    n_r = len(rsites)
    n_r_Ap = sum(1 for s in rsites if is_aplus(s["tier"]))
    n_r_A = sum(1 for s in rsites if is_a_or_better(s["tier"]))
    print(f"  {region} (N={n_r})")
    print(f"    A+ = {n_r_Ap}/{n_r} = {100*n_r_Ap/n_r:.1f}%  "
          f"(null: 4%,  enrichment: {(n_r_Ap/n_r)/P_NULL_AP:.2f}×)")
    if n_r >= 3:
        bt = binomtest(n_r_Ap, n_r, P_NULL_AP, alternative="greater")
        print(f"    binomial p = {bt.pvalue:.4f}  {sig(bt.pvalue)}")
    print(f"    A  = {n_r_A}/{n_r} = {100*n_r_A/n_r:.1f}%")
    print()

# ── Harmonic density: do Buddhist-occupied harmonics have more sites? ───────
print("=" * 95)
print("  HARMONIC DENSITY (same test as cluster_asymmetry_test.py)")
print("  Do harmonics with A+ Buddhist sites host more total UNESCO sites?")
print("=" * 95)
print()

# Use all UNESCO sites for density measurement
all_at_harmonic = defaultdict(list)
for s in all_sites:
    if s["dev"] <= TIER_B_MAX:
        all_at_harmonic[s["nearest"]].append(s)

bud_ap_harmonics = set(s["nearest"] for s in buddhist_sites if is_aplus(s["tier"]))
all_harmonics = set(s["nearest"] for s in all_sites if s["dev"] <= TIER_B_MAX)

density_ap = [len(all_at_harmonic[h]) for h in all_harmonics if h in bud_ap_harmonics]
density_nonap = [len(all_at_harmonic[h]) for h in all_harmonics if h not in bud_ap_harmonics]

if density_ap and density_nonap:
    print(f"  Harmonics with ≥1 Buddhist A+ site: {len(density_ap)}")
    print(f"    Mean total UNESCO sites per harmonic: {np.mean(density_ap):.1f}")
    print(f"  Harmonics with 0 Buddhist A+ sites: {len(density_nonap)}")
    print(f"    Mean total UNESCO sites per harmonic: {np.mean(density_nonap):.1f}")

    if len(density_ap) >= 2:
        U, p = mannwhitneyu(density_ap, density_nonap, alternative="greater")
        print(f"  Mann-Whitney U: p = {p:.4f}  {sig(p)}")
        print(f"  Ratio: {np.mean(density_ap)/np.mean(density_nonap):.2f}×")
    print()

# ── Anti-node analysis (Mauryan lineage sub-group) ─────────────────────────
print("=" * 95)
print("  ANTI-NODE ANALYSIS — BUDDHIST HERITAGE SITES")
print("  (proximity to harmonic MIDPOINTS, not harmonics)")
print("=" * 95)
print()

# A site is anti-node-eligible if it is closer to a midpoint than a harmonic,
# i.e. its fractional phase puts it in the anti-node half of the grid cell.
# Tier-A+ at the anti-node = anti_dev <= TIER_APLUS (0.002 beru = 6.7 km)
# Tier-A  at the anti-node = anti_dev <= TIER_A_MAX  (0.010 beru = 33.3 km)

def is_antinode_eligible(s):
    """True when closer to a midpoint than to a harmonic."""
    return s["anti_dev"] < s["dev"]

def is_antinode_aplus(s):
    return s["anti_dev"] <= TIER_APLUS

def is_antinode_a(s):
    return s["anti_dev"] <= TIER_A_MAX

# Full Buddhist corpus
anti_eligible = [s for s in buddhist_sites if is_antinode_eligible(s)]
n_anti = len(anti_eligible)
n_anti_Ap = sum(1 for s in anti_eligible if is_antinode_aplus(s))
n_anti_A  = sum(1 for s in anti_eligible if is_antinode_a(s))

print(f"  Anti-node-eligible Buddhist sites (closer to midpoint than harmonic): {n_anti}")
if n_anti > 0:
    bt_anti_Ap = binomtest(n_anti_Ap, n_anti, P_NULL_AP, alternative="greater")
    bt_anti_A  = binomtest(n_anti_A,  n_anti, P_NULL_A,  alternative="greater")
    print(f"  Tier-A+ (anti_dev ≤ 0.002 beru): {n_anti_Ap}/{n_anti} = "
          f"{100*n_anti_Ap/n_anti:.1f}%  null 4%  p={bt_anti_Ap.pvalue:.4f}")
    print(f"  Tier-A  (anti_dev ≤ 0.010 beru): {n_anti_A}/{n_anti}  = "
          f"{100*n_anti_A/n_anti:.1f}%  null 20% p={bt_anti_A.pvalue:.4f}")
print()

# Node vs anti-node rate comparison (Fisher exact)
node_eligible = [s for s in buddhist_sites if not is_antinode_eligible(s)]
n_node = len(node_eligible)
n_node_Ap = sum(1 for s in node_eligible if is_aplus(s["tier"]))
if n_node > 0 and n_anti > 0:
    ct = [[n_node_Ap, n_node - n_node_Ap],
          [n_anti_Ap, n_anti - n_anti_Ap]]
    or_val, fe_p = fisher_exact(ct, alternative="greater")
    print(f"  Node-side A+ rate:      {n_node_Ap}/{n_node} = {100*n_node_Ap/n_node:.1f}%")
    print(f"  Anti-node-side A+ rate: {n_anti_Ap}/{n_anti} = {100*n_anti_Ap/n_anti:.1f}%")
    print(f"  Fisher exact (node > anti-node): OR={or_val:.2f}  p={fe_p:.4f}")
    print()

# Regional breakdown of anti-node signal
print("  Anti-node A+ sites by region:")
for region in ["South/SE Asian (Mauryan lineage)", "East Asian",
               "Silk Road/Central Asian", "Other"]:
    rsites = [s for s in anti_eligible if s.get("region") == region]
    if not rsites:
        continue
    n_r = len(rsites)
    n_r_Ap = sum(1 for s in rsites if is_antinode_aplus(s))
    print(f"    {region} (N={n_r}): anti-node A+ = {n_r_Ap}/{n_r} = "
          f"{100*n_r_Ap/n_r:.1f}%")
print()

# List all anti-node A+ Buddhist sites
anti_ap_sites = sorted(
    [s for s in buddhist_sites if is_antinode_aplus(s)],
    key=lambda x: x["anti_dev"]
)
print(f"  Buddhist sites at Tier-A+ anti-node precision (anti_dev ≤ 0.002 beru):")
if anti_ap_sites:
    print(f"  {'Site':<60s}  {'lon':>9}  {'anti_dev':>9}  {'km':>6}  {'region'}")
    print(f"  {'─'*60}  {'─'*9}  {'─'*9}  {'─'*6}  {'─'*25}")
    for s in anti_ap_sites:
        print(f"  {s['name'][:60]:<60s}  {s['lon']:>9.4f}  "
              f"{s['anti_dev']:>9.5f}  {s['anti_dev_km']:>6.1f}  "
              f"{s.get('region', 'Unknown')}")
else:
    print("  (none)")
print()

# Mauryan-specific anti-node test (pre-specified sub-group)
mauryan = [s for s in buddhist_sites if s.get("region") == "South/SE Asian (Mauryan lineage)"]
m_anti  = [s for s in mauryan if is_antinode_eligible(s)]
m_node  = [s for s in mauryan if not is_antinode_eligible(s)]
n_m_anti_Ap = sum(1 for s in m_anti if is_antinode_aplus(s))
n_m_node_Ap = sum(1 for s in m_node if is_aplus(s["tier"]))

print(f"  Mauryan-lineage sub-group (N={len(mauryan)}):")
print(f"    Node-side:      N={len(m_node)}, A+ = {n_m_node_Ap}/{len(m_node)} = "
      f"{100*n_m_node_Ap/len(m_node):.1f}%" if m_node else "    Node-side: N=0")
print(f"    Anti-node-side: N={len(m_anti)}, A+ = {n_m_anti_Ap}/{len(m_anti)} = "
      f"{100*n_m_anti_Ap/len(m_anti):.1f}%" if m_anti else "    Anti-node-side: N=0")
if m_node and m_anti:
    ct_m = [[n_m_node_Ap, len(m_node) - n_m_node_Ap],
            [n_m_anti_Ap, len(m_anti) - n_m_anti_Ap]]
    or_m, fe_m = fisher_exact(ct_m, alternative="greater")
    print(f"    Fisher exact (Mauryan node > anti-node): OR={or_m:.2f}  p={fe_m:.4f}")
print()

# ── Summary ──────────────────────────────────────────────────────────────────
print("=" * 95)
print("  SUMMARY")
print("=" * 95)
print(f"""
  UNESCO Buddhist Heritage sites (keyword-filtered):  N = {N_bud}
  Tier-A+ rate: {n_Ap}/{N_bud} = {100*n_Ap/N_bud:.1f}%  (null: 4%)
  Tier-A+ binomial p = {bt_Ap.pvalue:.4f}  {sig(bt_Ap.pvalue)}
  Tier-A  rate: {n_A}/{N_bud} = {100*n_A/N_bud:.1f}%  (null: 20%)
  Tier-A  binomial p = {bt_A.pvalue:.4f}  {sig(bt_A.pvalue)}

  INTERPRETATION:
  This test uses the full UNESCO WHC corpus, filtered by Buddhist heritage
  keywords applied to site names and short descriptions.  The keywords
  correspond to the UNESCO Buddhist Heritage thematic scope — the same
  sites that appear in UNESCO's own "Buddhist" search filter.

  The test is broader than Test 2 (stupas only) and includes temples,
  grottoes, monasteries, and other forms.  Under the beru-grid hypothesis,
  the signal should be STRONGEST for stupa-form monuments (Test 2) and
  diluted in a broader Buddhist corpus that includes many derivative and
  non-stupa sites.

  The key finding is not whether the overall binomial is significant,
  but whether the PATTERN matches the prediction: Mauryan-lineage sites
  enriched, East Asian sites at null, and the signal concentrated at
  Tier-A+ (not Tier-A), indicating precise rather than approximate
  alignment.
""")
