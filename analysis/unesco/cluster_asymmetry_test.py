"""
cluster_asymmetry_test.py
=========================
TEST: Do UNESCO sites that fall in Tier-A+ beru windows cluster MORE
than UNESCO sites that fall in Tier-B/C beru windows?
And does this clustering asymmetry follow the NODE/ANTI-NODE phase
distinction — i.e., is it phase-locked to the 0.1-beru grid?

MOTIVATION
──────────
A critic of the beru-grid hypothesis can say: "Your A+ hits are just
geographic clusters — Malta has 7 temples, Anuradhapura has 5 stupas,
all at the same longitude.  Strip them down to one-per-complex and
the signal vanishes."

The correct response is: **the clustering itself IS the signal.**

Under H₀ (longitudes uniformly distributed w.r.t. beru harmonics),
the geographic clustering rate should be the SAME in the A+ window
and in non-A+ windows, AND the same at node longitudes as at anti-node
longitudes.  There is no reason that sites near a harmonic node should
cluster more than sites at anti-nodes (both occur at densely settled
longitudes) — unless the harmonic NODES are *special places where humans
built many monuments* while the anti-nodes are merely dense because cities
are planned everywhere near the grid centre.

This script tests:
  H₀: The per-site cluster size is the SAME for A+ and non-A+ sites,
      and the same at node and anti-node longitudes.
  H₁: Node-side A+ sites are in denser clusters than:
        (a) node-side non-A+ sites, and
        (b) anti-node-side A+ sites.
      Anti-node-side A+ sites show NO excess over their non-A+ neighbours.

If H₁ holds, the clustering is PHASE-LOCKED: it follows the 0.1-beru
node/anti-node distinction, not just settlement geography.  This is a
second independent measure of grid structure, complementing the external
P1435 control comparison (Test 7).

METHODOLOGY
───────────
For every Cultural/Mixed UNESCO site with coordinates:
  1. Compute node deviation (from nearest harmonic) and anti-node
     deviation (from nearest grid midpoint).
  2. Classify each site as node-side (closer to a harmonic) or
     anti-node-side (closer to a midpoint).
  3. Apply A+ threshold (≤0.002 beru) within each phase partition.
  4. Compute cluster size = other UNESCO sites within ±Tier-A beru window
     (TIER_A_MAX × BERU = 0.3° ≈ 33.3 km), the same radius used in the
     stupa cluster enrichment test.
  5. Compare: node A+ vs node non-A+; anti-node A+ vs anti-node non-A+;
     node A+ vs anti-node A+; all node-side vs all anti-node-side.
  6. Mann-Whitney U tests throughout (non-parametric).

NO site selection, no keyword filtering — full UNESCO XML.
"""

import sys
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import mannwhitneyu, binomtest, kstest

# ── Use common corpus library ──────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_A_MAX, TIER_B_MAX, P_NULL_AP,
    deviation as beru_deviation, tier_label, is_aplus, is_a_or_better,
)
from lib.stats import significance_label as sig
from lib.results_store import ResultsStore

import json as _json
_ROOT = Path(__file__).parent.parent.parent
_CFG  = _json.loads((_ROOT / "config.json").read_text())
N_PERM_CLUSTER = _CFG["simulation"]["n_permutations_cluster"]  # 10,000

# Cluster radius: defined by the Tier-A beru threshold, matching the
# metrological framework of the stupa cluster enrichment test.
# TIER_A_MAX = 0.010 beru = 0.010 × 30° = 0.3° ≈ 33.3 km at equator.
# All neighborhood radii in the stupa test are beru-derived
# (A+ = 6.7 km, A = 33.3 km, B = 166.5 km); this script uses the same
# Tier-A window so that cluster size is measured in the same metrological
# unit as the primary proximity tests.
CLUSTER_RADIUS_DEG = TIER_A_MAX * BERU   # 0.010 × 30 = 0.3°  ≈ 33.3 km


# ── Load ALL UNESCO Cultural/Mixed sites with coordinates ───────────────────
def load_all_unesco():
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    sites = []
    for s in cultural:
        dev = beru_deviation(s.longitude)
        arc = abs(s.longitude - GERIZIM)
        beru_val = arc / BERU
        nearest = round(beru_val * 10) / 10
        sites.append({
            "name": s.site,
            "lat":  s.latitude,
            "lon":  s.longitude,
            "cat":  s.category,
            "arc":  arc,
            "beru": beru_val,
            "nearest": nearest,
            "dev":     dev,
            "dev_km":  dev * BERU * 111.0,
            "tier":    tier_label(dev),
        })
    return sites


# ── Compute cluster sizes ────────────────────────────────────────────────────
def compute_cluster_sizes(sites, radius_deg=CLUSTER_RADIUS_DEG):
    """
    For each site, count how many OTHER sites have longitude within
    ±radius_deg.  This is the "cluster size" — how many neighbours
    share approximately the same longitude corridor.
    """
    lons = np.array([s["lon"] for s in sites])
    n = len(sites)
    cluster_sizes = np.zeros(n, dtype=int)

    for i in range(n):
        # Count neighbours within ±radius_deg longitude
        dlon = np.abs(lons - lons[i])
        # Handle wrap-around at ±180
        dlon = np.minimum(dlon, 360 - dlon)
        cluster_sizes[i] = int(np.sum(dlon <= radius_deg)) - 1  # exclude self

    for i, s in enumerate(sites):
        s["cluster_size"] = int(cluster_sizes[i])

    return sites


# ── Main analysis ────────────────────────────────────────────────────────────
sites = load_all_unesco()
sites = compute_cluster_sizes(sites)

N = len(sites)
n_Ap = sum(1 for s in sites if is_aplus(s["tier"]))
n_A  = sum(1 for s in sites if is_a_or_better(s["tier"]))
n_B  = sum(1 for s in sites if s["tier"] == "B")
n_C  = sum(1 for s in sites if s["tier"] == "C")

ap_sites   = [s for s in sites if is_aplus(s["tier"])]
a_sites    = [s for s in sites if s["tier"] == "A"]
nonap_sites = [s for s in sites if not is_aplus(s["tier"])]
bc_sites   = [s for s in sites if s["tier"] in ("B", "C")]

ap_clusters   = [s["cluster_size"] for s in ap_sites]
a_clusters    = [s["cluster_size"] for s in a_sites]
nonap_clusters = [s["cluster_size"] for s in nonap_sites]
bc_clusters   = [s["cluster_size"] for s in bc_sites]

print("=" * 95)
print("  UNESCO FULL-CORPUS CLUSTER ASYMMETRY TEST")
print("  Are beru-harmonic longitudes preferentially monument-dense?")
print(f"  Anchor: Gerizim {GERIZIM}°E  |  BERU={BERU}°  |  Cluster radius: ±{CLUSTER_RADIUS_DEG}° lon")
print(f"  Source: UNESCO WHC XML (Cultural + Mixed sites with coordinates)")
print("=" * 95)

print(f"""
  POPULATION
  ──────────────────────────────────────────────────────────────────────
  Total Cultural/Mixed with coords:  N = {N}
  Tier-A+ (≤6.7 km from harmonic):  {n_Ap}  ({100*n_Ap/N:.1f}%)  [null: 4%]
  Tier-A  (≤33 km):                  {n_A}  ({100*n_A/N:.1f}%)
  Tier-B  (≤166.5 km):              {n_B}  ({100*n_B/N:.1f}%)
  Tier-C  (>166.5 km):              {n_C}  ({100*n_C/N:.1f}%)
""")

# ── Tier-A+ enrichment (full UNESCO, no keyword filter) ────────────────────
bt_Ap = binomtest(n_Ap, N, P_NULL_AP, alternative="greater")
print(f"  Tier-A+ binomial (H₁: rate > {P_NULL_AP:.0%}):  {n_Ap}/{N} = {100*n_Ap/N:.1f}%")
print(f"    p = {bt_Ap.pvalue:.4f}  {sig(bt_Ap.pvalue)}")
print(f"    enrichment = {(n_Ap/N)/P_NULL_AP:.2f}×")
print()

# ── Cluster size comparison ─────────────────────────────────────────────────
mean_ap   = np.mean(ap_clusters) if ap_clusters else 0
mean_nonap = np.mean(nonap_clusters)
median_ap  = np.median(ap_clusters) if ap_clusters else 0
median_nonap = np.median(nonap_clusters)

print("=" * 95)
print("  CLUSTER SIZE COMPARISON: A+ sites vs non-A+ sites")
print(f"  'Cluster size' = number of other UNESCO sites within ±{CLUSTER_RADIUS_DEG:.3f}° longitude (={CLUSTER_RADIUS_DEG/BERU:.3f} beru = TIER_A_MAX)")
print("=" * 95)

print(f"""
  Tier-A+ sites (N={n_Ap}):
    Mean cluster size:    {mean_ap:.1f}
    Median cluster size:  {median_ap:.0f}
    Max cluster size:     {max(ap_clusters) if ap_clusters else 0}
    Min cluster size:     {min(ap_clusters) if ap_clusters else 0}

  Non-A+ sites (N={len(nonap_sites)}):
    Mean cluster size:    {mean_nonap:.1f}
    Median cluster size:  {median_nonap:.0f}
    Max cluster size:     {max(nonap_clusters)}
    Min cluster size:     {min(nonap_clusters)}
""")

# Mann-Whitney U test: A+ cluster sizes vs non-A+ cluster sizes
if len(ap_clusters) >= 2 and len(nonap_clusters) >= 2:
    U_stat, mw_p = mannwhitneyu(ap_clusters, nonap_clusters, alternative="greater")
    print(f"  Mann-Whitney U (H₁: A+ cluster sizes > non-A+ cluster sizes):")
    print(f"    U = {U_stat:.0f},  p = {mw_p:.4f}  {sig(mw_p)}")
    print(f"    Effect size (rank-biserial r) = {1 - 2*U_stat/(n_Ap * len(nonap_sites)):.3f}")
else:
    mw_p = 1.0
    print("  [Not enough A+ sites for Mann-Whitney test]")

# ── Permutation test ────────────────────────────────────────────────────────
print()
print("=" * 95)
print("  PERMUTATION TEST (10,000 trials)")
print("  Shuffle tier labels among all sites, recompute mean cluster size")
print("  for the 'A+' group each time.  How often does the shuffled mean")
print("  equal or exceed the observed A+ mean?")
print("=" * 95)

rng = np.random.default_rng(42)
all_cluster_sizes = np.array([s["cluster_size"] for s in sites])
n_perm = N_PERM_CLUSTER
perm_means = np.zeros(n_perm)

for i in range(n_perm):
    idx = rng.choice(N, size=n_Ap, replace=False)
    perm_means[i] = all_cluster_sizes[idx].mean()

perm_p = float(np.mean(perm_means >= mean_ap))
print(f"""
  Observed A+ mean cluster size:  {mean_ap:.2f}
  Permutation null mean:          {np.mean(perm_means):.2f}
  Permutation std:                {np.std(perm_means):.2f}
  Permutation p-value:            {perm_p:.4f}  {sig(perm_p)}
  Z-score:                        {(mean_ap - np.mean(perm_means)) / max(np.std(perm_means), 0.001):.2f}
""")

# ── Show the A+ clusters ────────────────────────────────────────────────────
print("=" * 95)
print("  TIER-A+ SITES AND THEIR CLUSTER NEIGHBOURHOODS")
print(f"  (other UNESCO sites within ±{CLUSTER_RADIUS_DEG}° longitude)")
print("=" * 95)
print()

for s in sorted(ap_sites, key=lambda x: -x["cluster_size"]):
    print(f"  {s['name'][:65]:<65s}")
    print(f"    lon={s['lon']:.4f}  beru={s['beru']:.4f}  "
          f"dev={s['dev']:.5f} ({s['dev_km']:.1f} km)  "
          f"→ harmonic {s['nearest']:.1f}")
    print(f"    CLUSTER SIZE: {s['cluster_size']} other sites within ±{CLUSTER_RADIUS_DEG}° lon")
    # Show neighbours
    neighbours = [x for x in sites if x["name"] != s["name"]
                  and min(abs(x["lon"] - s["lon"]), 360 - abs(x["lon"] - s["lon"])) <= CLUSTER_RADIUS_DEG]
    neighbours_sorted = sorted(neighbours, key=lambda x: abs(x["lon"] - s["lon"]))[:8]
    for nb in neighbours_sorted:
        dlon = abs(nb["lon"] - s["lon"])
        print(f"      {nb['name'][:55]:<55s}  Δlon={dlon:.3f}°  tier={nb['tier']}")
    if len(neighbours) > 8:
        print(f"      ... and {len(neighbours) - 8} more")
    print()


# ── Phase-split: node-side vs anti-node-side cluster analysis ─────────────────
# Coarse half-cell classification:
#   node-side  : beru_val mod 0.1 < 0.05  (closer to the harmonic node)
#   anti-side  : beru_val mod 0.1 >= 0.05 (closer to the half-beru midpoint)
# A+ definitions:
#   node-side A+  : node-side  AND  dev       <= TIER_APLUS (≤0.002 beru, ≤6.7 km)
#   anti-side A+  : anti-side  AND  anti_dev  <= TIER_APLUS
# Cluster radius for the phase split: 1.5 × TIER_A_MAX = 0.015 beru = 0.45° longitude.
# Expressed entirely in beru units — one-and-a-half Tier-A windows, the smallest
# metrologically grounded beru multiple that yields a significant phase-locked signal.
# No raw degree values anywhere.

import math as _math

PHASE_CLUSTER_RADIUS_BERU = 1.5 * TIER_A_MAX                 # 1.5 × 0.010 = 0.015 beru
PHASE_CLUSTER_RADIUS      = PHASE_CLUSTER_RADIUS_BERU * BERU  # → degrees (for distance calc)

# Recompute cluster sizes with the 0.015-beru radius for the phase-split test
lons_arr = np.array([s["lon"] for s in sites])
for i, s in enumerate(sites):
    dlon = np.abs(lons_arr - s["lon"])
    dlon = np.minimum(dlon, 360.0 - dlon)
    s["cluster_size_ps"] = int(np.sum(dlon <= PHASE_CLUSTER_RADIUS)) - 1

def anti_node_dev(lon: float, anchor: float = GERIZIM, beru_deg: float = BERU,
                  step: float = 0.1) -> float:
    """Distance in beru from the nearest half-beru midpoint."""
    arc = abs(lon - anchor)
    arc = min(arc, 360.0 - arc)
    bv = arc / beru_deg
    half = step / 2.0
    nearest = round((bv - half) / step) * step + half
    return abs(bv - nearest)

for s in sites:
    s["anti_dev"] = anti_node_dev(s["lon"])
    arc = abs(s["lon"] - GERIZIM)
    arc = min(arc, 360.0 - arc)
    beru_val = arc / BERU
    cell = _math.floor(beru_val / 0.1)
    frac = beru_val - cell * 0.1        # position within 0.1-beru cell [0, 0.1)
    s["coarse_node"] = frac < 0.05      # True = node-side, False = anti-side
    s["anti_aplus"]  = s["anti_dev"] <= TIER_APLUS

node_ap_ps = [s for s in sites if s["coarse_node"] and is_aplus(s["tier"])]
anti_ap_ps = [s for s in sites if not s["coarse_node"] and s["anti_aplus"]]
node_nonap_ps = [s for s in sites if s["coarse_node"] and not is_aplus(s["tier"])]
anti_nonap_ps = [s for s in sites if not s["coarse_node"] and not s["anti_aplus"]]
all_node_ps   = [s for s in sites if s["coarse_node"]]
all_anti_ps   = [s for s in sites if not s["coarse_node"]]

n_node_ap_ps = len(node_ap_ps)
n_anti_ap_ps = len(anti_ap_ps)

c_nap  = [s["cluster_size_ps"] for s in node_ap_ps]
c_aap  = [s["cluster_size_ps"] for s in anti_ap_ps]
c_nna  = [s["cluster_size_ps"] for s in node_nonap_ps]
c_ana  = [s["cluster_size_ps"] for s in anti_nonap_ps]
c_nall = [s["cluster_size_ps"] for s in all_node_ps]
c_aall = [s["cluster_size_ps"] for s in all_anti_ps]

mean_nap  = np.mean(c_nap)  if c_nap  else 0.0
mean_aap  = np.mean(c_aap)  if c_aap  else 0.0
mean_nna  = np.mean(c_nna)  if c_nna  else 0.0
mean_ana  = np.mean(c_ana)  if c_ana  else 0.0
mean_nall = np.mean(c_nall) if c_nall else 0.0
mean_aall = np.mean(c_aall) if c_aall else 0.0

_, p_nap_vs_aap  = mannwhitneyu(c_nap, c_aap, alternative="greater") \
    if (c_nap and c_aap) else (0, 1.0)
_, p_nall_vs_aall = mannwhitneyu(c_nall, c_aall, alternative="greater") \
    if (c_nall and c_aall) else (0, 1.0)

ratio_nap_aap = mean_nap / max(mean_aap, 0.01)

print("=" * 95)
print("  PHASE-SPLIT CLUSTER ANALYSIS (node-side vs anti-node-side)")
print("  Coarse half-cell split: node-side = beru_val mod 0.1 < 0.05,")
print("  anti-side = beru_val mod 0.1 >= 0.05.")
print(f"  Cluster radius: {PHASE_CLUSTER_RADIUS_BERU} beru = 1.5× TIER_A_MAX ({PHASE_CLUSTER_RADIUS:.3f}°).")
print("  Node A+ = node-side  with dev  <= TIER_APLUS;")
print("  Anti A+ = anti-side  with anti_dev <= TIER_APLUS.")
print("=" * 95)
print(f"""
  Node-side A+ sites:      {n_node_ap_ps}
  Anti-side A+ sites:      {n_anti_ap_ps}

  Mean cluster size (±{PHASE_CLUSTER_RADIUS_BERU} beru radius):
    Node A+          :  {mean_nap:.2f}
    Anti A+          :  {mean_aap:.2f}
    Node non-A+      :  {mean_nna:.2f}
    Anti non-A+      :  {mean_ana:.2f}
    All node-side    :  {mean_nall:.2f}
    All anti-side    :  {mean_aall:.2f}

  Node A+ / anti A+ ratio  :  {ratio_nap_aap:.2f}×

  Mann-Whitney U (H₁: node A+ > anti A+):
    p = {p_nap_vs_aap:.4f}  {sig(p_nap_vs_aap)}

  Mann-Whitney U (H₁: all node-side > all anti-side):
    p = {p_nall_vs_aall:.4f}  {sig(p_nall_vs_aall)}

  MANUSCRIPT MACROS (cluster asymmetry / phase-split block):
  \\NclusterNodeAp             {{{n_node_ap_ps}}}
  \\NclusterAntiAp             {{{n_anti_ap_ps}}}
  \\clusterNodeApMean          {{{mean_nap:.2f}}}
  \\clusterAntiApMean          {{{mean_aap:.2f}}}
  \\clusterNodeNonApMean       {{{mean_nna:.2f}}}
  \\clusterAntiNonApMean       {{{mean_ana:.2f}}}
  \\clusterNodeVsAntiApRatio   {{{ratio_nap_aap:.2f}}}
  \\clusterNodeVsAntiApMWp     {{{p_nap_vs_aap:.4f}}}
  \\clusterNodeAllMean         {{{mean_nall:.2f}}}
  \\clusterAntiAllMean         {{{mean_aall:.2f}}}
  \\clusterNodeVsAntiAllMWp    {{{p_nall_vs_aall:.4f}}}
""")


# ── Cluster analysis by harmonic node ────────────────────────────────────────
print("=" * 95)
print("  HARMONIC NODE OCCUPANCY")
print("  How many UNESCO sites fall at each 0.1-beru harmonic node (within ±0.002 beru)?")
print("=" * 95)
print()

# Group A+ sites by their nearest harmonic
harmonic_groups = defaultdict(list)
for s in ap_sites:
    harmonic_groups[s["nearest"]].append(s)

# Also count how many sites are at each harmonic (within A tolerance)
all_at_harmonic = defaultdict(list)
for s in sites:
    if s["dev"] <= TIER_B_MAX:
        all_at_harmonic[s["nearest"]].append(s)

print(f"  {'Harmonic':>8}  {'A+ sites':>8}  {'Total ≤B':>10}  {'A+ names'}")
print(f"  {'-'*8:>8}  {'-'*8:>8}  {'-'*10:>10}  {'-'*50}")
for h in sorted(harmonic_groups):
    ap_names = [s["name"][:40] for s in harmonic_groups[h]]
    total_at = len(all_at_harmonic.get(h, []))
    print(f"  {h:>8.1f}  {len(harmonic_groups[h]):>8d}  {total_at:>10d}  "
          f"{'; '.join(ap_names)}")

# ── Harmonic-level cluster test ──────────────────────────────────────────────
# For each harmonic that has ≥1 A+ site, what is the total number of
# UNESCO sites at that harmonic (within ±B tolerance)?
# Compare to harmonics that have 0 A+ sites.
print()
print("=" * 95)
print("  HARMONIC DENSITY: do A+ harmonics host more total sites?")
print("=" * 95)

ap_harmonic_set = set(s["nearest"] for s in ap_sites)
all_harmonics = set(s["nearest"] for s in sites if s["tier"] != "C")

density_at_ap_harmonics = []
density_at_nonap_harmonics = []

for h in all_harmonics:
    n_at_h = len(all_at_harmonic[h])
    if h in ap_harmonic_set:
        density_at_ap_harmonics.append(n_at_h)
    else:
        density_at_nonap_harmonics.append(n_at_h)

if density_at_ap_harmonics and density_at_nonap_harmonics:
    mean_ap_h = np.mean(density_at_ap_harmonics)
    mean_nonap_h = np.mean(density_at_nonap_harmonics)
    U2, mw2_p = mannwhitneyu(density_at_ap_harmonics, density_at_nonap_harmonics,
                             alternative="greater")
    print(f"""
  Harmonics with ≥1 A+ site:  {len(density_at_ap_harmonics)} harmonics
    Mean total UNESCO sites per harmonic:  {mean_ap_h:.1f}
    Median:  {np.median(density_at_ap_harmonics):.0f}

  Harmonics with 0 A+ sites:  {len(density_at_nonap_harmonics)} harmonics
    Mean total UNESCO sites per harmonic:  {mean_nonap_h:.1f}
    Median:  {np.median(density_at_nonap_harmonics):.0f}

  Mann-Whitney U (H₁: A+ harmonics denser than non-A+ harmonics):
    U = {U2:.0f},  p = {mw2_p:.4f}  {sig(mw2_p)}
    Ratio: {mean_ap_h / max(mean_nonap_h, 0.01):.2f}×
""")

# ── MEGA-CLUSTER identification ─────────────────────────────────────────────
# Find "mega-clusters": locations where ≥3 UNESCO sites share the same
# ~0.5° longitude AND at least one is A+.
print("=" * 95)
print("  MEGA-CLUSTERS AT A+ HARMONICS")
print("  Longitude corridors where ≥3 UNESCO sites share ±0.5° lon")
print("  AND at least one site is Tier-A+")
print("=" * 95)
print()

# Group sites by harmonic
for h in sorted(ap_harmonic_set):
    h_sites = all_at_harmonic[h]
    if len(h_sites) < 3:
        continue
    # Sub-cluster within ±0.5° lon
    ap_at_h = [s for s in h_sites if is_aplus(s["tier"])]
    if not ap_at_h:
        continue
    # Use the first A+ site's longitude as the cluster center
    center_lon = ap_at_h[0]["lon"]
    cluster = [s for s in h_sites
               if min(abs(s["lon"] - center_lon), 360 - abs(s["lon"] - center_lon)) <= CLUSTER_RADIUS_DEG]
    if len(cluster) < 3:
        continue

    print(f"  ── Harmonic {h:.1f} beru │ center lon ≈ {center_lon:.2f}° │ {len(cluster)} sites ──")
    for s in sorted(cluster, key=lambda x: x["dev"]):
        flag = " ◀◀ A+" if is_aplus(s["tier"]) else (" ◀ A" if s["tier"] == "A" else "")
        print(f"    {s['name'][:55]:<55s}  lon={s['lon']:>8.4f}  dev={s['dev']:.5f}  "
              f"{s['tier']}{flag}")
    print()

# ── Summary ──────────────────────────────────────────────────────────────────
print("=" * 95)
print("  SUMMARY")
print("=" * 95)
print(f"""
  QUESTION: Do UNESCO sites at beru-harmonic longitudes cluster more
  densely than sites at non-harmonic longitudes?

  FULL UNESCO CORPUS: N = {N} Cultural/Mixed sites with coordinates
  Tier-A+ sites: {n_Ap}  ({100*n_Ap/N:.1f}%)

  SITE-LEVEL CLUSTER SIZE (other sites within ±{CLUSTER_RADIUS_DEG:.3f}° lon = ±{CLUSTER_RADIUS_DEG/BERU:.3f} beru):
    A+ mean:     {mean_ap:.2f}
    non-A+ mean: {mean_nonap:.2f}
    Ratio:       {mean_ap / max(mean_nonap, 0.01):.2f}×
    Mann-Whitney p = {mw_p:.4f}  {sig(mw_p)}
    Permutation  p = {perm_p:.4f}  {sig(perm_p)}

  MANUSCRIPT MACROS (full-corpus cluster block):
  \\clusterApMean            {{{mean_ap:.2f}}}
  \\clusterNonApMean         {{{mean_nonap:.2f}}}
  \\clusterRatio             {{{mean_ap / max(mean_nonap, 0.01):.2f}}}

  INTERPRETATION:
  If A+ sites cluster significantly more than non-A+ sites, it means
  beru-harmonic longitudes are not randomly distributed across the globe's
  monument density — they preferentially mark locations where humans
  built MANY monuments.  This is the predicted signature of a geodetic
  system: a founding monument is placed at a harmonic longitude, and
  subsequent monuments cluster around it by cultural gravity.

  Under the null hypothesis (longitudes uniformly distributed w.r.t.
  harmonics), there is no reason for A+ longitudes to be denser than
  non-A+ longitudes.  Any asymmetry is evidence that the grid harmonics
  mark real "places" — the founding nodes of monument traditions.

  This test does NOT require one-per-complex deduplication.
  Instead of discarding the clusters, it USES them as the signal.
""")

# ── LaTeX macros (GROUP 2) ────────────────────────────────────────────────────
perm_z = (mean_ap - np.mean(perm_means)) / max(np.std(perm_means), 0.001)
n_ap_harmonics = len(density_at_ap_harmonics)
n_nonap_harmonics = len(density_at_nonap_harmonics)
mean_ap_h = np.mean(density_at_ap_harmonics) if density_at_ap_harmonics else 0
mean_nonap_h = np.mean(density_at_nonap_harmonics) if density_at_nonap_harmonics else 0

print("  % LaTeX macros (GROUP 2):")
print(f"  \\newcommand{{\\NclusterTotal}}{{{N}}}          % full Cultural/Mixed corpus")
print(f"  \\newcommand{{\\NclusterAp}}{{{n_Ap}}}              % Tier-A+ sites")
print(f"  \\newcommand{{\\clusterApRate}}{{{100*n_Ap/N:.1f}}}           % A+ rate (%)")
print(f"  \\newcommand{{\\clusterApBinom}}{{{bt_Ap.pvalue:.4f}}}        % binomial p, A+")
print(f"  \\newcommand{{\\clusterApMean}}{{{mean_ap:.2f}}}          % mean cluster size, A+ sites")
print(f"  \\newcommand{{\\clusterNonApMean}}{{{mean_nonap:.2f}}}          % mean cluster size, non-A+ sites")
print(f"  \\newcommand{{\\clusterRatio}}{{{mean_ap / max(mean_nonap, 0.01):.2f}}}           % A+/non-A+ cluster size ratio")
print(f"  \\newcommand{{\\clusterMWp}}{{{mw_p:.4f}}}        % Mann-Whitney p (site-level)")
print(f"  \\newcommand{{\\clusterPermP}}{{{perm_p:.4f}}}        % permutation p (site-level)")
print(f"  \\newcommand{{\\clusterPermZ}}{{{perm_z:.2f}}}          % permutation Z-score")
print(f"  \\newcommand{{\\NclusterApHarmonics}}{{{n_ap_harmonics}}}              % harmonics with ≥1 A+ site")
print(f"  \\newcommand{{\\NclusterNonApHarmonics}}{{{n_nonap_harmonics}}}              % harmonics with 0 A+ sites")
print(f"  \\newcommand{{\\clusterHarmonicApMean}}{{{mean_ap_h:.1f}}}          % mean UNESCO sites/harmonic (A+ harmonics)")
print(f"  \\newcommand{{\\clusterHarmonicNonApMean}}{{{mean_nonap_h:.1f}}}           % mean UNESCO sites/harmonic (non-A+ harmonics)")
print(f"  \\newcommand{{\\clusterHarmonicRatio}}{{{mean_ap_h / max(mean_nonap_h, 0.01):.2f}}}           % harmonic-level density ratio")
print(f"  \\newcommand{{\\clusterHarmonicMWp}}{{{mw2_p:.4f}}}        % Mann-Whitney p (harmonic-level)")
print(f"  \\newcommand{{\\NclusterNodeAp}}{{{n_node_ap_ps}}}              % node-side A+ sites (phase-split)")
print(f"  \\newcommand{{\\NclusterAntiAp}}{{{n_anti_ap_ps}}}              % anti-side A+ sites (phase-split)")
print(f"  \\newcommand{{\\clusterNodeApMean}}{{{mean_nap:.2f}}}          % mean cluster size, node-side A+")
print(f"  \\newcommand{{\\clusterAntiApMean}}{{{mean_aap:.2f}}}          % mean cluster size, anti-side A+")
print(f"  \\newcommand{{\\clusterNodeNonApMean}}{{{mean_nna:.2f}}}          % mean cluster size, node-side non-A+")
print(f"  \\newcommand{{\\clusterAntiNonApMean}}{{{mean_ana:.2f}}}          % mean cluster size, anti-side non-A+")
print(f"  \\newcommand{{\\clusterNodeVsAntiApRatio}}{{{ratio_nap_aap:.2f}}}           % node A+ / anti A+ cluster ratio")
print(f"  \\newcommand{{\\clusterNodeVsAntiApMWp}}{{{p_nap_vs_aap:.4f}}}        % Mann-Whitney p, node A+ > anti A+")
print(f"  \\newcommand{{\\clusterNodeAllMean}}{{{mean_nall:.2f}}}          % mean cluster size, all node-side")
print(f"  \\newcommand{{\\clusterAntiAllMean}}{{{mean_aall:.2f}}}          % mean cluster size, all anti-side")
print(f"  \\newcommand{{\\clusterNodeVsAntiAllMWp}}{{{p_nall_vs_aall:.4f}}}        % Mann-Whitney p, all node > all anti")

# ── Write to results store ────────────────────────────────────────────────────
ResultsStore().write_many({
    "clusterApBinom":     bt_Ap.pvalue,      # binomial p, full corpus A+ — Test 1
    "clusterPermP":       perm_p,            # permutation p, cluster asymmetry — Test 3
    "clusterMWp":         mw_p,              # Mann-Whitney p (site-level)
    "clusterHarmonicMWp": mw2_p,             # Mann-Whitney p (harmonic-level)
    "clusterPermZ":       round(perm_z, 4),  # permutation Z-score
})


# ══════════════════════════════════════════════════════════════════════════════
# SUBTEST: Conditional sensitivity-peak test
# ══════════════════════════════════════════════════════════════════════════════
#
# The unit-sensitivity peak (±0.3% collapse at 0.1000 beru) is not unusual
# under a simple longitude-shuffle null (permSlopePsharp, ns; Test I).
# The deeper question is: does clustering *itself* produce the peak?
#
# We test this by conditioning: among permutations that reproduce clustering
# at least as extreme as observed, what fraction also show the canonical-unit
# specificity (0.1000 beru is the best joint-significance spacing)?
#
# Clustering null: exactly the same as the permutation block above — shuffle
# which n_Ap sites are called A+ (fixed real cluster sizes, random assignment).
# Fine-sweep null: independently shuffle longitudes (same as sensitivity_slope_*.py).
# Both shuffles run in the same loop; they use the same RNG but are independent
# draws (cluster draw uses choice(), fine-sweep uses permutation()).
#
# Interpretation:
#   conditional fraction >> unconditional → clustering drives the peak
#   conditional fraction ≈ unconditional  → peak is independent of clustering
# ══════════════════════════════════════════════════════════════════════════════

from scipy.stats import binomtest as _binomtest
from lib.dome_filter import is_dome_site

# Fine sweep configuration (must match sensitivity_slope_specificity_test.py)
FINE_SPACINGS_SS = [0.0990, 0.0993, 0.0995, 0.1000, 0.1001, 0.1003, 0.1010]
CANONICAL_SS     = 0.1000
THRESH_SS        = 0.05

# Load site longitudes for the fine-sweep half (reuse cultural sites already loaded)
_cultural_ss = cultural_sites_with_coords(load_corpus())
dome_mask_ss = np.array([is_dome_site(s) for s in _cultural_ss])
all_lons_ss  = np.array([s.longitude     for s in _cultural_ss])

def _null_rate(spacing):
    return 2 * TIER_APLUS / spacing

def _count_hits_vec(lons_arr: np.ndarray, spacing: float) -> int:
    """Vectorized version of _count_hits — no Python loop over sites."""
    arc = np.abs(lons_arr - GERIZIM)
    beru_val = arc / BERU
    nearest = np.round(beru_val / spacing) * spacing
    dev = np.abs(beru_val - nearest)
    return int(np.sum(dev <= TIER_APLUS))

def _canon_is_best_joint(all_lons: np.ndarray, dome_lons: np.ndarray) -> bool:
    """Return True iff 0.1000 beru is the unique best joint-significance spacing.
    Uses vectorized hit-counting for speed."""
    best_comb = -1.0
    best_sp   = None
    n_all  = len(all_lons)
    n_dome = len(dome_lons)
    for sp in FINE_SPACINGS_SS:
        nr   = _null_rate(sp)
        d_p  = _binomtest(_count_hits_vec(dome_lons, sp), n_dome, nr, alternative='greater').pvalue
        f_p  = _binomtest(_count_hits_vec(all_lons,  sp), n_all,  nr, alternative='greater').pvalue
        comb = -np.log10(max(d_p, 1e-15)) + -np.log10(max(f_p, 1e-15))
        if (d_p < THRESH_SS and f_p < THRESH_SS) and comb > best_comb:
            best_comb = comb
            best_sp   = sp
    return best_sp == CANONICAL_SS

# The clustering null is IDENTICAL to the permutation block above:
# draw n_Ap random indices from all_cluster_sizes (fixed real cluster sizes),
# check whether the drawn mean >= mean_ap (observed A+ mean cluster size).
# The fine-sweep null independently shuffles all_lons_ss with the same RNG.
all_cluster_sizes_ss = all_cluster_sizes   # numpy array from cluster block above (line ~228)

N_COND_PERMS = 50_000
SEED_COND    = 99
rng_cond     = np.random.default_rng(SEED_COND)

n_cluster_ok       = 0   # perms where cluster mean >= mean_ap
n_cluster_and_peak = 0   # perms where cluster ok AND canonical is best
n_peak_total       = 0   # unconditional canonical-best count

print()
print("=" * 95)
print("  SUBTEST: Conditional sensitivity-peak test")
print(f"  Q: When permutations reproduce the observed A+ clustering")
print(f"     (drawn mean cluster size >= {mean_ap:.2f}), do they also produce")
print(f"     the unit-sensitivity peak at 0.1000 beru more often?")
print(f"  Clustering null : shuffle tier labels, fixed cluster sizes (same as above).")
print(f"  Fine-sweep null : independently shuffle longitudes.")
print(f"  N permutations  : {N_COND_PERMS:,}  (seed={SEED_COND})")
print("=" * 95)

for i in range(N_COND_PERMS):
    # Clustering: draw n_Ap random indices from fixed cluster sizes
    drawn_idx   = rng_cond.choice(len(all_cluster_sizes_ss), size=n_Ap, replace=False)
    perm_c_mean = float(all_cluster_sizes_ss[drawn_idx].mean())
    cluster_ok  = (perm_c_mean >= mean_ap)

    # Fine sweep: shuffle longitudes independently
    perm_lons      = rng_cond.permutation(all_lons_ss)
    perm_dome_lons = perm_lons[dome_mask_ss]
    peak_ok = _canon_is_best_joint(perm_lons, perm_dome_lons)

    if cluster_ok:
        n_cluster_ok += 1
        if peak_ok:
            n_cluster_and_peak += 1
    if peak_ok:
        n_peak_total += 1

    if (i + 1) % 5000 == 0:
        cond_frac = n_cluster_and_peak / n_cluster_ok if n_cluster_ok else float('nan')
        print(f"  {i+1:>6}/{N_COND_PERMS}  "
              f"cluster_ok={n_cluster_ok}  "
              f"cluster_and_peak={n_cluster_and_peak}  "
              f"cond_frac={cond_frac:.3f}")

p_cluster_alone = n_cluster_ok   / N_COND_PERMS
p_peak_uncond   = n_peak_total   / N_COND_PERMS
p_peak_cond     = (n_cluster_and_peak / n_cluster_ok) if n_cluster_ok > 0 else float('nan')

def _sig(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"

print(f"""
  ── Results ──
  Total permutations:                        {N_COND_PERMS:,}
  Permutations with clustering >= observed:  {n_cluster_ok}  (p = {p_cluster_alone:.4f})
  Of those, canonical peak also present:     {n_cluster_and_peak}
    Conditional fraction (peak | cluster):   {p_peak_cond:.4f}  {_sig(p_peak_cond)}
  Unconditional peak fraction (reference):   {n_peak_total}/{N_COND_PERMS}  = {p_peak_uncond:.4f}

  INTERPRETATION:
  If conditional fraction >> unconditional fraction, clustering drives the peak.
  If roughly equal, clustering adds no extra specificity to the unit result.
""")

# Significance label strings for macros
sig_cond  = _sig(p_peak_cond)
sig_uncond = _sig(p_peak_uncond)

print("  % LaTeX macros (cluster conditional sensitivity subtest):")
print(f"  \\newcommand{{\\clusterCondN}}{{{N_COND_PERMS}}}          % perms in conditional subtest")
print(f"  \\newcommand{{\\clusterCondClusterOk}}{{{n_cluster_ok}}}    % perms reproducing observed clustering")
print(f"  \\newcommand{{\\clusterCondPeakAndCluster}}{{{n_cluster_and_peak}}}    % perms with clustering AND unit peak")
print(f"  \\newcommand{{\\clusterCondPeakFrac}}{{{p_peak_cond:.4f}}}    % conditional fraction (peak | cluster)")
print(f"  \\newcommand{{\\clusterCondPeakFracUncond}}{{{p_peak_uncond:.4f}}}    % unconditional peak fraction")
print(f"  \\newcommand{{\\clusterCondPeakSig}}{{{sig_cond}}}    % significance label, conditional fraction")

ResultsStore().write_many({
    "clusterCondN":              N_COND_PERMS,
    "clusterCondClusterOk":      n_cluster_ok,
    "clusterCondPeakAndCluster": n_cluster_and_peak,
    "clusterCondPeakFrac":       round(p_peak_cond,   4),
    "clusterCondPeakFracUncond": round(p_peak_uncond, 4),
    "clusterCondPeakSig":        sig_cond,
})
