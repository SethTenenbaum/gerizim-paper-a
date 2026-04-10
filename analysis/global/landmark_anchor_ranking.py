
import sys
from pathlib import Path
import numpy as np
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_B_MAX, P_NULL_AP,
    CONFIG, load_notable_anchors,
)

GERIZIM_TEMPLE = GERIZIM  # Tentative list ref. 5706, loaded from config.json
TEL_MEGIDDO    = CONFIG["anchors"]["megiddo"]["longitude"]

# ── Load sites ──────────────────────────────────────────────────────────────
corpus = load_corpus()
_cultural = cultural_sites_with_coords(corpus)
sites = [{"name": s.site, "lon": s.longitude, "lat": s.latitude} for s in _cultural]

N = len(sites)
lons = np.array([s["lon"] for s in sites])

def count_ap_at_anchor(anchor):
    arcs = np.abs(lons - anchor)
    arcs = np.minimum(arcs, 360 - arcs)
    bvs = arcs / BERU
    nearests = np.round(bvs * 10) / 10
    devs = np.abs(bvs - nearests)
    return int(np.sum(devs <= TIER_APLUS))

# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 110)
print("  LANDMARK-AS-ANCHOR RANKING")
print(f"  Use each of {N} UNESCO Cultural/Mixed sites as anchor, count A+ sites")
print(f"  Tier-A+ threshold: ≤ {TIER_APLUS} beru (≤ 6.7 km)")
print("=" * 110)

results = []
for s in sites:
    ap = count_ap_at_anchor(s["lon"])
    phase = s["lon"] % 3.0
    results.append({
        "name": s["name"],
        "lon": s["lon"],
        "ap": ap,
        "phase": phase,
    })

results.sort(key=lambda x: -x["ap"])

ap_vals = np.array([r["ap"] for r in results])
mean_ap = float(np.mean(ap_vals))
std_ap = float(np.std(ap_vals))
max_ap = int(np.max(ap_vals))
min_ap = int(np.min(ap_vals))

# ── Compute specific anchor stats ──────────────────────────────────────────
gerizim_temple_ap = count_ap_at_anchor(GERIZIM_TEMPLE)
gerizim_summit_ap = count_ap_at_anchor(GERIZIM)
megiddo_ap = count_ap_at_anchor(TEL_MEGIDDO)

# Rank: count how many sites have A+ >= this anchor's A+ count
gerizim_rank = int(np.sum(ap_vals >= gerizim_temple_ap))
megiddo_rank = int(np.sum(ap_vals >= megiddo_ap))

# Count sites at x.18° artifact and in Gerizim corridor among top-50
n_artifact = sum(1 for r in results[:50] if abs(r["phase"] - 2.18) < 0.05)
n_gerizim_corridor = sum(1 for r in results[:50] if 35.15 <= r["lon"] <= 35.30)

print(f"""
  DISTRIBUTION SUMMARY
  ──────────────────────────────────────────────────────────────
  Mean A+ count:    {mean_ap:.1f} ± {std_ap:.1f}
  Max A+ count:     {max_ap}
  Min A+ count:     {min_ap}
  Expected (4%):    {0.04 * N:.1f}
  Gerizim temple ({GERIZIM_TEMPLE}°E):
    A+ count:   {gerizim_temple_ap}
    Rank:       #{gerizim_rank} of {N} (top {100*gerizim_rank/N:.1f}%)
    Phase:      {GERIZIM_TEMPLE % 3.0:.3f}° (mod 3°)
    
  Gerizim summit ({GERIZIM}°E):
    A+ count:   {gerizim_summit_ap}
    Phase:      {GERIZIM % 3.0:.3f}° (mod 3°)
    
  Tel Megiddo ({TEL_MEGIDDO}°E):
    A+ count:   {megiddo_ap}
    Rank:       #{megiddo_rank} of {N} (top {100*megiddo_rank/N:.1f}%)
    Phase:      {TEL_MEGIDDO % 3.0:.3f}° (mod 3°)

  Among top 50 sites by A+ count:
    Sites at x.18° artifact (phase within 0.05° of 2.18):  {n_artifact}
    Sites in Gerizim corridor (35.15–35.30°E):              {n_gerizim_corridor}
    Other:                                                  {50 - n_artifact - n_gerizim_corridor}
    
  Note: Many high-A+ sites are at the x.18° phase artifact which spans
  continents — they are mathematical, not archaeological.
""")

# ── Top-50 table ──────────────────────────────────────────────────────────
print(f"\n  {'Rank':>4}  {'A+':>4}  {'Phase':>6}  {'Lon':>8}  {'Site':<60}")
print("  " + "─" * 90)
for i, r in enumerate(results[:50]):
    marker = " ★" if abs(r["lon"] - GERIZIM) < 0.1 else ""
    art = " [x.18°]" if abs(r["phase"] - 2.18) < 0.05 else ""
    print(f"  {i+1:>4}  {r['ap']:>4}  {r['phase']:>6.3f}  {r['lon']:>8.3f}  {r['name'][:58]:<60}{marker}{art}")

# ── Gerizim rank among non-artifact sites ─────────────────────────────────
non_artifact = [r for r in results if abs(r["phase"] - 2.18) >= 0.05]
non_art_ap = np.array([r["ap"] for r in non_artifact])
gerizim_rank_noart = int(np.sum(non_art_ap >= gerizim_temple_ap))
print(f"\n  Gerizim rank excluding x.18° artifact sites: #{gerizim_rank_noart} of {len(non_artifact)}")
print(f"  (top {100*gerizim_rank_noart/len(non_artifact):.1f}%)")

# ── LaTeX macros (GROUP 13) ───────────────────────────────────────────────────
print("  % LaTeX macros (GROUP 13):")
print(f"  \\newcommand{{\\GerizimTempleAp}}{{{gerizim_temple_ap}}}          % A+ count at Gerizim Temple anchor")
print(f"  \\newcommand{{\\GerizimTempleRank}}{{{gerizim_rank}}}          % rank of Gerizim Temple among all landmarks")

print(f"\n{'='*110}")
print("  DONE — landmark_anchor_ranking.py")
print(f"{'='*110}")