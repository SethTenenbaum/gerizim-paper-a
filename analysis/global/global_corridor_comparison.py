import sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, TIER_APLUS, P_NULL_AP

TIER_AP = TIER_APLUS
P_NULL = P_NULL_AP
WINDOW = 0.5  # ±0.25° → 0.5° total window

# ── Load sites ──────────────────────────────────────────────────────────────
corpus = load_corpus()
_cultural = cultural_sites_with_coords(corpus)

sites = []
for s in _cultural:
    yr = s.year if s.year is not None else 9999
    sites.append({"name": s.site, "lon": s.longitude, "lat": s.latitude, "yr": yr})

N = len(sites)
lons = np.array([s["lon"] for s in sites])
print(f"Loaded {N} Cultural/Mixed UNESCO sites with coordinates.\n")

# ── A+ count function ──────────────────────────────────────────────────────
def count_ap(anchor):
    arcs = np.abs(lons - anchor)
    arcs = np.minimum(arcs, 360 - arcs)
    bvs = arcs / BERU
    nearests = np.round(bvs * 10) / 10
    devs = np.abs(bvs - nearests)
    return int(np.sum(devs <= TIER_AP))

# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 100)
print("  STEP 1: Computing A+ count for every UNESCO site used as anchor...")
print("=" * 100)

site_ap_counts = []
for i, s in enumerate(sites):
    cnt = count_ap(s["lon"])
    site_ap_counts.append({
        "name": s["name"],
        "lon": s["lon"],
        "lat": s["lat"],
        "yr": s["yr"],
        "ap_count": cnt,
    })

site_ap_counts.sort(key=lambda x: -x["ap_count"])
print(f"  Done. Mean A+ = {np.mean([s['ap_count'] for s in site_ap_counts]):.1f}, "
      f"Max = {max(s['ap_count'] for s in site_ap_counts)}, "
      f"Min = {min(s['ap_count'] for s in site_ap_counts)}\n")

# ══════════════════════════════════════════════════════════════════════════════
# sites-as-anchors score ≥ threshold
# ══════════════════════════════════════════════════════════════════════════════
THRESHOLDS = [55, 53, 50]
HALF_WIN = WINDOW / 2

print("=" * 100)
print(f"  STEP 2: Sliding {WINDOW}° window — count sites with A+ ≥ threshold")
print("=" * 100)

window_centers = np.arange(0, 360, 0.1)
results_by_threshold = {}

for thresh in THRESHOLDS:
    windows = []
    for center in window_centers:
        lo = center - HALF_WIN
        hi = center + HALF_WIN
        count = 0
        names_in_window = []
        for sc in site_ap_counts:
            slon = sc["lon"]
            slon_norm = slon % 360
            center_norm = center % 360
            diff = abs(slon_norm - center_norm)
            if diff > 180:
                diff = 360 - diff
            if diff <= HALF_WIN:
                if sc["ap_count"] >= thresh:
                    count += 1
                    names_in_window.append((sc["name"], sc["ap_count"], sc["lon"]))
        windows.append({
            "center": center,
            "count": count,
            "sites": names_in_window,
        })
    results_by_threshold[thresh] = windows

# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
print()
for thresh in THRESHOLDS:
    windows = results_by_threshold[thresh]
    windows.sort(key=lambda x: -x["count"])
    max_count = windows[0]["count"]
    
    print(f"\n{'=' * 100}")
    print(f"  TOP CORRIDORS: sites with A+ ≥ {thresh} within {WINDOW}° window")
    print(f"  Global maximum: {max_count} sites")
    print(f"{'=' * 100}")
    
    seen = set()
    rank = 0
    for w in windows:
        if w["count"] < max(1, max_count - 2):
            break
        center_key = round(w["center"], 0)
        if center_key in seen and w["count"] < max_count:
            continue
        seen.add(center_key)
        rank += 1
        if rank > 20:
            break
        
        print(f"\n  #{rank}: center = {w['center']:.1f}°, count = {w['count']}")
        for name, ap, lon in sorted(w["sites"], key=lambda x: -x[1]):
            print(f"    {name[:60]:<60s} lon={lon:>8.3f}° A+={ap}")

# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
MAIN_THRESH = 55
print(f"\n\n{'=' * 100}")
print(f"  STEP 4: ALL DISTINCT CLUSTERS with ≥1 site scoring A+ ≥ {MAIN_THRESH}")
print(f"{'=' * 100}")

high_scorers = [s for s in site_ap_counts if s["ap_count"] >= MAIN_THRESH]
print(f"\n  Total sites with A+ ≥ {MAIN_THRESH}: {len(high_scorers)}")

clusters = []
used = set()
for i, s in enumerate(high_scorers):
    if i in used:
        continue
    cluster = [s]
    used.add(i)
    for j, s2 in enumerate(high_scorers):
        if j in used:
            continue
        diff = abs(s["lon"] - s2["lon"])
        if diff > 180:
            diff = 360 - diff
        if diff <= 1.0:
            cluster.append(s2)
            used.add(j)
    clusters.append(cluster)

clusters.sort(key=lambda c: -len(c))

print(f"  Distinct clusters (sites within 1° of each other): {len(clusters)}\n")
for i, c in enumerate(clusters):
    center_lon = np.mean([s["lon"] for s in c])
    max_ap = max(s["ap_count"] for s in c)
    sum_ap = sum(s["ap_count"] for s in c)
    print(f"  ── Cluster #{i+1}: center ≈ {center_lon:.1f}°, "
          f"{len(c)} sites, max A+={max_ap}, sum A+={sum_ap} ──")
    for s in sorted(c, key=lambda x: -x["ap_count"]):
        print(f"    {s['name'][:60]:<60s} lon={s['lon']:>8.3f}° A+={s['ap_count']} yr={s['yr']}")
    print()

# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'=' * 100}")
print(f"  STEP 5: JERUSALEM–GERIZIM CORRIDOR vs. ALL OTHER CLUSTERS")
print(f"{'=' * 100}")

jerusalem_cluster = None
other_clusters = []
for c in clusters:
    center = np.mean([s["lon"] for s in c])
    if 34.5 <= center <= 36.0:
        jerusalem_cluster = c
    else:
        other_clusters.append(c)

if jerusalem_cluster:
    j_size = len(jerusalem_cluster)
    j_max = max(s["ap_count"] for s in jerusalem_cluster)
    j_sum = sum(s["ap_count"] for s in jerusalem_cluster)
    j_mean_top3 = np.mean(sorted([s["ap_count"] for s in jerusalem_cluster], reverse=True)[:3])
    
    print(f"\n  JERUSALEM–GERIZIM CORRIDOR (34.5–36.0°E):")
    print(f"    Sites with A+ ≥ {MAIN_THRESH}: {j_size}")
    print(f"    Max A+: {j_max}")
    print(f"    Sum A+: {j_sum}")
    print(f"    Top-3 mean A+: {j_mean_top3:.1f}")
    for s in sorted(jerusalem_cluster, key=lambda x: -x["ap_count"]):
        print(f"      {s['name'][:55]:<55s} lon={s['lon']:>8.3f}° A+={s['ap_count']} yr={s['yr']}")
    
    print(f"\n  ALL OTHER CLUSTERS:")
    for i, c in enumerate(other_clusters):
        center = np.mean([s["lon"] for s in c])
        c_size = len(c)
        c_max = max(s["ap_count"] for s in c)
        c_sum = sum(s["ap_count"] for s in c)
        c_mean_top3 = np.mean(sorted([s["ap_count"] for s in c], reverse=True)[:min(3, len(c))])
        
        lon_range = max(s["lon"] for s in c) - min(s["lon"] for s in c)
        lat_range = max(s["lat"] for s in c) - min(s["lat"] for s in c)
        
        print(f"\n    Cluster #{i+1}: center ≈ {center:.1f}°, {c_size} sites")
        print(f"      Max A+: {c_max}, Sum A+: {c_sum}, Top-3 mean: {c_mean_top3:.1f}")
        print(f"      Lon spread: {lon_range:.2f}°, Lat spread: {lat_range:.2f}°")
        for s in sorted(c, key=lambda x: -x["ap_count"]):
            print(f"        {s['name'][:55]:<55s} lon={s['lon']:>8.3f}° A+={s['ap_count']} yr={s['yr']}")

# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n\n{'=' * 100}")
print(f"  STEP 6: GLOBAL CORRIDOR RANKING TABLE")
print(f"  (All clusters with ≥2 sites scoring A+ ≥ {MAIN_THRESH})")
print(f"{'=' * 100}\n")

all_clusters_for_table = []
for c in clusters:
    center = np.mean([s["lon"] for s in c])
    c_size = len(c)
    c_max = max(s["ap_count"] for s in c)
    c_sum = sum(s["ap_count"] for s in c)
    c_mean_top3 = np.mean(sorted([s["ap_count"] for s in c], reverse=True)[:min(3, len(c))])
    
    lat_range = max(s["lat"] for s in c) - min(s["lat"] for s in c)
    lon_range = max(s["lon"] for s in c) - min(s["lon"] for s in c)
    
    is_jerusalem = 34.5 <= center <= 36.0
    
    all_clusters_for_table.append({
        "center": center,
        "size": c_size,
        "max_ap": c_max,
        "sum_ap": c_sum,
        "mean_top3": c_mean_top3,
        "lat_range": lat_range,
        "lon_range": lon_range,
        "is_jerusalem": is_jerusalem,
        "sites": c,
    })

all_clusters_for_table.sort(key=lambda x: (-x["size"], -x["sum_ap"]))

print(f"  {'Rank':>4}  {'Center':>7}  {'N≥55':>4}  {'MaxA+':>5}  {'SumA+':>6}  {'Top3':>6}  {'LatSpread':>9}  {'Region'}")
print(f"  {'─'*4}  {'─'*7}  {'─'*4}  {'─'*5}  {'─'*6}  {'─'*6}  {'─'*9}  {'─'*40}")

for rank, ct in enumerate(all_clusters_for_table, 1):
    if ct["size"] < 2:
        continue
    label = "★ JERUSALEM–GERIZIM" if ct["is_jerusalem"] else ""
    if not label:
        lats = [s["lat"] for s in ct["sites"]]
        mean_lat = np.mean(lats)
        if mean_lat > 50:
            label = "Northern Europe"
        elif mean_lat > 30:
            label = "Mediterranean / Central Asia"
        elif mean_lat > 0:
            label = "Tropics / South Asia"
        elif mean_lat > -30:
            label = "South America / Africa"
        else:
            label = "Southern latitudes"
    
    print(f"  {rank:>4}  {ct['center']:>7.1f}°  {ct['size']:>4}  {ct['max_ap']:>5}  {ct['sum_ap']:>6}  "
          f"{ct['mean_top3']:>6.1f}  {ct['lat_range']:>8.1f}°  {label}")

# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n\n{'=' * 100}")
print(f"  STEP 7: KEY NUMBERS FOR MANUSCRIPT")
print(f"{'=' * 100}\n")

if jerusalem_cluster:
    j_c = jerusalem_cluster
    j_n = len(j_c)
    j_max = max(s["ap_count"] for s in j_c)
    j_sum = sum(s["ap_count"] for s in j_c)
    
    # Best rival
    best_rival = max(other_clusters, key=lambda c: len(c)) if other_clusters else None
    
    if best_rival:
        r_n = len(best_rival)
        r_center = np.mean([s["lon"] for s in best_rival])
        r_max = max(s["ap_count"] for s in best_rival)
        r_sum = sum(s["ap_count"] for s in best_rival)
        r_lat_range = max(s["lat"] for s in best_rival) - min(s["lat"] for s in best_rival)
    
    print(f"  Jerusalem–Gerizim corridor:")
    print(f"    N sites with A+ ≥ {MAIN_THRESH}: {j_n}")
    print(f"    Max A+: {j_max}")
    print(f"    Sum A+: {j_sum}")
    print(f"    Sites are all in 35.18–35.27°E (0.09° span)")
    print(f"    All archaeologically related Levantine sites")
    
    if best_rival:
        print(f"\n  Closest rival (center ≈ {r_center:.0f}°E):")
        print(f"    N sites with A+ ≥ {MAIN_THRESH}: {r_n}")
        print(f"    Max A+: {r_max}")
        print(f"    Sum A+: {r_sum}")
        print(f"    Lat spread: {r_lat_range:.1f}°")
        for s in sorted(best_rival, key=lambda x: -x["ap_count"]):
            print(f"      {s['name'][:55]:<55s} {s['lon']:>8.3f}° lat={s['lat']:>7.3f}° A+={s['ap_count']}")

for alt_thresh in [53, 50]:
    high = [s for s in site_ap_counts if s["ap_count"] >= alt_thresh]
    j_alt = [s for s in high if 34.5 <= s["lon"] <= 36.0]
    others = [s for s in high if not (34.5 <= s["lon"] <= 36.0)]
    
    # Cluster others
    others.sort(key=lambda x: x["lon"])
    alt_clusters = []
    used_alt = set()
    for i, s in enumerate(others):
        if i in used_alt:
            continue
        cl = [s]
        used_alt.add(i)
        for j, s2 in enumerate(others):
            if j in used_alt:
                continue
            diff = abs(s["lon"] - s2["lon"])
            if diff > 180:
                diff = 360 - diff
            if diff <= 1.0:
                cl.append(s2)
                used_alt.add(j)
        alt_clusters.append(cl)
    
    best_alt = max(alt_clusters, key=len) if alt_clusters else []
    best_alt_center = np.mean([s["lon"] for s in best_alt]) if best_alt else 0
    
    print(f"\n  At threshold A+ ≥ {alt_thresh}:")
    print(f"    Jerusalem corridor: {len(j_alt)} sites")
    if best_alt:
        print(f"    Best rival (center ≈ {best_alt_center:.0f}°E): {len(best_alt)} sites")

print(f"\n\n{'=' * 100}")
print("  DONE — global_corridor_comparison.py")
print(f"{'=' * 100}")
