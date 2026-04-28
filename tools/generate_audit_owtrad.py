"""
generate_audit_owtrad.py
========================
Produce supplementary/audit/owtrad_audit.txt

Full reproducibility audit for the OWTRAD Silk Road harmonic-tier analysis.
Lists every A++ and A+ vertex, all population counts, binomial/Z/permutation
test results, and cluster asymmetry test — enough for independent replication.

Run from repo root:
    python3 tools/generate_audit_owtrad.py
"""

import csv
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from math import sqrt
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from lib.beru import (
    GERIZIM, BERU, HARMONIC_STEP,
    TIER_APP, TIER_APLUS, TIER_A_MAX,
    P_NULL_APP, P_NULL_AP, P_NULL_A,
)

SILK_ROAD_DIR = ROOT / "data" / "store" / "silk_road"
OUT_PATH      = ROOT / "supplementary" / "audit" / "owtrad_audit.txt"
N_PERMS       = 100_000
RNG_SEED      = 42
BONF_K        = 4


# ── helpers ───────────────────────────────────────────────────────────────────

def _sig(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"

def _zscore(obs, n, p_null):
    mu = n * p_null
    sd = sqrt(n * p_null * (1 - p_null))
    return (obs - mu) / sd if sd > 0 else 0.0

def _beru_dev(lon, anchor=GERIZIM):
    bv   = (lon - anchor) / BERU
    near = round(bv / HARMONIC_STEP) * HARMONIC_STEP
    return abs(bv - near)

def _beru_devs_np(lons, anchor=GERIZIM):
    bv   = (lons - anchor) / BERU
    near = np.round(bv / HARMONIC_STEP) * HARMONIC_STEP
    return np.abs(bv - near)

def _tier(dev):
    if dev <= TIER_APP:   return "A++"
    if dev <= TIER_APLUS: return "A+"
    if dev <= TIER_A_MAX: return "A"
    return "-"


# ── load data ─────────────────────────────────────────────────────────────────

def _load_nodes():
    """Return list of dicts with name, lon, lat, country, dataset."""
    nodes = []
    with open(SILK_ROAD_DIR / "owtrad_nodes.csv", newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            try:
                nodes.append({
                    "name":    row.get("name", ""),
                    "lon":     float(row["lon"]),
                    "lat":     float(row["lat"]),
                    "country": row.get("country", ""),
                    "dataset": row.get("dataset", ""),
                })
            except (ValueError, KeyError):
                continue
    return nodes

def _load_edges():
    edges = []
    with open(SILK_ROAD_DIR / "owtrad_routes.csv", newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            try:
                edges.append({
                    "node1": row.get("node1", ""), "lon1": float(row["lon1"]),
                    "node2": row.get("node2", ""), "lon2": float(row["lon2"]),
                })
            except (ValueError, KeyError):
                continue
    return edges


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    rng   = np.random.default_rng(RNG_SEED)
    nodes = _load_nodes()
    edges = _load_edges()

    # Deduplicated vertices
    seen, verts_dedup = set(), []
    for n in nodes:
        key = (round(n["lon"], 4), round(n["lat"], 4))
        if key not in seen:
            seen.add(key)
            verts_dedup.append(n)

    n_edges    = len(edges)
    n_verts    = len(verts_dedup)

    mid_lons   = np.array([(e["lon1"] + e["lon2"]) / 2 for e in edges])
    v_lons     = np.array([v["lon"] for v in verts_dedup])

    mid_devs   = _beru_devs_np(mid_lons)
    v_devs     = _beru_devs_np(v_lons)

    # Degree map
    degree = Counter()
    for e in edges:
        degree[e["node1"]] += 1
        degree[e["node2"]] += 1

    # Degree-weighted vertex list (each city repeated by its degree)
    name_to_lon = {n["name"]: n["lon"] for n in nodes}
    dw_lons = np.array([
        lon for name, lon in name_to_lon.items()
        for _ in range(max(degree.get(name, 0), 1))
    ])
    n_dw = len(dw_lons)
    dw_devs = _beru_devs_np(dw_lons)

    # Observed counts
    obs_v_app  = int((v_devs   <= TIER_APP).sum())
    obs_v_ap   = int((v_devs   <= TIER_APLUS).sum())
    obs_v_a    = int((v_devs   <= TIER_A_MAX).sum())
    obs_m_app  = int((mid_devs <= TIER_APP).sum())
    obs_m_ap   = int((mid_devs <= TIER_APLUS).sum())
    obs_m_a    = int((mid_devs <= TIER_A_MAX).sum())
    obs_dw_app = int((dw_devs  <= TIER_APP).sum())
    obs_dw_ap  = int((dw_devs  <= TIER_APLUS).sum())
    obs_dw_a   = int((dw_devs  <= TIER_A_MAX).sum())

    # Binomial p-values
    def _bp(obs, n, null): return binomtest(obs, n, null, alternative="greater").pvalue

    p_v_app  = _bp(obs_v_app,  n_verts, P_NULL_APP)
    p_v_ap   = _bp(obs_v_ap,   n_verts, P_NULL_AP)
    p_v_a    = _bp(obs_v_a,    n_verts, P_NULL_A)
    p_m_app  = _bp(obs_m_app,  n_edges, P_NULL_APP)
    p_m_ap   = _bp(obs_m_ap,   n_edges, P_NULL_AP)
    p_m_a    = _bp(obs_m_a,    n_edges, P_NULL_A)
    p_dw_app = _bp(obs_dw_app, n_dw,    P_NULL_APP)
    p_dw_ap  = _bp(obs_dw_ap,  n_dw,    P_NULL_AP)
    p_dw_a   = _bp(obs_dw_a,   n_dw,    P_NULL_A)

    # Bonferroni-adjusted (k=4: 2 populations × 2 tiers)
    p_v_app_adj  = min(p_v_app  * BONF_K, 1.0)
    p_v_ap_adj   = min(p_v_ap   * BONF_K, 1.0)
    p_m_app_adj  = min(p_m_app  * BONF_K, 1.0)
    p_m_ap_adj   = min(p_m_ap   * BONF_K, 1.0)
    p_dw_app_adj = min(p_dw_app * BONF_K, 1.0)
    p_dw_ap_adj  = min(p_dw_ap  * BONF_K, 1.0)

    # Cluster asymmetry (A++-bearing harmonics vs non-A++ harmonics)
    harm_deg  = defaultdict(list)
    harm_app  = defaultdict(bool)
    for name, lon in name_to_lon.items():
        bv   = (lon - GERIZIM) / BERU
        hidx = round(bv / HARMONIC_STEP)
        dev  = abs(bv - hidx * HARMONIC_STEP)
        deg  = max(degree.get(name, 0), 1)
        harm_deg[hidx].append(deg)
        if dev <= TIER_APP:
            harm_app[hidx] = True

    app_totals    = [sum(harm_deg[h]) for h in harm_deg if harm_app[h]]
    nonapp_totals = [sum(harm_deg[h]) for h in harm_deg if not harm_app[h]]
    n_app_harm    = len(app_totals)
    n_nonapp_harm = len(nonapp_totals)
    mean_app_deg    = float(np.mean(app_totals))  if app_totals    else 0
    mean_nonapp_deg = float(np.mean(nonapp_totals)) if nonapp_totals else 0
    clust_ratio     = mean_app_deg / max(mean_nonapp_deg, 0.01)
    _, p_clust_deg  = mannwhitneyu(app_totals, nonapp_totals, alternative="greater")

    # Phase-shift permutation (vertices)
    print("  Running phase-shift permutation…", file=sys.stderr)
    v_bv     = (v_lons - GERIZIM) / BERU
    offsets  = rng.uniform(0, HARMONIC_STEP, size=N_PERMS)
    BATCH    = 5_000
    phase_app_null = np.zeros(N_PERMS, dtype=np.int32)
    phase_ap_null  = np.zeros(N_PERMS, dtype=np.int32)
    for start in range(0, N_PERMS, BATCH):
        end    = min(start + BATCH, N_PERMS)
        shifted = v_bv[None, :] + offsets[start:end, None]
        near   = np.round(shifted / HARMONIC_STEP) * HARMONIC_STEP
        devs   = np.abs(shifted - near)
        phase_app_null[start:end] = (devs <= TIER_APP).sum(axis=1)
        phase_ap_null[start:end]  = (devs <= TIER_APLUS).sum(axis=1)
    p_phase_app = float((phase_app_null >= obs_v_app).mean())
    p_phase_ap  = float((phase_ap_null  >= obs_v_ap).mean())

    # Node-label permutation (edge midpoints)
    print("  Running node-label permutation…", file=sys.stderr)
    pool_names = list({e["node1"] for e in edges} | {e["node2"] for e in edges})
    pool_lons  = np.array([name_to_lon.get(n, 0.0) for n in pool_names])
    name_idx   = {n: i for i, n in enumerate(pool_names)}
    idx1 = np.array([name_idx.get(e["node1"], 0) for e in edges])
    idx2 = np.array([name_idx.get(e["node2"], 0) for e in edges])

    perm_app = np.zeros(N_PERMS, dtype=np.int32)
    perm_ap  = np.zeros(N_PERMS, dtype=np.int32)
    for i in range(N_PERMS):
        shuf     = pool_lons[rng.permutation(len(pool_lons))]
        mid_p    = (shuf[idx1] + shuf[idx2]) / 2
        devs_p   = _beru_devs_np(mid_p)
        perm_app[i] = (devs_p <= TIER_APP).sum()
        perm_ap[i]  = (devs_p <= TIER_APLUS).sum()
    p_perm_app = float((perm_app >= obs_m_app).mean())
    p_perm_ap  = float((perm_ap  >= obs_m_ap).mean())

    # A++ and A+ vertex list for the audit
    v_tier_rows = []
    for v, dev in zip(verts_dedup, v_devs):
        t = _tier(dev)
        if t in ("A++", "A+"):
            v_tier_rows.append({**v, "dev": float(dev),
                                 "dev_km": float(dev * BERU * 111),
                                 "tier": t,
                                 "degree": degree.get(v["name"], 0)})
    v_tier_rows.sort(key=lambda r: r["dev"])

    # ── Write output ──────────────────────────────────────────────────────────
    now = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")
    W   = 100

    lines = []
    def L(s=""): lines.append(s)
    def HR(): lines.append("─" * W)
    def SIG(p): return f"{p:.4f}  {_sig(p)}"

    L("OWTRAD SILK ROAD HARMONIC ALIGNMENT AUDIT")
    L(f"Generated : {now}")
    L(f"Script    : tools/generate_audit_owtrad.py")
    L(f"Source    : Ciolek, T.M. (2004+). OWTRAD Project. www.ciolek.com/OWTRAD/")
    L(f"Data      : data/store/silk_road/owtrad_routes.csv  ({n_edges} edges)")
    L(f"          : data/store/silk_road/owtrad_nodes.csv   ({len(nodes)} raw nodes, {n_verts} deduplicated)")
    L(f"Anchor    : Gerizim  {GERIZIM}°E  |  BERU = {BERU}°")
    L(f"Thresholds: A++ ≤ {TIER_APP} beru ({TIER_APP*BERU*111:.1f} km)  |  "
      f"A+ ≤ {TIER_APLUS} beru ({TIER_APLUS*BERU*111:.1f} km)  |  "
      f"A ≤ {TIER_A_MAX} beru ({TIER_A_MAX*BERU*111:.1f} km)")
    L(f"Null rates: A++ {P_NULL_APP:.4f}  |  A+ {P_NULL_AP:.4f}  |  A {P_NULL_A:.4f}  (geometric, uniform longitude)")
    L(f"Perms     : {N_PERMS:,}  (seed {RNG_SEED})  |  Bonferroni k = {BONF_K}")
    HR()

    L("POPULATION SIZES")
    L()
    L(f"  Deduplicated vertices          : {n_verts}")
    L(f"  Route edges                    : {n_edges}")
    L(f"  Degree-weighted vertex list    : {n_dw}  (each city repeated by its edge count)")
    L()
    HR()

    L("BINOMIAL TEST RESULTS  (one-tailed, H0: uniform longitude)")
    L()

    def pop_block(label, n, obs_app, obs_ap, obs_a,
                  p_app, p_ap, p_a, p_app_adj, p_ap_adj):
        L(f"  Population: {label}  (N = {n})")
        L(f"  {'Tier':<6}  {'Obs':>5}  {'Rate':>6}  {'Null':>6}  {'Enrich':>7}  {'Z':>6}  {'p (raw)':>10}  {'sig':>3}  {'p (Bonf k=4)':>14}  {'sig':>3}")
        L("  " + "─" * 86)
        for tier_nm, obs, null in [
            ("A++", obs_app, P_NULL_APP),
            ("A+",  obs_ap,  P_NULL_AP),
            ("A",   obs_a,   P_NULL_A),
        ]:
            rate   = obs / n
            enrich = rate / null
            z      = _zscore(obs, n, null)
            p_raw  = _bp(obs, n, null)
            p_adj  = min(p_raw * BONF_K, 1.0) if tier_nm in ("A++", "A+") else None
            adj_str = f"{p_adj:.4f}  {_sig(p_adj)}" if p_adj is not None else "  —"
            L(f"  {tier_nm:<6}  {obs:>5}  {rate:>5.1%}  {null:>5.1%}  {enrich:>6.2f}x  {z:>6.2f}  {p_raw:>10.4f}  {_sig(p_raw):>3}  {adj_str:>17}")
        L()

    pop_block("Deduplicated vertices", n_verts,
              obs_v_app, obs_v_ap, obs_v_a,
              p_v_app, p_v_ap, p_v_a, p_v_app_adj, p_v_ap_adj)
    pop_block("Route edge midpoints", n_edges,
              obs_m_app, obs_m_ap, obs_m_a,
              p_m_app, p_m_ap, p_m_a, p_m_app_adj, p_m_ap_adj)
    pop_block("Degree-weighted vertices", n_dw,
              obs_dw_app, obs_dw_ap, obs_dw_a,
              p_dw_app, p_dw_ap, p_dw_a, p_dw_app_adj, p_dw_ap_adj)
    HR()

    L("PHASE-SHIFT PERMUTATION TEST  (vertices, N_perms = {:,})".format(N_PERMS))
    L()
    L("  Vertex beru-values shifted by a uniform random offset in [0, HARMONIC_STEP).")
    L("  Tests whether the Gerizim-anchored 3° grid specifically drives enrichment.")
    L()
    L(f"  Tier   Observed   Null mean ± SD          p          sig")
    L("  " + "─" * 60)
    for tier_nm, obs, null_arr in [("A++", obs_v_app, phase_app_null),
                                    ("A+",  obs_v_ap,  phase_ap_null)]:
        mu  = float(null_arr.mean())
        sd  = float(null_arr.std(ddof=1))
        p   = float((null_arr >= obs).mean())
        L(f"  {tier_nm:<6} {obs:>8}   {mu:>6.2f} ± {sd:>5.2f}          {p:.4f}  {_sig(p)}")
    L()
    HR()

    L("NODE-LABEL PERMUTATION TEST  (edge midpoints, N_perms = {:,})".format(N_PERMS))
    L()
    L("  City longitudes kept fixed; edge assignments randomly shuffled.")
    L("  Tests whether routing choices (not city positions) drive the signal.")
    L()
    L(f"  Tier   Observed   Null mean ± SD          p          sig")
    L("  " + "─" * 60)
    for tier_nm, obs, null_arr in [("A++", obs_m_app, perm_app),
                                    ("A+",  obs_m_ap,  perm_ap)]:
        mu  = float(null_arr.mean())
        sd  = float(null_arr.std(ddof=1))
        p   = float((null_arr >= obs).mean())
        L(f"  {tier_nm:<6} {obs:>8}   {mu:>6.2f} ± {sd:>5.2f}          {p:.4f}  {_sig(p)}")
    L()
    HR()

    L("CLUSTER ASYMMETRY TEST  (degree-weighted, Mann-Whitney one-tailed)")
    L()
    L("  Harmonics bearing at least one A++ vertex vs harmonics without one.")
    L("  Tests whether A++-bearing harmonics are more heavily trafficked.")
    L()
    L(f"  A++ harmonics   : {n_app_harm}   mean connectivity = {mean_app_deg:.1f} edges/harmonic")
    L(f"  Non-A++ harmonics: {n_nonapp_harm}   mean connectivity = {mean_nonapp_deg:.1f} edges/harmonic")
    L(f"  Ratio           : {clust_ratio:.2f}x")
    L(f"  Mann-Whitney p  : {p_clust_deg:.4f}  {_sig(p_clust_deg)}")
    L()
    HR()

    L("A++ AND A+ DEDUPLICATED VERTICES  (sorted by deviation, tightest first)")
    L()
    L(f"  {'Tier':<5}  {'Name':<35}  {'Lon°E':>8}  {'Lat°N':>7}  {'Dev (beru)':>11}  {'Dev (km)':>9}  {'Deg':>4}  Country")
    L("  " + "─" * 100)
    for r in v_tier_rows:
        L(f"  {r['tier']:<5}  {r['name']:<35}  {r['lon']:>8.4f}  {r['lat']:>7.4f}  "
          f"{r['dev']:>11.6f}  {r['dev_km']:>9.2f}  {r['degree']:>4}  {r['country']}")
    L()
    L(f"  Total A++ vertices: {sum(1 for r in v_tier_rows if r['tier']=='A++')}")
    L(f"  Total A+  vertices: {sum(1 for r in v_tier_rows if r['tier']=='A+')}")
    HR()

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Written: {OUT_PATH}", file=sys.stderr)


if __name__ == "__main__":
    main()
