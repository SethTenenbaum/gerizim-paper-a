#!/usr/bin/env python3
"""
owtrad_route_alignment.py — Harmonic-tier analysis of OWTRAD Silk Road
edge midpoints with three levels of robustness testing.

Tests
-----
  Binomial       — observed A++/A+ midpoint counts vs. uniform null
  Z-score        — standardised effect size from binomial mean/SD
  Node-label     — keep the 194 actual city longitudes, randomly re-pair them
                   as edges; tests whether the *routing choices* are harmonic
                   rather than the city distribution
  Phase-shift    — shift the Gerizim anchor by random offsets in [0, 3°);
                   tests whether the Gerizim-anchored grid specifically drives
                   the signal vs. any regular 3° grid

Writes to ResultsStore (and emits \\newcommand lines)
-----------------------------------------------------
  Counts / sizes
    NOwtradEdges, NOwtradVertices
    NOwtradMidApp, NOwtradMidAp
    NOwtradVertexApp, NOwtradVertexAp
    NOwtradPerms
  Enrichment ratios
    owtradMidAppEnrich, owtradMidApEnrich
  Z-scores
    zOwtradMidApp, zOwtradMidAp
  Binomial p-values (uncorrected + Bonferroni k=4)
    pOwtradMidApp,    pOwtradMidAp
    pOwtradMidAppAdj, pOwtradMidApAdj
    pOwtradVertexApp, pOwtradVertexAp
  Node-label permutation p-values
    pOwtradMidAppPerm, pOwtradMidApPerm
  Phase-shift permutation p-values
    pOwtradMidAppPhase, pOwtradMidApPhase
"""

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from lib.beru import (
    GERIZIM, BERU, HARMONIC_STEP,
    TIER_APP, TIER_APLUS,
    P_NULL_APP, P_NULL_AP,
    deviation as beru_dev,
)
from lib.results_store import ResultsStore

SILK_ROAD_DIR = ROOT / "data" / "store" / "silk_road"
BONF_K    = 4         # 2 populations × 2 tiers (planned comparisons)
N_PERMS   = 100_000
RNG_SEED  = 42


# ── Data loading ───────────────────────────────────────────────────────────────

def _load_edges_raw():
    """Return list of (node1_name, lon1, node2_name, lon2)."""
    rows = []
    with open(SILK_ROAD_DIR / "owtrad_routes.csv", newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            try:
                rows.append((
                    row.get("node1", ""),  float(row["lon1"]),
                    row.get("node2", ""),  float(row["lon2"]),
                ))
            except (ValueError, KeyError):
                continue
    return rows


def _load_vertices_dedup():
    """Return deduplicated list of node longitudes from owtrad_nodes.csv."""
    lons, seen = [], set()
    with open(SILK_ROAD_DIR / "owtrad_nodes.csv", newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            try:
                lon, lat = float(row["lon"]), float(row["lat"])
            except (ValueError, KeyError):
                continue
            key = (round(lon, 4), round(lat, 4))
            if key not in seen:
                seen.add(key)
                lons.append(lon)
    return lons


# ── Helpers ────────────────────────────────────────────────────────────────────

def _beru_devs_np(lons_arr, anchor=GERIZIM):
    """Vectorised beru deviation for a numpy array of longitudes."""
    bv   = (lons_arr - anchor) / BERU
    near = np.round(bv / HARMONIC_STEP) * HARMONIC_STEP
    return np.abs(bv - near)


def _zscore(obs, n, p_null):
    mu  = n * p_null
    sd  = np.sqrt(n * p_null * (1 - p_null))
    return (obs - mu) / sd if sd > 0 else 0.0


# ── Main ───────────────────────────────────────────────────────────────────────

def run():
    rng = np.random.default_rng(RNG_SEED)

    edges_raw    = _load_edges_raw()
    vertex_lons  = _load_vertices_dedup()

    n_edges   = len(edges_raw)
    n_vertices = len(vertex_lons)

    mid_lons  = np.array([(r[1] + r[3]) / 2 for r in edges_raw])
    v_lons    = np.array(vertex_lons)

    # ── Observed counts ────────────────────────────────────────────────────────
    mid_devs  = _beru_devs_np(mid_lons)
    v_devs    = _beru_devs_np(v_lons)

    obs_m_app = int((mid_devs <= TIER_APP).sum())
    obs_m_ap  = int((mid_devs <= TIER_APLUS).sum())
    obs_v_app = int((v_devs   <= TIER_APP).sum())
    obs_v_ap  = int((v_devs   <= TIER_APLUS).sum())

    # ── Binomial tests ─────────────────────────────────────────────────────────
    p_m_app = binomtest(obs_m_app, n_edges, P_NULL_APP, alternative="greater").pvalue
    p_m_ap  = binomtest(obs_m_ap,  n_edges, P_NULL_AP,  alternative="greater").pvalue
    p_v_app = binomtest(obs_v_app, n_vertices, P_NULL_APP, alternative="greater").pvalue
    p_v_ap  = binomtest(obs_v_ap,  n_vertices, P_NULL_AP,  alternative="greater").pvalue

    p_m_app_adj = min(p_m_app * BONF_K, 1.0)
    p_m_ap_adj  = min(p_m_ap  * BONF_K, 1.0)

    # ── Z-scores ───────────────────────────────────────────────────────────────
    z_m_app = _zscore(obs_m_app, n_edges, P_NULL_APP)
    z_m_ap  = _zscore(obs_m_ap,  n_edges, P_NULL_AP)

    # ── Enrichment ratios ──────────────────────────────────────────────────────
    enrich_m_app = (obs_m_app / n_edges) / P_NULL_APP
    enrich_m_ap  = (obs_m_ap  / n_edges) / P_NULL_AP

    # ── Node-label permutation ─────────────────────────────────────────────────
    # Build node-name → index map using coordinates from the edge file directly,
    # so that every node referenced by an edge is present in the pool.
    node_name_to_lon = {}
    for n1, l1, n2, l2 in edges_raw:
        if n1 not in node_name_to_lon:
            node_name_to_lon[n1] = l1
        if n2 not in node_name_to_lon:
            node_name_to_lon[n2] = l2

    pool_names = list(node_name_to_lon.keys())
    pool_lons  = np.array([node_name_to_lon[n] for n in pool_names])
    name_to_idx = {n: i for i, n in enumerate(pool_names)}

    edge_idx1 = np.array([name_to_idx[r[0]] for r in edges_raw])
    edge_idx2 = np.array([name_to_idx[r[2]] for r in edges_raw])

    print(f"  Node-label permutation: pool={len(pool_lons)} nodes, "
          f"{n_edges} edges, {N_PERMS:,} perms…", file=sys.stderr)

    perm_app = np.zeros(N_PERMS, dtype=np.int32)
    perm_ap  = np.zeros(N_PERMS, dtype=np.int32)

    for i in range(N_PERMS):
        shuffled        = pool_lons[rng.permutation(len(pool_lons))]
        mid_perm        = (shuffled[edge_idx1] + shuffled[edge_idx2]) / 2
        devs_perm       = _beru_devs_np(mid_perm)
        perm_app[i]     = (devs_perm <= TIER_APP).sum()
        perm_ap[i]      = (devs_perm <= TIER_APLUS).sum()

    p_perm_app = float(np.mean(perm_app >= obs_m_app))
    p_perm_ap  = float(np.mean(perm_ap  >= obs_m_ap))

    # ── Phase-shift permutation ────────────────────────────────────────────────
    # Shift all midpoint beru-values by a uniform random offset in [0, HARMONIC_STEP).
    # Equivalent to testing a random regular 3° grid instead of Gerizim's grid.
    print(f"  Phase-shift permutation: {N_PERMS:,} perms…", file=sys.stderr)

    beru_vals   = (mid_lons - GERIZIM) / BERU                   # shape (n_edges,)
    offsets     = rng.uniform(0, HARMONIC_STEP, size=N_PERMS)   # shape (N_PERMS,)

    # Vectorised: (N_PERMS, n_edges)
    BATCH = 5_000
    phase_app = np.zeros(N_PERMS, dtype=np.int32)
    phase_ap  = np.zeros(N_PERMS, dtype=np.int32)
    for start in range(0, N_PERMS, BATCH):
        end  = min(start + BATCH, N_PERMS)
        sh   = beru_vals[np.newaxis, :] + offsets[start:end, np.newaxis]
        near = np.round(sh / HARMONIC_STEP) * HARMONIC_STEP
        devs = np.abs(sh - near)
        phase_app[start:end] = (devs <= TIER_APP).sum(axis=1)
        phase_ap[start:end]  = (devs <= TIER_APLUS).sum(axis=1)

    p_phase_app = float(np.mean(phase_app >= obs_m_app))
    p_phase_ap  = float(np.mean(phase_ap  >= obs_m_ap))

    # ── Store & emit ───────────────────────────────────────────────────────────
    store = ResultsStore()

    results = {
        # Counts
        "NOwtradEdges":        n_edges,
        "NOwtradVertices":     n_vertices,
        "NOwtradMidApp":       obs_m_app,
        "NOwtradMidAp":        obs_m_ap,
        "NOwtradVertexApp":    obs_v_app,
        "NOwtradVertexAp":     obs_v_ap,
        "NOwtradPerms":        N_PERMS,
        # Enrichment
        "owtradMidAppEnrich":  round(enrich_m_app, 2),
        "owtradMidApEnrich":   round(enrich_m_ap,  2),
        # Z-scores
        "zOwtradMidApp":       round(z_m_app, 2),
        "zOwtradMidAp":        round(z_m_ap,  2),
        # Binomial
        "pOwtradMidApp":       p_m_app,
        "pOwtradMidAp":        p_m_ap,
        "pOwtradMidAppAdj":    p_m_app_adj,
        "pOwtradMidApAdj":     p_m_ap_adj,
        "pOwtradVertexApp":    p_v_app,
        "pOwtradVertexAp":     p_v_ap,
        # Node-label permutation
        "pOwtradMidAppPerm":   p_perm_app,
        "pOwtradMidApPerm":    p_perm_ap,
        # Phase-shift permutation
        "pOwtradMidAppPhase":  p_phase_app,
        "pOwtradMidApPhase":   p_phase_ap,
    }

    for key, val in results.items():
        store.write(key, val)

    def _fmt(v):
        if isinstance(v, int):
            return str(v)
        if v >= 100:
            return f"{v:.0f}"
        if v >= 1:
            return f"{v:.2f}"
        if v < 0.0001:
            return "< 0.0001"
        return f"{v:.4f}"

    for name, val in results.items():
        print(f"\\newcommand{{\\{name}}}{{{_fmt(val)}}}")

    # Summary to stderr
    print(f"\n  ── Results ──────────────────────────────────", file=sys.stderr)
    print(f"  Midpoints A++: {obs_m_app}/{n_edges}  "
          f"enrich={enrich_m_app:.2f}×  z={z_m_app:.2f}  "
          f"binom p={p_m_app:.4f}  perm p={p_perm_app:.4f}  phase p={p_phase_app:.4f}",
          file=sys.stderr)
    print(f"  Midpoints A+:  {obs_m_ap}/{n_edges}  "
          f"enrich={enrich_m_ap:.2f}×  z={z_m_ap:.2f}  "
          f"binom p={p_m_ap:.4f}  perm p={p_perm_ap:.4f}  phase p={p_phase_ap:.4f}",
          file=sys.stderr)
    print(f"  Vertices  A++: {obs_v_app}/{n_vertices}  "
          f"binom p={p_v_app:.4f}", file=sys.stderr)
    print(f"  Vertices  A+:  {obs_v_ap}/{n_vertices}  "
          f"binom p={p_v_ap:.4f}", file=sys.stderr)


if __name__ == "__main__":
    run()
