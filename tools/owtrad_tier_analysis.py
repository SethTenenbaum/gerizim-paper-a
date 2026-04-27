#!/usr/bin/env python3
"""
owtrad_tier_analysis.py — Harmonic-tier analysis of the OWTRAD Silk Road network.

Tests whether OWTRAD network vertices (cities/waypoints) and edge midpoints
align with beru harmonics at rates exceeding the uniform null expectation,
using the same tier definitions and binomial tests as the main UNESCO corpus
analysis.

Populations analysed
--------------------
  vertices    — 194 OWTRAD nodes (deduplicated by lon/lat), each tested by
                its own longitude.
  midpoints   — 251 edge midpoint longitudes, i.e. (lon1+lon2)/2 for each
                attested route segment.

Output
------
  Console table + saved text report: data/store/silk_road/owtrad_tier_report.txt
  Figure: manuscript/figures/fig_owtrad_tiers.pdf / .png

Usage
-----
  cd gerizim-paper-a
  python3 tools/owtrad_tier_analysis.py
"""

import csv
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import binomtest

from lib.beru import (
    GERIZIM, BERU, HARMONIC_STEP,
    TIER_APP, TIER_APLUS, TIER_A_MAX,
    P_NULL_APP, P_NULL_AP, P_NULL_A,
    TIER_APP_KM, TIER_APLUS_KM, TIER_A_KM,
    deviation as beru_dev, tier_label as tier_label,
)
from lib.stats import significance_label as sig_label

SILK_ROAD_DIR = PROJECT_ROOT / "data" / "store" / "silk_road"
OUTDIR        = PROJECT_ROOT / "manuscript" / "figures"
REPORT_PATH   = SILK_ROAD_DIR / "owtrad_tier_report.txt"

plt.rcParams.update({
    "font.family": "serif", "font.size": 10,
    "axes.labelsize": 11, "axes.titlesize": 12,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 9, "figure.dpi": 300,
    "savefig.dpi": 300, "savefig.bbox": "tight",
    "savefig.pad_inches": 0.08,
    "axes.spines.top": False, "axes.spines.right": False,
})

C_PRIMARY   = "#2c6fbb"
C_HIGHLIGHT = "#d62728"
C_NULL      = "#888888"
C_APP       = "#922b21"
C_AP        = "#c0392b"
C_A         = "#e67e22"


# ── Data loading ──────────────────────────────────────────────────────────────

def load_nodes():
    """Return list of dicts for deduplicated OWTRAD nodes."""
    rows = []
    seen = set()
    with open(SILK_ROAD_DIR / "owtrad_nodes.csv", newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            try:
                lon = float(row["lon"])
                lat = float(row["lat"])
            except (ValueError, KeyError):
                continue
            key = (round(lon, 4), round(lat, 4))
            if key in seen:
                continue
            seen.add(key)
            rows.append({
                "name": row.get("name", ""),
                "lon":  lon, "lat": lat,
                "dataset": row.get("dataset", ""),
                "country": row.get("country", ""),
                "dev":  beru_dev(lon),
                "tier": tier_label(beru_dev(lon)),
            })
    return rows


def load_edges():
    """Return list of dicts for OWTRAD edges, including midpoint longitude."""
    rows = []
    with open(SILK_ROAD_DIR / "owtrad_routes.csv", newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            try:
                lon1 = float(row["lon1"]); lat1 = float(row["lat1"])
                lon2 = float(row["lon2"]); lat2 = float(row["lat2"])
            except (ValueError, KeyError):
                continue
            mid_lon = (lon1 + lon2) / 2
            mid_lat = (lat1 + lat2) / 2
            rows.append({
                "node1": row.get("node1", ""), "node2": row.get("node2", ""),
                "lon1": lon1, "lat1": lat1, "lon2": lon2, "lat2": lat2,
                "mid_lon": mid_lon, "mid_lat": mid_lat,
                "dataset": row.get("dataset", ""),
                "dev":  beru_dev(mid_lon),
                "tier": tier_label(beru_dev(mid_lon)),
                "dev1": beru_dev(lon1), "tier1": tier_label(beru_dev(lon1)),
                "dev2": beru_dev(lon2), "tier2": tier_label(beru_dev(lon2)),
            })
    return rows


# ── Statistical summary ───────────────────────────────────────────────────────

def tier_stats(lons, label):
    """Return dict of counts and binomial p-values for a set of longitudes."""
    n = len(lons)
    devs = [beru_dev(lon) for lon in lons]

    n_app  = sum(1 for d in devs if d <= TIER_APP)
    n_ap   = sum(1 for d in devs if d <= TIER_APLUS)
    n_a    = sum(1 for d in devs if d <= TIER_A_MAX)

    p_app = binomtest(n_app, n, P_NULL_APP, alternative="greater").pvalue
    p_ap  = binomtest(n_ap,  n, P_NULL_AP,  alternative="greater").pvalue
    p_a   = binomtest(n_a,   n, P_NULL_A,   alternative="greater").pvalue

    return {
        "label": label, "n": n,
        "n_app": n_app, "rate_app": 100*n_app/n if n else 0,
        "null_app": 100*P_NULL_APP, "p_app": p_app,
        "n_ap":  n_ap,  "rate_ap":  100*n_ap/n  if n else 0,
        "null_ap":  100*P_NULL_AP,  "p_ap":  p_ap,
        "n_a":   n_a,   "rate_a":   100*n_a/n   if n else 0,
        "null_a":   100*P_NULL_A,   "p_a":   p_a,
        "devs": devs,
    }


def print_stats(s, file=None):
    def _w(line):
        print(line)
        if file:
            file.write(line + "\n")

    _w(f"\n{'─'*60}")
    _w(f"  Population: {s['label']}  (N = {s['n']})")
    _w(f"{'─'*60}")
    _w(f"  {'Tier':<6} {'Obs':>5} {'Rate':>7}  {'Null':>7}  {'p':>9}  Sig")
    _w(f"  {'─'*55}")
    for tier, k, rate, null, p in [
        ("A++",  s['n_app'], s['rate_app'], s['null_app'], s['p_app']),
        ("A+",   s['n_ap'],  s['rate_ap'],  s['null_ap'],  s['p_ap']),
        ("A",    s['n_a'],   s['rate_a'],   s['null_a'],   s['p_a']),
    ]:
        _w(f"  {tier:<6} {k:>5}  {rate:>6.1f}%  {null:>6.1f}%  {p:>9.4f}  {sig_label(p)}")


# ── Figure ────────────────────────────────────────────────────────────────────

def make_figure(vertex_stats, midpoint_stats, all_node_lons, all_mid_lons):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    def _panel(ax, s, lons, title):
        devs = np.array(s["devs"])
        bin_width = 0.005
        max_dev = 0.05
        bins = np.arange(0, max_dev + bin_width, bin_width)
        n_bins = len(bins) - 1

        counts, edges, patches = ax.hist(
            devs, bins=bins,
            color=C_PRIMARY, edgecolor="white", linewidth=0.6,
            alpha=0.85, zorder=3,
        )

        ax.axvspan(0, TIER_APP,   color=C_APP, alpha=0.10, zorder=1)
        ax.axvspan(0, TIER_APLUS, color=C_AP,  alpha=0.08, zorder=1)
        ax.axvspan(0, TIER_A_MAX, color=C_A,   alpha=0.05, zorder=1)

        expected = s["n"] / n_bins
        ax.axhline(expected, color=C_NULL, linewidth=1.5, linestyle="--", zorder=4)
        # Label the null line at the far-right edge of the plot
        ax.text(max_dev * 0.99, expected * 1.04,
                f"null ({expected:.1f}/bin)",
                color=C_NULL, fontsize=7.5, ha="right", va="bottom", zorder=5)

        ax.axvline(TIER_APP,   color=C_APP, linewidth=0.8, linestyle=":", alpha=0.7, zorder=4)
        ax.axvline(TIER_APLUS, color=C_AP,  linewidth=0.8, linestyle=":", alpha=0.7, zorder=4)
        ax.axvline(TIER_A_MAX, color=C_A,   linewidth=0.8, linestyle=":", alpha=0.7, zorder=4)

        y_max = counts.max()
        stats_txt = (
            f"A++: {s['n_app']}/{s['n']}  {sig_label(s['p_app'])}\n"
            f"A+:  {s['n_ap']}/{s['n']}  {sig_label(s['p_ap'])}\n"
            f"A:   {s['n_a']}/{s['n']}  {sig_label(s['p_a'])}"
        )
        ax.text(0.97, 0.97, stats_txt,
                transform=ax.transAxes, fontsize=8.5,
                ha="right", va="top", family="monospace",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f0f0",
                          edgecolor="#cccccc", alpha=0.9))

        ax.set_xlabel("Beru deviation (δ) from nearest 0.1-beru harmonic")
        ax.set_ylabel("Count")
        ax.set_title(title)
        ax.set_xlim(0, max_dev)
        ax.set_ylim(0, y_max * 1.25)

    _panel(axes[0], vertex_stats,   all_node_lons,
           f"OWTRAD vertices (N = {vertex_stats['n']}, deduplicated)")
    _panel(axes[1], midpoint_stats, all_mid_lons,
           f"OWTRAD vertices degree-weighted (Σdeg = {midpoint_stats['n']})")

    fig.suptitle(
        "OWTRAD Silk Road network — beru deviation distribution\n"
        "Left: each city once.  Right: each city weighted by its route-count.\n"
        "Shaded: A++ / A+ / A windows.  Dashed: uniform null expectation.",
        fontsize=10, y=1.02,
    )
    fig.tight_layout()
    out = OUTDIR / "fig_owtrad_tiers.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    print(f"\n  ✓ {out.name} ({out.with_suffix('.png').name})")
    return out


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("  OWTRAD Silk Road — harmonic tier analysis")
    print("=" * 60)

    # Load data
    nodes = load_nodes()
    edges = load_edges()
    print(f"\n  Vertices (deduplicated):  {len(nodes)}")
    print(f"  Edges:                    {len(edges)}")

    # Per-dataset vertex breakdown
    from collections import Counter
    ds_counts = Counter(n["dataset"] for n in nodes)
    print(f"\n  Nodes by dataset:")
    for ds, cnt in sorted(ds_counts.items()):
        print(f"    {ds:<20} {cnt}")

    # Tier distributions
    tier_counts_nodes = Counter(n["tier"] for n in nodes)
    tier_counts_edges = Counter(e["tier"] for e in edges)
    print(f"\n  Vertex tier distribution:")
    for t in ("A++","A+","A","B","C","C-","C--"):
        c = tier_counts_nodes.get(t, 0)
        print(f"    {t:<5}  {c:>3}  ({100*c/len(nodes):.1f}%)")
    print(f"\n  Edge midpoint tier distribution:")
    for t in ("A++","A+","A","B","C","C-","C--"):
        c = tier_counts_edges.get(t, 0)
        print(f"    {t:<5}  {c:>3}  ({100*c/len(edges):.1f}%)")

    # Statistical tests
    node_lons = [n["lon"] for n in nodes]
    mid_lons  = [e["mid_lon"] for e in edges]

    vstats = tier_stats(node_lons, "OWTRAD vertices (deduplicated)")
    mstats = tier_stats(mid_lons,  "OWTRAD edge midpoints")

    # ── Vertex degree analysis ────────────────────────────────────────────────
    from collections import Counter as _Counter
    degree = _Counter()
    for e in edges:
        degree[e["node1"]] += 1
        degree[e["node2"]] += 1

    for n in nodes:
        n["degree"] = degree.get(n["name"], 0)

    # Degree-weighted: count each vertex proportional to its network degree
    # (equivalent to sampling edges uniformly and picking an endpoint at random)
    weighted_lons = []
    for n in nodes:
        weighted_lons.extend([n["lon"]] * n["degree"])
    wstats = tier_stats(weighted_lons, f"Vertices degree-weighted (Σdeg = {len(weighted_lons)})")

    # Hub vs. non-hub split (degree ≥ 3 = junction city)
    hub_threshold = 3
    hub_nodes  = [n for n in nodes if n["degree"] >= hub_threshold]
    leaf_nodes = [n for n in nodes if n["degree"] <  hub_threshold]
    hstats = tier_stats([n["lon"] for n in hub_nodes],
                        f"Hub vertices (degree ≥ {hub_threshold}, N = {len(hub_nodes)})")
    lstats = tier_stats([n["lon"] for n in leaf_nodes],
                        f"Leaf/non-hub vertices (degree < {hub_threshold}, N = {len(leaf_nodes)})")

    print(f"\n  Degree distribution:")
    deg_vals = sorted(degree.values(), reverse=True)
    print(f"    max={deg_vals[0]}  median={sorted(deg_vals)[len(deg_vals)//2]}  "
          f"mean={sum(deg_vals)/len(deg_vals):.1f}  min={deg_vals[-1]}")
    print(f"    Hub nodes (deg ≥ {hub_threshold}): {len(hub_nodes)}")
    print(f"    Leaf nodes (deg <  {hub_threshold}): {len(leaf_nodes)}")
    print(f"\n  Top hub nodes:")
    for n in sorted(nodes, key=lambda x: -x["degree"])[:12]:
        print(f"    deg={n['degree']}  {n['tier']:<4}  δ={n['dev']:.4f}"
              f"  {n['lon']:>8.3f}°E  {n['name']}")

    all_pop_stats = (vstats, mstats, wstats, hstats, lstats)

    with open(REPORT_PATH, "w", encoding="utf-8") as fh:
        fh.write("OWTRAD Silk Road — harmonic tier analysis\n")
        fh.write(f"Generated by tools/owtrad_tier_analysis.py\n\n")
        fh.write(f"Null rates (from lib/beru.py):\n")
        fh.write(f"  A++: {P_NULL_APP:.4f} ({100*P_NULL_APP:.2f}%)  [{TIER_APP_KM:.1f} km window]\n")
        fh.write(f"  A+:  {P_NULL_AP:.4f}  ({100*P_NULL_AP:.2f}%)  [{TIER_APLUS_KM:.1f} km window]\n")
        fh.write(f"  A:   {P_NULL_A:.4f}  ({100*P_NULL_A:.2f}%)  [{TIER_A_KM:.0f} km window]\n\n")
        fh.write("Binomial test: one-tailed (greater), H0 = uniform null\n")
        for s in all_pop_stats:
            print_stats(s, file=fh)

    print(f"\n  Report saved → {REPORT_PATH}")

    for s in all_pop_stats:
        print_stats(s)

    # A+ city listing
    print(f"\n{'─'*60}")
    print("  A+ and A++ OWTRAD vertices (all):")
    print(f"{'─'*60}")
    ap_nodes = sorted(
        [n for n in nodes if n["tier"] in ("A+", "A++")],
        key=lambda x: x["dev"],
    )
    for n in ap_nodes:
        print(f"  {n['tier']:<4}  δ={n['dev']:.4f}  deg={n['degree']}  "
              f"{n['lon']:>9.3f}°E  {n['name']}  ({n['country']})")

    print(f"\n{'─'*60}")
    print("  A+ and A++ edge midpoints:")
    print(f"{'─'*60}")
    ap_edges = sorted(
        [e for e in edges if e["tier"] in ("A+", "A++")],
        key=lambda x: x["dev"],
    )
    for e in ap_edges:
        print(f"  {e['tier']:<4}  δ={e['dev']:.4f}  mid={e['mid_lon']:>9.3f}°E"
              f"  {e['node1']} → {e['node2']}")

    # Figure: vertices (raw) vs degree-weighted vertices
    make_figure(vstats, wstats, node_lons, weighted_lons)

    print()
    print("=" * 60)
    print("  Done.")
    print("=" * 60)
