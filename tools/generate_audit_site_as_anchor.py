"""
generate_audit_site_as_anchor.py
=================================
Produce supplementary/audit/site_as_anchor_audit.txt

Use every site in the extended corpus (1011 inscribed Cultural/Mixed UNESCO
sites + Gerizim synthetic) as its own anchor.  For each site, count A+, A++,
and A hits against the remaining N = 1011 sites (self-excluded), then rank
all 1012 sites by A+ count.

This directly answers: "Is the ~35°E corridor anomalous when measured by how
well each inscribed site serves as an anchor for the rest of the corpus?"

Run from repo root:
    python3 tools/generate_audit_site_as_anchor.py
"""

import sys
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent))
from data.unesco_corpus import (
    load_corpus, cultural_sites_with_coords_extended,
)
from lib.beru import GERIZIM, BERU, TIER_APP, TIER_APLUS, TIER_A_MAX, P_NULL_AP, CONFIG

OUT = Path(__file__).parent.parent / "supplementary" / "audit" / "site_as_anchor_audit.txt"
SEP = "─" * 110

JERUSALEM_LON = CONFIG["anchors"]["jerusalem"]["longitude"]   # 35.2317°E

# ── Build extended corpus ─────────────────────────────────────────────────────
corpus      = load_corpus()
sites_ext   = cultural_sites_with_coords_extended(corpus)   # 1012 sites
lons_all    = np.array([s.longitude for s in sites_ext])
N_EXT       = len(sites_ext)    # 1012
N_FOCAL     = N_EXT - 1         # 1011 (working set after self-exclusion)


def count_tiers(anchor_lon: float) -> tuple[int, int, int]:
    """Return (n_app, n_ap, n_a) for this anchor against the rest of the corpus."""
    mask = np.abs(lons_all - anchor_lon) > 0.001
    arcs = np.abs(lons_all[mask] - anchor_lon)
    arcs = np.minimum(arcs, 360 - arcs)
    bvs  = arcs / BERU
    devs = np.abs(bvs - np.round(bvs * 10) / 10)
    return (int(np.sum(devs <= TIER_APP)),
            int(np.sum(devs <= TIER_APLUS)),
            int(np.sum(devs <= TIER_A_MAX)))


print("Computing A+ count for each of 1012 sites as anchor…", flush=True)
results = []
for s in sites_ext:
    app, ap, a = count_tiers(s.longitude)
    results.append({
        "name":    s.site,
        "lon":     s.longitude,
        "app":     app,
        "ap":      ap,
        "a":       a,
        "is_ger":  abs(s.longitude - GERIZIM) < 0.001,
        "is_jer":  abs(s.longitude - JERUSALEM_LON) < 0.001,
    })
print("Done.", flush=True)

# Sort by A+ descending, then A++ descending, then longitude
results.sort(key=lambda x: (-x["ap"], -x["app"], x["lon"]))

# Assign ranks (tied sites share the same rank)
ap_vals = np.array([r["ap"] for r in results])
for i, r in enumerate(results):
    r["rank"] = int(np.sum(ap_vals > r["ap"])) + 1

# Stats
ap_mean  = float(np.mean(ap_vals))
ap_std   = float(np.std(ap_vals))
ap_max   = int(np.max(ap_vals))
ap_min   = int(np.min(ap_vals))
n_sites  = len(results)

# Focal anchors
ger_row  = next(r for r in results if r["is_ger"])
jer_row  = next(r for r in results if r["is_jer"])

def sig(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"

p_ger = binomtest(ger_row["ap"], N_FOCAL, P_NULL_AP, alternative="greater").pvalue
p_jer = binomtest(jer_row["ap"], N_FOCAL, P_NULL_AP, alternative="greater").pvalue

ts = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")

lines = [
    "SITE-AS-ANCHOR RANKING AUDIT",
    f"Generated : {ts}",
    f"Script    : tools/generate_audit_site_as_anchor.py",
    f"Corpus    : {N_EXT} sites (1011 inscribed Cultural/Mixed + Gerizim synthetic)",
    f"           Each site used as anchor; tests against remaining N = {N_FOCAL} sites",
    f"Threshold : A+ ≤ {TIER_APLUS} beru ({TIER_APLUS*BERU*111:.1f} km)",
    f"Null rate : {P_NULL_AP:.0%} geometric (4% of 360° in A+ bands)",
    "",
]

# ── Distribution ──────────────────────────────────────────────────────────────
lines += [
    SEP,
    "DISTRIBUTION SUMMARY",
    "",
    f"  Sites ranked:  {n_sites}",
    f"  Mean A+:       {ap_mean:.2f} ± {ap_std:.2f}",
    f"  Max A+:        {ap_max}",
    f"  Min A+:        {ap_min}",
    f"  Expected (4%): {0.04 * N_FOCAL:.1f}",
    "",
]

# ── Focal anchors ─────────────────────────────────────────────────────────────
pct_above_ger = 100.0 * sum(1 for r in results if r["ap"] > ger_row["ap"]) / n_sites
pct_above_jer = 100.0 * sum(1 for r in results if r["ap"] > jer_row["ap"]) / n_sites

lines += [
    SEP,
    "FOCAL ANCHOR RESULTS",
    "",
    f"  {'Site':<45}  {'Lon°E':>7}  {'A++':>4}  {'A+':>4}  {'A':>5}  "
    f"{'Rank':>5}  {'% above':>8}  {'p-value':>10}  sig",
    "  " + "-" * 100,
    f"  {'Mount Gerizim (Tentative List)':<45}  {ger_row['lon']:>7.3f}  "
    f"{ger_row['app']:>4}  {ger_row['ap']:>4}  {ger_row['a']:>5}  "
    f"{ger_row['rank']:>5}  {pct_above_ger:>7.1f}%  {p_ger:>10.4e}  {sig(p_ger)}",
    f"  {'Old City of Jerusalem':<45}  {jer_row['lon']:>7.4f}  "
    f"{jer_row['app']:>4}  {jer_row['ap']:>4}  {jer_row['a']:>5}  "
    f"{jer_row['rank']:>5}  {pct_above_jer:>7.1f}%  {p_jer:>10.4e}  {sig(p_jer)}",
    "",
]

# ── Full ranked listing ───────────────────────────────────────────────────────
lines += [
    SEP,
    f"FULL RANKING  (all {n_sites} sites, sorted by A+ count)",
    "",
    f"  {'Rank':>5}  {'A++':>4}  {'A+':>4}  {'A':>5}  {'Lon°E':>9}  Site",
    "  " + "-" * 105,
]

prev_ap = None
for r in results:
    # Insert blank line between A+ tiers for readability
    if prev_ap is not None and r["ap"] != prev_ap:
        lines.append("")
    prev_ap = r["ap"]

    flag = ""
    if r["is_ger"]:
        flag = "  ◄ GERIZIM (synthetic)"
    elif r["is_jer"]:
        flag = "  ◄ JERUSALEM"

    lines.append(
        f"  {r['rank']:>5}  {r['app']:>4}  {r['ap']:>4}  {r['a']:>5}  "
        f"{r['lon']:>9.4f}  {r['name'][:65]}{flag}"
    )

lines.append("")

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text("\n".join(lines), encoding="utf-8")
print(f"Written -> {OUT}")
print(f"  {n_sites} sites ranked")
print(f"  Gerizim  : A+={ger_row['ap']}  rank {ger_row['rank']}/{n_sites}  ({pct_above_ger:.1f}% above)")
print(f"  Jerusalem: A+={jer_row['ap']}  rank {jer_row['rank']}/{n_sites}  ({pct_above_jer:.1f}% above)")
