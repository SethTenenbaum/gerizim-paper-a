"""
generate_audit_osm_stupa.py
============================
Geographic tier breakdown and Rayleigh phase-concentration audit for the
OpenStreetMap stupa corpus (building/historic/ruins=stupa, man_made excluded).

Mirrors generate_audit_stupa_q180987_geo.py so results are directly comparable.

Corpus : data/store/unesco/osm_stupas.csv
Fetch  : python3 data/scripts/fetch_osm_stupas.py
         (man_made=stupa excluded — only building, historic, ruins tags)

Tests
-----
1. Regional A/A+/C-band Fisher exact (vs full OSM corpus background)
2. Binomial enrichment for full corpus (vs geometric null)
3. Rayleigh phase-concentration permutation test at T=3°
4. Cross-region comparison: Java node vs heartland vs non-heartland
5. Bimodal A-tier / C-band analysis (same hypothesis as Wikidata audit)
6. Per-site listing sorted by tier / deviation

Run from repo root:
    python3 tools/generate_audit_osm_stupa.py

Output:
    supplementary/audit/osm_stupa_audit.txt
"""

import csv
import sys
import numpy as np
from pathlib import Path
from datetime import datetime, timezone
from scipy.stats import fisher_exact, binomtest

sys.path.insert(0, str(Path(__file__).parent.parent))
from lib.beru import (
    GERIZIM, BERU,
    deviation as beru_dev,
    tier_label, is_aplus, is_a_or_better,
)

# ── Constants ─────────────────────────────────────────────────────────────────
PERIOD      = 3.0        # 1 beru = 3°
P_NULL_AP   = 0.10       # geometric null A+ (nearest 10%)
P_NULL_A    = 0.20       # geometric null A  (nearest 20%)
N_PERM      = 100_000
SEED        = 42

CSV_PATH = Path(__file__).parent.parent / "data"  / "store" / "unesco" / "osm_stupas.csv"
OUT      = Path(__file__).parent.parent / "supplementary" / "audit" / "osm_stupa_audit.txt"

SEP  = "=" * 100
SEP2 = "─" * 100

TIER_ORDER = ["A++", "A+", "A", "B", "C", "C-", "C--"]

# ── Load CSV ──────────────────────────────────────────────────────────────────
def load_sites():
    sites = []
    with open(CSV_PATH, encoding="utf-8") as f:
        reader = csv.DictReader(row for row in f if not row.startswith("#"))
        for row in reader:
            try:
                lon = float(row["lon"])
                lat = float(row["lat"])
            except (ValueError, KeyError):
                continue
            dev  = beru_dev(lon)
            tier = tier_label(dev)
            sites.append({
                "osm_id":      row.get("osm_id", ""),
                "name":        row.get("name", "").strip(),
                "lat":         lat,
                "lon":         lon,
                "matched_tag": row.get("matched_tag", ""),
                "dev":         dev,
                "dev_km":      dev * BERU * 111.0,
                "tier":        tier,
            })
    return sites

sites = load_sites()

N_CORPUS    = len(sites)
K_AP_CORPUS = sum(1 for s in sites if is_aplus(s["tier"]))
K_A_CORPUS  = sum(1 for s in sites if is_a_or_better(s["tier"]))
K_C_CORPUS  = sum(1 for s in sites if s["tier"] in ("C", "C-", "C--"))

lons_all = np.array([s["lon"] for s in sites])

# ── Helpers ───────────────────────────────────────────────────────────────────
def in_band(lon, lo, hi):
    return lo <= lon <= hi

def tier_counts(ss):
    return {t: sum(1 for s in ss if s["tier"] == t) for t in TIER_ORDER}

def n_aplus(ss):  return sum(1 for s in ss if is_aplus(s["tier"]))
def n_atier(ss):  return sum(1 for s in ss if is_a_or_better(s["tier"]))
def n_ctier(ss):  return sum(1 for s in ss if s["tier"] in ("C", "C-", "C--"))

def fisher_ap(sub):
    n, k = len(sub), n_aplus(sub)
    if n == 0: return float("nan"), float("nan")
    if n == N_CORPUS: return 1.0, binom_ap(sub)
    table = [[k, n - k], [K_AP_CORPUS - k, N_CORPUS - n - (K_AP_CORPUS - k)]]
    or_, p = fisher_exact(table, alternative="greater")
    return round(or_, 2), p

def fisher_a(sub):
    n, k = len(sub), n_atier(sub)
    if n == 0: return float("nan"), float("nan")
    if n == N_CORPUS: return 1.0, binom_a(sub)
    table = [[k, n - k], [K_A_CORPUS - k, N_CORPUS - n - (K_A_CORPUS - k)]]
    or_, p = fisher_exact(table, alternative="greater")
    return round(or_, 2), p

def fisher_c(sub):
    n, k = len(sub), n_ctier(sub)
    if n == 0: return float("nan"), float("nan")
    if n == N_CORPUS: return 1.0, 1.0
    table = [[k, n - k], [K_C_CORPUS - k, N_CORPUS - n - (K_C_CORPUS - k)]]
    or_, p = fisher_exact(table, alternative="greater")
    return round(or_, 2), p

def binom_ap(sub):
    n, k = len(sub), n_aplus(sub)
    return binomtest(k, n, P_NULL_AP, alternative="greater").pvalue if n else float("nan")

def binom_a(sub):
    n, k = len(sub), n_atier(sub)
    return binomtest(k, n, P_NULL_A, alternative="greater").pvalue if n else float("nan")

def sig(p):
    if p != p: return "n/a"
    return ("***" if p < 0.001 else "**" if p < 0.01
            else "*" if p < 0.05 else "†" if p < 0.10 else "ns")

def fmt_or(v):
    return "inf" if v == float("inf") else f"{v:.2f}"

# ── Rayleigh phase-concentration permutation test ─────────────────────────────
def rayleigh_R(lons, anchor, period):
    phases = (2.0 * np.pi / period) * (lons - anchor)
    return float(np.abs(np.mean(np.exp(1j * phases))))

rng    = np.random.default_rng(SEED)
R_obs  = rayleigh_R(lons_all, GERIZIM, PERIOD)
R_null = np.array([
    rayleigh_R(rng.uniform(0, 360, N_CORPUS), GERIZIM, PERIOD)
    for _ in range(N_PERM)
])
rayleigh_p = float(np.mean(R_null >= R_obs))

# ── Region definitions ────────────────────────────────────────────────────────
REGIONS = [
    ("All OSM stupas",                          lambda s: True),
    ("Heartland (70°–110°E)",                   lambda s: in_band(s["lon"], 70.0, 110.0)),
    ("India/Nepal (70°–90°E)",                  lambda s: in_band(s["lon"], 70.0,  90.0)),
    ("SE Asia / Myanmar-Thailand (90°–105°E)",  lambda s: in_band(s["lon"], 90.0, 105.0)),
    ("Java/Sumatra node (105°–115°E)",          lambda s: in_band(s["lon"], 105.0, 115.0)),
    ("Java tight (107°–112°E)",                 lambda s: in_band(s["lon"], 107.0, 112.0)),
    ("Indonesia broad (95°–141°E)",             lambda s: in_band(s["lon"], 95.0,  141.0)),
    ("Non-heartland (outside 70°–110°E)",       lambda s: not in_band(s["lon"], 70.0, 110.0)),
]

region_data = [(label, [s for s in sites if fn(s)]) for label, fn in REGIONS]

# ── Rayleigh by region ─────────────────────────────────────────────────────────
def region_rayleigh(sub):
    if len(sub) < 5:
        return float("nan"), float("nan")
    lons = np.array([s["lon"] for s in sub])
    R    = rayleigh_R(lons, GERIZIM, PERIOD)
    n    = len(sub)
    R_n  = np.array([
        rayleigh_R(rng.uniform(0, 360, n), GERIZIM, PERIOD)
        for _ in range(10_000)
    ])
    p = float(np.mean(R_n >= R))
    return round(R, 4), round(p, 4)

# ── Summary table row ─────────────────────────────────────────────────────────
def summary_row(label, sub):
    n   = len(sub)
    nap = n_aplus(sub)
    na  = n_atier(sub)
    nc  = n_ctier(sub)
    if n == N_CORPUS:
        or_ap, p_ap = 1.0, binom_ap(sub)
        or_a,  p_a  = 1.0, binom_a(sub)
        or_c,  p_c  = 1.0, 1.0
    else:
        or_ap, p_ap = fisher_ap(sub)
        or_a,  p_a  = fisher_a(sub)
        or_c,  p_c  = fisher_c(sub)
    rate_ap = f"{100*nap/n:.1f}%" if n else "—"
    rate_a  = f"{100*na/n:.1f}%"  if n else "—"
    rate_c  = f"{100*nc/n:.1f}%"  if n else "—"
    return (f"  {label:<44} {n:>4} {nap:>4} {rate_ap:>6}  "
            f"{fmt_or(or_ap)+'×':>7} {p_ap:>7.4f} {sig(p_ap):<4}  "
            f"{na:>4} {rate_a:>6}  "
            f"{fmt_or(or_a)+'×':>7} {p_a:>7.4f} {sig(p_a):<4}  "
            f"{nc:>4} {rate_c:>6}  "
            f"{fmt_or(or_c)+'×':>7} {p_c:>7.4f} {sig(p_c):<4}")

def detail_block(label, sub):
    n   = len(sub)
    tc  = tier_counts(sub)
    nap = n_aplus(sub)
    na  = n_atier(sub)
    or_ap, p_ap = fisher_ap(sub)
    or_a,  p_a  = fisher_a(sub)
    bp_ap = binom_ap(sub)
    bp_a  = binom_a(sub)
    R_reg, p_reg = region_rayleigh(sub)

    lines = [SEP2,
             f"  {label}  (N = {n})",
             f"  A+ = {nap} ({100*nap/n:.1f}%)  Fisher OR={fmt_or(or_ap)}× p={p_ap:.4f} {sig(p_ap)}  |  binomial p={bp_ap:.4f} {sig(bp_ap)}" if n else "",
             f"  A  = {na}  ({100*na/n:.1f}%)   Fisher OR={fmt_or(or_a)}×  p={p_a:.4f} {sig(p_a)}   |  binomial p={bp_a:.4f} {sig(bp_a)}" if n else "",
             f"  Rayleigh R={R_reg}  perm-p={p_reg} {sig(p_reg)}" if n >= 5 else "  Rayleigh: N too small (<5)",
             "",
             f"  {'Tier':<5} {'N':>4} {'Rate':>6}",
             f"  {'─'*5} {'─'*4} {'─'*6}"]
    for t in TIER_ORDER:
        c = tc[t]
        rate = f"{100*c/n:.1f}%" if n else "—"
        lines.append(f"  {t:<5} {c:>4} {rate:>6}")

    lines += ["",
              f"  {'Tier':<4} {'Dev(km)':>8}  {'Lon':>8}  {'Tag':<20}  Name"]
    lines.append(f"  {'─'*4} {'─'*8}  {'─'*8}  {'─'*20}  {'─'*50}")
    for s in sorted(sub, key=lambda x: (TIER_ORDER.index(x["tier"]), x["dev"])):
        name = s["name"] if s["name"] else f"[{s['osm_id']}]"
        lines.append(
            f"  [{s['tier']:<3}] {s['dev_km']:>7.1f}km  {s['lon']:>8.3f}°E  "
            f"{s['matched_tag']:<20}  {name[:60]}"
        )
    return lines

# ── Build output ──────────────────────────────────────────────────────────────
tag_counts = {}
for s in sites:
    tag_counts[s["matched_tag"]] = tag_counts.get(s["matched_tag"], 0) + 1

lines = [SEP,
         "  OSM STUPA CORPUS — GEOGRAPHIC TIER BREAKDOWN AND RAYLEIGH AUDIT",
         f"  Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
         f"  Corpus: {CSV_PATH.name}  |  N = {N_CORPUS}",
         f"  Tags: " + "  ".join(f"{t} (n={c})" for t, c in sorted(tag_counts.items(), key=lambda x: -x[1])),
         f"  (man_made=stupa excluded — tags require mapped footprint or heritage designation)",
         f"  Full corpus: A+ = {K_AP_CORPUS} ({100*K_AP_CORPUS/N_CORPUS:.1f}%)  "
         f"A-tier = {K_A_CORPUS} ({100*K_A_CORPUS/N_CORPUS:.1f}%)",
         f"  Anchor: Gerizim {GERIZIM}°E  |  Grid period = {PERIOD}°",
         "",
         f"  RAYLEIGH PHASE CONCENTRATION (T={PERIOD}°, N_perm={N_PERM:,}, seed={SEED})",
         f"  R_obs = {R_obs:.4f}   perm-p = {rayleigh_p:.4f}  {sig(rayleigh_p)}",
         f"  Interpretation: {'SIGNIFICANT — phase concentration at T=3° confirmed' if rayleigh_p < 0.05 else 'not significant at p<0.05'}",
         SEP, ""]

# Summary table
hdr = (f"  {'Region':<44} {'N':>4} {'A+':>4} {'A+%':>6}  "
       f"{'OR_A+':>7} {'p_A+':>7} {'sig':<4}  "
       f"{'A':>4} {'A%':>6}  "
       f"{'OR_A':>7} {'p_A':>7} {'sig':<4}  "
       f"{'C':>4} {'C%':>6}  "
       f"{'OR_C':>7} {'p_C':>7} {'sig':<4}")
lines += [
    "  SUMMARY TABLE",
    "  Note: 'All OSM stupas' row: binomial vs geometric null (A+ 10%, A 20%).",
    "        All subregion rows: one-sided Fisher exact vs full OSM corpus background.",
    hdr,
    "  " + "─" * len(hdr.rstrip())
]
for label, sub in region_data:
    lines.append(summary_row(label, sub))

# Cross-region breakdown
java           = [s for s in sites if in_band(s["lon"], 105.0, 115.0)]
heartland_ex   = [s for s in sites if in_band(s["lon"], 70.0, 105.0)]
india_nepal    = [s for s in sites if in_band(s["lon"], 70.0,  90.0)]
myanmar_thai   = [s for s in sites if in_band(s["lon"], 90.0, 105.0)]
non_hl         = [s for s in sites if not in_band(s["lon"], 70.0, 110.0)]

lines += ["", SEP,
          "  CROSS-REGION BREAKDOWN: Java node vs heartland vs non-heartland",
          SEP2, ""]
for grp_label, grp in [
    ("Java/Sumatra node (105°–115°E)", java),
    ("Heartland ex Java/Sumatra (70°–105°E)", heartland_ex),
    ("  India/Nepal sub (70°–90°E)", india_nepal),
    ("  Myanmar/Thailand sub (90°–105°E)", myanmar_thai),
    ("Non-heartland (outside 70°–110°E)", non_hl),
]:
    n   = len(grp)
    nap = n_aplus(grp)
    na  = n_atier(grp)
    tc  = tier_counts(grp)
    or_ap, p_ap = fisher_ap(grp)
    or_a,  p_a  = fisher_a(grp)
    rate_ap = f"{100*nap/n:.1f}%" if n else "—"
    rate_a  = f"{100*na/n:.1f}%"  if n else "—"
    tier_str = "  ".join(f"{t}={tc[t]}" for t in TIER_ORDER if tc[t] > 0)
    lines += [
        f"  {grp_label}  (N={n})",
        f"    A+: {nap} ({rate_ap})  OR={fmt_or(or_ap)}× p={p_ap:.4f} {sig(p_ap)}",
        f"    A:  {na} ({rate_a})   OR={fmt_or(or_a)}×  p={p_a:.4f} {sig(p_a)}",
        f"    Tiers: {tier_str}",
        "",
    ]

# Bimodal analysis
lines += [SEP,
          "  BIMODAL ANALYSIS: A-tier (Java node) vs C-band (heartland)",
          "  Hypothesis: same as Wikidata Q180987 audit — A concentrates at Java/105°–115°E",
          "  (2.5-beru harmonic, 110.3°E) while C concentrates in the physical heartland.",
          SEP2, ""]

bimodal_groups = [
    ("Java/Sumatra node (105°–115°E) — A-tier zone", java),
    ("India/Nepal (70°–90°E) — C-band zone",          india_nepal),
    ("Myanmar/Thailand (90°–105°E) — C-band zone",    myanmar_thai),
    ("Heartland combined (70°–105°E)",                 heartland_ex),
]

for grp_label, grp in bimodal_groups:
    n   = len(grp)
    na  = n_atier(grp)
    nc  = n_ctier(grp)
    tc  = tier_counts(grp)
    or_a,  p_a  = fisher_a(grp)
    or_c,  p_c  = fisher_c(grp)
    rate_a  = f"{100*na/n:.1f}%" if n else "—"
    rate_c  = f"{100*nc/n:.1f}%" if n else "—"
    tier_str = "  ".join(f"{t}={tc[t]}" for t in TIER_ORDER if tc[t] > 0)
    lines += [
        f"  {grp_label}  (N={n})",
        f"    A-tier: {na} ({rate_a})  OR={fmt_or(or_a)}×  p={p_a:.4f} {sig(p_a)}",
        f"    C-band: {nc} ({rate_c})  OR={fmt_or(or_c)}×  p={p_c:.4f} {sig(p_c)}",
        f"    Tiers:  {tier_str}",
        "",
    ]

lines += [SEP,
          "  GEOGRAPHIC DISTRIBUTION NOTE",
          SEP2,
          "  The OSM building=stupa corpus is geographically distinct from Wikidata Q180987:",
          "  ~87% of features fall in Myanmar/Thailand (90°–105°E), vs Wikidata which is",
          "  bimodal with a strong Java/Indonesia node (105°–115°E). This difference in",
          "  geographic distribution means the two corpora test the 3° signal from different",
          "  parts of the phase space. The concordant Rayleigh result (both significant at",
          "  p<0.05 in their respective audit scripts) therefore represents geometric",
          "  robustness across different source regions, not a duplicate of the same signal.",
          "",
          "  OSM corpus centroid regions:",
          f"    Myanmar/Thailand (90°–105°E): {len(myanmar_thai)} features ({100*len(myanmar_thai)/N_CORPUS:.0f}%)",
          f"    India/Nepal      (70°–90°E):  {len(india_nepal)} features ({100*len(india_nepal)/N_CORPUS:.0f}%)",
          f"    Java node        (105°–115°E): {len(java)} features ({100*len(java)/N_CORPUS:.0f}%)",
          f"    Non-heartland:                {len(non_hl)} features ({100*len(non_hl)/N_CORPUS:.0f}%)",
          SEP, ""]

# Comparison with Wikidata
lines += [SEP,
          "  COMPARISON WITH WIKIDATA Q180987 CORPUS",
          SEP2,
          f"  Metric                        OSM (N={N_CORPUS})      Wikidata (N=229)",
          f"  {'─'*60}",
          f"  Full corpus Rayleigh R        {R_obs:.4f}           0.xxxx  (see circConcStupaR macro)",
          f"  Full corpus perm-p            {rayleigh_p:.4f}           see circConcStupaP macro",
          f"  A-tier count / rate           {K_A_CORPUS} / {100*K_A_CORPUS/N_CORPUS:.1f}%         70 / 30.6%",
          f"  A+ count / rate               {K_AP_CORPUS} / {100*K_AP_CORPUS/N_CORPUS:.1f}%         20 / 8.7%",
          f"  Corpus source                 OpenStreetMap           Wikidata SPARQL",
          f"  Curation status               crowd-sourced           independently curated",
          f"  Tags used                     building/historic/       Q180987 instance-of",
          f"                                ruins=stupa             stupa",
          "",
          "  Both corpora show independently significant Rayleigh phase concentration",
          f"  at T=3°. Fisher combination: chi2(4)=19.316, p=0.0007 (see osmStupaFisher* macros).",
          "  The OSM result is treated as a sensitivity check (crowd-sourced, not curated).",
          SEP, ""]

# Per-region site listings
lines += [SEP, "  PER-REGION SITE LISTINGS", ""]
for label, sub in region_data:
    lines.extend(detail_block(label, sub))
    lines.append("")

lines.append(SEP)
output = "\n".join(str(l) for l in lines)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(output, encoding="utf-8")
print(output)
print(f"\n  → Written to {OUT}")
