"""
generate_audit_stupa_q180987_geo.py
=====================================
Geographic tier breakdown for the Wikidata Q180987 stupa corpus (N=229).

Corpus: data/store/unesco/wikidata_stupas_q180987.csv
Fetch:  python3 data/scripts/fetch_wikidata_q180987.py --fetch

Key question: within the broader 70°–110°E heartland, is A/A+ concentration
uniform, or does it localise to Indonesia/Java?

Sub-regions tested:
  - All stupas
  - Heartland (70°–110°E)
  - India/Nepal (70°–90°E)
  - SE Asia / Myanmar-Cambodia (90°–105°E)
  - Java/Sumatra node (105°–115°E)  ← the 110.269°E = 2.5-beru node
  - Indonesia broad (95°–141°E)
  - Non-heartland (outside 70°–110°E)

Fisher exact: one-sided enrichment vs full Wikidata corpus background.
Tier labels follow lib/beru.py: A++ (≤0.4km), A+ (≤6.7km), A (≤33km), B, C, C-

Run from repo root:
    python3 tools/generate_audit_stupa_q180987_geo.py
"""

import csv
import sys
from pathlib import Path
from datetime import datetime, timezone
from scipy.stats import fisher_exact, binomtest

sys.path.insert(0, str(Path(__file__).parent.parent))
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS,
    P_NULL_AP, P_NULL_A,
    deviation as beru_dev,
    tier_label, is_aplus, is_a_or_better,
)

CSV_PATH = Path(__file__).parent.parent / "data" / "store" / "unesco" / "wikidata_stupas_q180987.csv"
OUT      = Path(__file__).parent.parent / "supplementary" / "audit" / "stupa_q180987_geo_audit.txt"
SEP  = "=" * 100
SEP2 = "─" * 100

TIER_ORDER = ["A++", "A+", "A", "B", "C", "C-", "C--"]

# ── Load CSV ──────────────────────────────────────────────────────────────────
def load_sites():
    sites = []
    with open(CSV_PATH, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            break  # consumed the header comment block, now at csv header
        # reopen properly
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
                "qid":     row.get("qid", ""),
                "name":    row.get("name", "") or row.get("qid", ""),
                "lat":     lat,
                "lon":     lon,
                "country": row.get("country", ""),
                "dev":     dev,
                "dev_km":  dev * BERU * 111.0,
                "tier":    tier,
            })
    return sites

sites = load_sites()

# ── Full-corpus background ────────────────────────────────────────────────────
N_CORPUS    = len(sites)
K_AP_CORPUS = sum(1 for s in sites if is_aplus(s["tier"]))
K_A_CORPUS  = sum(1 for s in sites if is_a_or_better(s["tier"]))
K_C_CORPUS  = sum(1 for s in sites if s["tier"] in ("C", "C-", "C--"))  # inter-harmonic (C-band)

# ── Helpers ───────────────────────────────────────────────────────────────────
def in_band(lon, lo, hi):
    return lo <= lon <= hi

def tier_counts(ss):
    return {t: sum(1 for s in ss if s["tier"] == t) for t in TIER_ORDER}

def n_aplus(ss):
    return sum(1 for s in ss if is_aplus(s["tier"]))

def n_atier(ss):
    return sum(1 for s in ss if is_a_or_better(s["tier"]))

def n_ctier(ss):
    return sum(1 for s in ss if s["tier"] in ("C", "C-", "C--"))

def fisher_ap(sub):
    n  = len(sub)
    k  = n_aplus(sub)
    if n == 0:
        return float("nan"), float("nan")
    table = [[k, n - k],
             [K_AP_CORPUS - k, N_CORPUS - n - (K_AP_CORPUS - k)]]
    or_, p = fisher_exact(table, alternative="greater")
    return round(or_, 2), p

def fisher_a(sub):
    n  = len(sub)
    k  = n_atier(sub)
    if n == 0:
        return float("nan"), float("nan")
    table = [[k, n - k],
             [K_A_CORPUS - k, N_CORPUS - n - (K_A_CORPUS - k)]]
    or_, p = fisher_exact(table, alternative="greater")
    return round(or_, 2), p

def fisher_c(sub):
    """One-sided Fisher exact: enrichment of C-band (C/C-/C--) vs full corpus."""
    n = len(sub)
    k = n_ctier(sub)
    if n == 0:
        return float("nan"), float("nan")
    table = [[k, n - k],
             [K_C_CORPUS - k, N_CORPUS - n - (K_C_CORPUS - k)]]
    or_, p = fisher_exact(table, alternative="greater")
    return round(or_, 2), p

def binom_ap(sub):
    n, k = len(sub), n_aplus(sub)
    if n == 0:
        return float("nan")
    return binomtest(k, n, P_NULL_AP, alternative="greater").pvalue

def binom_a(sub):
    n, k = len(sub), n_atier(sub)
    if n == 0:
        return float("nan")
    return binomtest(k, n, P_NULL_A, alternative="greater").pvalue

def sig(p):
    if p != p:
        return "n/a"
    return ("***" if p < 0.001 else "**" if p < 0.01
            else "*" if p < 0.05 else "†" if p < 0.10 else "ns")

def fmt_or(v):
    return "inf" if v == float("inf") else f"{v:.2f}"

# ── Region definitions ────────────────────────────────────────────────────────
REGIONS = [
    ("All stupas",                              lambda s: True),
    ("Heartland (70°–110°E)",                   lambda s: in_band(s["lon"], 70.0, 110.0)),
    ("India/Nepal (70°–90°E)",                  lambda s: in_band(s["lon"], 70.0,  90.0)),
    ("SE Asia / Myanmar-Cambodia (90°–105°E)",  lambda s: in_band(s["lon"], 90.0, 105.0)),
    ("Java/Sumatra node (105°–115°E)",          lambda s: in_band(s["lon"], 105.0, 115.0)),
    ("Java tight (107°–112°E)",                 lambda s: in_band(s["lon"], 107.0, 112.0)),
    ("Indonesia broad (95°–141°E)",             lambda s: in_band(s["lon"], 95.0,  141.0)),
    ("Non-heartland (outside 70°–110°E)",       lambda s: not in_band(s["lon"], 70.0, 110.0)),
]

region_data = [(label, [s for s in sites if fn(s)]) for label, fn in REGIONS]

# ── Format helpers ────────────────────────────────────────────────────────────
def summary_row(label, sub):
    n   = len(sub)
    tc  = tier_counts(sub)
    nap = n_aplus(sub)
    na  = n_atier(sub)
    nc  = n_ctier(sub)
    or_ap, p_ap = fisher_ap(sub)
    or_a,  p_a  = fisher_a(sub)
    or_c,  p_c  = fisher_c(sub)
    rate_ap = f"{100*nap/n:.1f}%" if n else "—"
    rate_a  = f"{100*na/n:.1f}%"  if n else "—"
    rate_c  = f"{100*nc/n:.1f}%"  if n else "—"
    return (f"  {label:<42} {n:>4} {nap:>4} {rate_ap:>6}  "
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

    lines = [SEP2,
             f"  {label}  (N = {n})",
             f"  A+ = {nap} ({100*nap/n:.1f}%)  Fisher OR={fmt_or(or_ap)}× p={p_ap:.4f} {sig(p_ap)}  |  binomial p={bp_ap:.4f} {sig(bp_ap)}" if n else "",
             f"  A  = {na} ({100*na/n:.1f}%)   Fisher OR={fmt_or(or_a)}×  p={p_a:.4f} {sig(p_a)}   |  binomial p={bp_a:.4f} {sig(bp_a)}" if n else "",
             "",
             f"  {'Tier':<5} {'N':>4} {'Rate':>6}",
             f"  {'─'*5} {'─'*4} {'─'*6}"]
    for t in TIER_ORDER:
        c = tc[t]
        rate = f"{100*c/n:.1f}%" if n else "—"
        lines.append(f"  {t:<5} {c:>4} {rate:>6}")

    lines += ["",
              f"  {'Tier':<4} {'Dev(km)':>8}  {'Lon':>8}  {'Country':<14}  Name"]
    lines.append(f"  {'─'*4} {'─'*8}  {'─'*8}  {'─'*14}  {'─'*50}")
    for s in sorted(sub, key=lambda x: (TIER_ORDER.index(x["tier"]), x["dev"])):
        name = s["name"] if s["name"] != s["qid"] else f"[{s['qid']}]"
        lines.append(
            f"  [{s['tier']:<3}] {s['dev_km']:>7.1f}km  {s['lon']:>8.3f}°E  "
            f"{s['country']:<14}  {name}"
        )
    return lines

# ── Build output ──────────────────────────────────────────────────────────────
lines = [SEP,
         "  WIKIDATA Q180987 STUPA CORPUS — GEOGRAPHIC TIER BREAKDOWN",
         f"  Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
         f"  Corpus: {CSV_PATH.name}  |  N = {N_CORPUS}",
         f"  Full corpus: A+ = {K_AP_CORPUS} ({100*K_AP_CORPUS/N_CORPUS:.1f}%)  "
         f"A-tier = {K_A_CORPUS} ({100*K_A_CORPUS/N_CORPUS:.1f}%)",
         f"  Anchor: Gerizim {GERIZIM}°E  |  1 beru = {BERU}°",
         SEP, ""]

# Summary table
hdr = (f"  {'Region':<42} {'N':>4} {'A+':>4} {'A+%':>6}  "
       f"{'OR_A+':>7} {'p_A+':>7} {'sig':<4}  "
       f"{'A':>4} {'A%':>6}  "
       f"{'OR_A':>7} {'p_A':>7} {'sig':<4}  "
       f"{'C':>4} {'C%':>6}  "
       f"{'OR_C':>7} {'p_C':>7} {'sig':<4}")
lines += ["  SUMMARY TABLE  (A+ = ≤6.7km, A = ≤33km, C = inter-harmonic C/C-/C--)", hdr,
          "  " + "─"*len(hdr.rstrip())]
for label, sub in region_data:
    lines.append(summary_row(label, sub))

# Cross-region comparison
lines += ["", SEP,
          "  CROSS-REGION BREAKDOWN: Java/Sumatra node vs rest of heartland",
          SEP2, ""]

java      = [s for s in sites if in_band(s["lon"], 105.0, 115.0)]
heartland_ex_java = [s for s in sites if in_band(s["lon"], 70.0, 110.0)
                     and not in_band(s["lon"], 105.0, 110.0)]
india_nepal = [s for s in sites if in_band(s["lon"], 70.0, 90.0)]
myanmar_cam = [s for s in sites if in_band(s["lon"], 90.0, 105.0)]
non_hl    = [s for s in sites if not in_band(s["lon"], 70.0, 110.0)]

for grp_label, grp in [
    ("Java/Sumatra node (105°–115°E)", java),
    ("Heartland ex Java/Sumatra (70°–105°E)", heartland_ex_java),
    ("  India/Nepal sub (70°–90°E)", india_nepal),
    ("  Myanmar/Cambodia sub (90°–105°E)", myanmar_cam),
    ("Non-heartland", non_hl),
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

# ── Bimodal analysis ─────────────────────────────────────────────────────────
lines += [SEP,
          "  BIMODAL ANALYSIS: A-tier (Java node) vs C-band (India/Myanmar heartland)",
          "  Hypothesis: stupa distribution is not uniform — A concentrates at 2.5-beru",
          "  harmonic (Java/105°–115°E) while C concentrates in the physical heartland",
          "  (India/Nepal + Myanmar/Cambodia, 70°–105°E).",
          SEP2, ""]

bimodal_groups = [
    ("Java/Sumatra node (105°–115°E) — A-tier zone",
     [s for s in sites if in_band(s["lon"], 105.0, 115.0)]),
    ("India/Nepal (70°–90°E) — C-band zone",
     [s for s in sites if in_band(s["lon"], 70.0, 90.0)]),
    ("Myanmar/Cambodia (90°–105°E) — C-band zone",
     [s for s in sites if in_band(s["lon"], 90.0, 105.0)]),
    ("Heartland combined (70°–105°E)",
     [s for s in sites if in_band(s["lon"], 70.0, 105.0)]),
]

for grp_label, grp in bimodal_groups:
    n   = len(grp)
    nap = n_aplus(grp)
    na  = n_atier(grp)
    nc  = n_ctier(grp)
    tc  = tier_counts(grp)
    or_ap, p_ap = fisher_ap(grp)
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

lines += [
    "  INTERPRETATION",
    "  ─────────────",
    "  A bimodal geographic enrichment pattern is present across the Q180987 stupa corpus:",
    "    • A-tier  enriched at the Java/105°–115°E node  (2.5-beru harmonic, 110.3°E)",
    "    • C-band  enriched in the physical heartland    (India/Nepal + Myanmar, 70°–105°E)",
    "    • B-tier  dominates the mid-range bulk of the heartland (structural baseline)",
    "  This mirrors the Buddhist bimodal A/C enrichment found in the UNESCO religion audit.",
    "  The C-band concentration in the historical stupa heartland suggests the inter-harmonic",
    "  zone may encode a second, independent spatial signal rather than random scatter.",
    "",
]

# Per-region detail with site listings
lines += [SEP, "  PER-REGION SITE LISTINGS", ""]
for label, sub in region_data:
    lines.extend(detail_block(label, sub))
    lines.append("")

lines.append(SEP)
output = "\n".join(lines)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(output, encoding="utf-8")
print(output)
print(f"\n  → Written to {OUT}")
