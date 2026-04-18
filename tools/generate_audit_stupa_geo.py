"""
generate_audit_stupa_geo.py
===========================
Produce supplementary/audit/stupa_geo_audit.txt

Geographic tier breakdown for stupa sites:
  - Full stupa population
  - Heartland band (70°–110°E, South/Southeast Asia)
  - Indonesia sub-region (95°–141°E, Indonesian archipelago)
  - Java sub-region (105°–112°E)
  - Non-heartland stupas

For each sub-region:
  - Tier breakdown table (A+, A, B, C, C-)
  - Fisher exact vs full corpus background (one-sided, enrichment)
  - Site listings sorted by tier then deviation

Key observation being probed:
  Within the heartland band, C/C- appears enriched overall, but
  the Indonesia/Java sub-region concentrates A and A+ sites.

Run from repo root:
    python3 tools/generate_audit_stupa_geo.py
"""

import re
import sys
import json
from pathlib import Path
from datetime import datetime, timezone
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU,
    P_NULL_AP, P_NULL_A,
    deviation as beru_dev,
    tier_label, is_aplus, is_a_or_better,
)

OUT = Path(__file__).parent.parent / "supplementary" / "audit" / "stupa_geo_audit.txt"
SEP  = "=" * 100
SEP2 = "─" * 100

# ── Keyword setup ─────────────────────────────────────────────────────────────
_KW_PATH = Path(__file__).parent.parent / "keywords.json"
with open(_KW_PATH) as f:
    _KW = json.load(f)

STUPA_KEYWORDS = _KW["mound_evolution"]["stupa"]
STUPA_RES = {kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
             for kw in STUPA_KEYWORDS}

def is_stupa(site) -> bool:
    return any(rx.search(site.full_text) for rx in STUPA_RES.values())

def matched_keywords(site) -> list:
    return [kw for kw, rx in STUPA_RES.items() if rx.search(site.full_text)]

# ── Geographic sub-regions ────────────────────────────────────────────────────
REGIONS = [
    ("All stupas",                        (None,  None)),
    ("Heartland (70°–110°E)",             (70.0,  110.0)),
    ("Non-heartland",                     "invert"),     # special: outside 70–110
    ("Indonesia (95°–141°E)",             (95.0,  141.0)),
    ("Java (105°–112°E)",                 (105.0, 112.0)),
    ("India/Nepal (70°–90°E)",            (70.0,   90.0)),
    ("SE Asia ex-Indonesia (95°–110°E)",  (95.0,  110.0)),
]

def in_band(lon, lo, hi):
    return lo <= lon <= hi

# ── Load corpus ───────────────────────────────────────────────────────────────
corpus   = load_corpus()
cultural = cultural_sites_with_coords(corpus)

# Full corpus background counts
N_CORPUS    = len(cultural)
K_AP_CORPUS = sum(1 for s in cultural if is_aplus(tier_label(beru_dev(s.longitude))))
K_A_CORPUS  = sum(1 for s in cultural if is_a_or_better(tier_label(beru_dev(s.longitude))))

# Stupa population
stupa_sites = []
for s in cultural:
    if not is_stupa(s):
        continue
    dev  = beru_dev(s.longitude)
    tier = tier_label(dev)
    stupa_sites.append({
        "name":     s.site,
        "lon":      s.longitude,
        "dev":      dev,
        "dev_km":   dev * BERU * 111.0,
        "tier":     tier,
        "keywords": matched_keywords(s),
    })

stupa_sites.sort(key=lambda x: x["dev"])

# ── Helpers ───────────────────────────────────────────────────────────────────
TIER_ORDER = ["A+", "A", "B", "C", "C-"]

def tier_counts(sites):
    return {t: sum(1 for s in sites if s["tier"] == t) for t in TIER_ORDER}

def fisher_ap(n_sub, n_ap_sub):
    """One-sided Fisher exact: enrichment of A+ in sub-pop vs full corpus."""
    if n_sub == 0:
        return float("nan"), float("nan")
    table = [
        [n_ap_sub,              n_sub - n_ap_sub],
        [K_AP_CORPUS - n_ap_sub, N_CORPUS - n_sub - (K_AP_CORPUS - n_ap_sub)],
    ]
    or_, p = fisher_exact(table, alternative="greater")
    return round(or_, 2), p

def fisher_a(n_sub, n_a_sub):
    """One-sided Fisher exact: enrichment of A-tier in sub-pop vs full corpus."""
    if n_sub == 0:
        return float("nan"), float("nan")
    table = [
        [n_a_sub,              n_sub - n_a_sub],
        [K_A_CORPUS - n_a_sub, N_CORPUS - n_sub - (K_A_CORPUS - n_a_sub)],
    ]
    or_, p = fisher_exact(table, alternative="greater")
    return round(or_, 2), p

def sig(p):
    if p != p:  # nan
        return "n/a"
    return ("***" if p < 0.001 else "**" if p < 0.01
            else "*" if p < 0.05 else "†" if p < 0.10 else "ns")

def get_region_sites(label, band):
    """Return sites for a named region."""
    if band == "invert":
        return [s for s in stupa_sites if not in_band(s["lon"], 70.0, 110.0)]
    lo, hi = band
    if lo is None:
        return list(stupa_sites)
    return [s for s in stupa_sites if in_band(s["lon"], lo, hi)]

def tier_block(sites, label):
    """Format a tier-breakdown block for a region."""
    n   = len(sites)
    tc  = tier_counts(sites)
    n_ap = tc["A+"]
    n_a  = tc["A+"] + tc["A"]   # A-or-better
    or_ap, p_ap = fisher_ap(n, n_ap)
    or_a,  p_a  = fisher_a(n, n_a)
    lines = []
    lines.append(f"  {label}  (N = {n})")
    lines.append(f"  {'Tier':<6} {'Count':>5} {'Rate':>7}  {'vs corpus'}")
    lines.append(f"  {'─'*6} {'─'*5} {'─'*7}  {'─'*30}")
    for t in TIER_ORDER:
        c = tc[t]
        rate = f"{100*c/n:.1f}%" if n else "—"
        lines.append(f"  {t:<6} {c:>5} {rate:>7}")
    lines.append(f"  {'─'*50}")
    lines.append(f"  A+ Fisher exact: OR = {or_ap}×  p = {p_ap:.4f}  {sig(p_ap)}")
    lines.append(f"  A   Fisher exact: OR = {or_a}×  p = {p_a:.4f}  {sig(p_a)}")
    return lines

def site_listing(sites):
    lines = []
    lines.append(f"  {'Tier':<4} {'Dev(km)':>8}  {'Lon':>8}  Keywords / Name")
    lines.append(f"  {'─'*4} {'─'*8}  {'─'*8}  {'─'*60}")
    for s in sorted(sites, key=lambda x: (TIER_ORDER.index(x["tier"]), x["dev"])):
        kws = ", ".join(s["keywords"])
        lines.append(
            f"  [{s['tier']:<3}] {s['dev_km']:>7.1f}km  {s['lon']:>8.3f}°E  "
            f"({kws})  {s['name']}"
        )
    return lines

# ── Build region list ─────────────────────────────────────────────────────────
region_data = []
for label, spec in REGIONS:
    sites = get_region_sites(label, spec)
    region_data.append((label, sites))

# ── Write output ──────────────────────────────────────────────────────────────
lines = []
lines.append(SEP)
lines.append("  STUPA GEOGRAPHIC TIER AUDIT")
lines.append(f"  Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}")
lines.append(f"  Full corpus background: N = {N_CORPUS},  A+ = {K_AP_CORPUS} ({100*K_AP_CORPUS/N_CORPUS:.1f}%),  A = {K_A_CORPUS} ({100*K_A_CORPUS/N_CORPUS:.1f}%)")
lines.append(f"  Stupa keywords: {', '.join(STUPA_KEYWORDS)}")
lines.append(SEP)

# Summary table across all regions
lines.append("")
lines.append("  SUMMARY TABLE")
lines.append(f"  {'Region':<40} {'N':>4} {'A+':>4} {'A+%':>6} {'OR_A+':>6} {'p_A+':>8} {'A':>4} {'A%':>6} {'OR_A':>6} {'p_A':>8}")
lines.append(f"  {'─'*40} {'─'*4} {'─'*4} {'─'*6} {'─'*6} {'─'*8} {'─'*4} {'─'*6} {'─'*6} {'─'*8}")

for label, sites in region_data:
    n    = len(sites)
    tc   = tier_counts(sites)
    n_ap = tc["A+"]
    n_a  = tc["A+"] + tc["A"]
    or_ap, p_ap = fisher_ap(n, n_ap)
    or_a,  p_a  = fisher_a(n, n_a)
    rate_ap = f"{100*n_ap/n:.1f}%" if n else "—"
    rate_a  = f"{100*n_a/n:.1f}%"  if n else "—"
    p_ap_s  = f"{p_ap:.4f} {sig(p_ap)}" if p_ap == p_ap else "—"
    p_a_s   = f"{p_a:.4f} {sig(p_a)}"   if p_a  == p_a  else "—"
    lines.append(
        f"  {label:<40} {n:>4} {n_ap:>4} {rate_ap:>6} {str(or_ap)+'×':>6} {p_ap_s:>12} "
        f"{n_a:>4} {rate_a:>6} {str(or_a)+'×':>6} {p_a_s:>12}"
    )

# Per-region detail
for label, sites in region_data:
    lines.append("")
    lines.append(SEP2)
    lines.extend(tier_block(sites, label))
    lines.append("")
    lines.extend(site_listing(sites))

lines.append("")
lines.append(SEP)
lines.append("  CROSS-REGION COMPARISON: A+ AND A CONCENTRATION")
lines.append(SEP2)
lines.append("")
lines.append("  Key question: does the Java/Indonesia sub-region drive the stupa A+ signal,")
lines.append("  or is it spread across the full heartland?")
lines.append("")

# Java vs heartland-ex-Java comparison
java   = [s for s in stupa_sites if in_band(s["lon"], 105.0, 112.0)]
non_java_heartland = [s for s in stupa_sites
                      if in_band(s["lon"], 70.0, 110.0)
                      and not in_band(s["lon"], 105.0, 112.0)]
outside = [s for s in stupa_sites if not in_band(s["lon"], 70.0, 110.0)]

for grp_label, grp in [("Java (105°–112°E)", java),
                        ("Heartland ex-Java (70°–110°E, excl. 105°–112°E)", non_java_heartland),
                        ("Outside heartland", outside)]:
    n    = len(grp)
    tc   = tier_counts(grp)
    n_ap = tc["A+"]
    n_a  = tc["A+"] + tc["A"]
    or_ap, p_ap = fisher_ap(n, n_ap)
    or_a,  p_a  = fisher_a(n, n_a)
    rate_ap = f"{100*n_ap/n:.1f}%" if n else "—"
    rate_a  = f"{100*n_a/n:.1f}%"  if n else "—"
    lines.append(f"  {grp_label}")
    lines.append(f"    N={n}  A+={n_ap} ({rate_ap})  OR_A+={or_ap}× p={p_ap:.4f} {sig(p_ap)}")
    lines.append(f"          A={n_a}  ({rate_a})       OR_A={or_a}×  p={p_a:.4f}  {sig(p_a)}")
    lines.append(f"    Tier breakdown: " +
                 "  ".join(f"{t}={tier_counts(grp)[t]}" for t in TIER_ORDER))
    lines.append("")

lines.append(SEP)

output = "\n".join(lines)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(output, encoding="utf-8")
print(output)
print(f"\n  → Written to {OUT}")
