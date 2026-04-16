"""
generate_audit_antinode.py
==========================
Produce supplementary/audit/antinode_audit.txt

Analysis of UNESCO Cultural/Mixed sites that fall near the midpoints
between 0.1-beru harmonics — the "anti-nodes" of the harmonic grid.

A site at the midpoint between two adjacent harmonics has deviation = 0.05 beru
(the maximum possible distance, ~166 km). Sites with dev >= 0.04 beru are
within 6.6 km of a midpoint.

This audit documents the anti-node distribution, compares it across religion
keyword sub-corpora, and gives the full listing of deep anti-node sites.
Buddhist heritage is highlighted because it shows a bimodal pattern:
tight harmonic hits (Lumbini at the 1.6 beru harmonic) alongside a
disproportionate cluster of anti-node sites (Nara, Ajanta, Bodh Gaya, etc.).

Run from repo root:
    python3 tools/generate_audit_antinode.py
"""

import sys
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_A_MAX, TIER_B_MAX,
    P_NULL_AP, P_NULL_A, P_NULL_C, P_NULL_CMINUS,
    MIDPOINT, TIER_C_MAX, TIER_CMINUS, TIER_CMINUS2,
    deviation as beru_dev, tier_label, is_aplus, is_a_or_better,
    is_c_or_better, is_cminus_or_better,
    load_religion_sets,
)
from scipy.stats import binomtest, fisher_exact

OUT = Path(__file__).parent.parent / "supplementary" / "audit" / "antinode_audit.txt"
SEP = "─" * 96

# C-tier thresholds (mirror of A-tiers, measured from inter-harmonic midpoint)
# C   = dist_mid <= 0.010 beru  (~33 km)  ←→  A  = dev <= 0.010 beru
# C-  = dist_mid <= 0.002 beru  (~6.7 km) ←→  A+ = dev <= 0.002 beru
# C-- = dist_mid <= 0.0002 beru (~0.67 km)←→  A++ = dev <= 0.0002 beru
# These are equivalent to: dev >= 0.040, 0.048, 0.0498 respectively.
ANTINODE_DEV = MIDPOINT - TIER_C_MAX    # = 0.040  (C-tier boundary, for display)
DEEP_DEV     = MIDPOINT - TIER_CMINUS  # = 0.048  (C- tier boundary)


def sig(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"


def run():
    corpus    = load_corpus()
    all_sites = cultural_sites_with_coords(corpus)

    RELIGION_SETS = load_religion_sets()
    any_relig_kws = [k for _, kws in RELIGION_SETS for k in kws]

    records = []
    for s in all_sites:
        if s.year is None:
            continue
        dev  = beru_dev(s.longitude)
        bv   = abs(s.longitude - GERIZIM) / BERU
        near = round(bv / 0.1) * 0.1
        text = (s.short_description + " " + s.site).lower()
        tags = [rname for rname, kws in RELIGION_SETS if any(k in text for k in kws)]
        dist_mid = 0.05 - dev          # 0 = perfect anti-node, 0.05 = perfect harmonic
        records.append({
            "name":      s.site,
            "lon":       s.longitude,
            "dev":       dev,
            "dist_mid":  dist_mid,
            "dist_mid_km": dist_mid * BERU * 111.0,
            "beru_val":  bv,
            "harmonic":  near,
            "km":        dev * BERU * 111.0,
            "tier":      tier_label(dev),
            "text":      text,
            "tags":      tags,
            "direction": "E" if s.longitude >= GERIZIM else "W",
        })

    N              = len(records)
    antinode       = [r for r in records if is_c_or_better(r["tier"])]
    deep_antinode  = [r for r in records if is_cminus_or_better(r["tier"])]
    n_anti         = len(antinode)
    n_deep         = len(deep_antinode)

    ts = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")

    lines = [
        "UNESCO HARMONIC ANTI-NODE AUDIT",
        f"Generated  : {ts}",
        f"Script     : tools/generate_audit_antinode.py",
        f"Corpus     : {N} Cultural/Mixed UNESCO sites with year",
        f"Anchor     : Gerizim {GERIZIM}°E   |   Beru = {BERU}°",
        f"Harmonics  : 0.0, 0.1, 0.2 ... beru (step = 0.1 beru = {0.1*BERU*111:.1f} km)",
        f"Midpoints  : 0.05, 0.15, 0.25 ... beru (inter-harmonic nodes)",
        f"Max dev    : 0.050 beru = {0.05*BERU*111:.1f} km  (perfect inter-harmonic midpoint = C-- if dist_mid=0)",
        "",
        f"C-tier  (dist_mid <= {TIER_C_MAX} beru): within {TIER_C_MAX*BERU*111:.1f} km of a midpoint  [symmetric null rate = {P_NULL_C*100:.0f}%]",
        f"C- tier (dist_mid <= {TIER_CMINUS} beru): within {TIER_CMINUS*BERU*111:.1f} km of a midpoint  [symmetric null rate = {P_NULL_CMINUS*100:.0f}%]",
        f"C-- tier(dist_mid <= {TIER_CMINUS2} beru): within {TIER_CMINUS2*BERU*111:.2f} km of a midpoint",
        "",
    ]

    # ── Overview counts ───────────────────────────────────────────────────────
    from collections import Counter
    tier_counts = Counter(r["tier"] for r in records)

    lines += [
        SEP,
        "TIER DISTRIBUTION (all sites)",
        "",
        f"  {'Tier':<6}  {'N':>5}  {'%':>6}  {'Description'}",
        "  " + "-" * 70,
        f"  {'A++':<6}  {tier_counts['A++']:>5}  {100*tier_counts['A++']/N:>5.1f}%  dev <= 0.0002 beru  (~0.67 km from harmonic)",
        f"  {'A+':<6}  {tier_counts['A+']:>5}  {100*tier_counts['A+']/N:>5.1f}%  dev <= 0.002  beru  (~6.7 km from harmonic)",
        f"  {'A':<6}  {tier_counts['A']:>5}  {100*tier_counts['A']/N:>5.1f}%  dev <= 0.010  beru  (~33 km from harmonic)",
        f"  {'B':<6}  {tier_counts['B']:>5}  {100*tier_counts['B']/N:>5.1f}%  middle zone",
        f"  {'C':<6}  {tier_counts['C']:>5}  {100*tier_counts['C']/N:>5.1f}%  dist_mid <= 0.010 beru  (~33 km from midpoint)",
        f"  {'C-':<6}  {tier_counts['C-']:>5}  {100*tier_counts['C-']/N:>5.1f}%  dist_mid <= 0.002 beru  (~6.7 km from midpoint)",
        f"  {'C--':<6}  {tier_counts['C--']:>5}  {100*tier_counts['C--']/N:>5.1f}%  dist_mid <= 0.0002 beru (~0.67 km from midpoint)",
        "",
        f"  C-tier total (C + C- + C--):  {n_anti:>5}  {100*n_anti/N:>5.1f}%  (symmetric null rate = {P_NULL_C*100:.0f}%)",
        f"  C-/C-- total:                 {n_deep:>5}  {100*n_deep/N:>5.1f}%  (symmetric null rate = {P_NULL_CMINUS*100:.0f}%)",
        "",
    ]

    # ── Religion comparison ───────────────────────────────────────────────────
    lines += [
        SEP,
        f"RELIGION SUB-CORPUS VS C-TIER BAND  (dist_mid <= {TIER_C_MAX} beru = ~33 km from midpoint)",
        "",
        f"  H0: religion sub-corpus has same C-tier rate as full corpus ({100*n_anti/N:.1f}%)",
        f"  Fisher exact (one-tailed, greater): religion over-represented in C-tier band",
        "",
        f"  {'Religion':<18}  {'N':>5}  {'C-tier':>10}  {'%':>6}  {'OR':>6}  {'Fisher p':>9}  {''}",
        "  " + "-" * 72,
    ]

    for rname, kws in RELIGION_SETS:
        rel    = [r for r in records if rname in r["tags"]]
        n_rel  = len(rel)
        if n_rel < 8:
            continue
        n_mid  = sum(1 for r in rel if is_c_or_better(r["tier"]))
        n_not  = n_rel - n_mid
        table  = [[n_mid, n_not],
                  [n_anti - n_mid, (N - n_anti) - n_not]]
        OR, p  = fisher_exact(table, alternative="greater")
        lines.append(
            f"  {rname:<18}  {n_rel:>5}  {n_mid:>10}  {100*n_mid/n_rel:>5.1f}%  "
            f"{OR:>6.2f}  {p:>9.4f}  {sig(p)}"
        )

    # Full corpus row
    lines += [
        "  " + "-" * 72,
        f"  {'Full corpus':<18}  {N:>5}  {n_anti:>10}  {100*n_anti/N:>5.1f}%",
        "",
    ]

    # ── Buddhist detail ───────────────────────────────────────────────────────
    bkws     = next(kws for name, kws in RELIGION_SETS if "buddhis" in name.lower())
    buddhist = [r for r in records if any(k in r["text"] for k in bkws)]
    b_ap     = [r for r in buddhist if r["tier"] in ("A++", "A+")]
    b_a      = [r for r in buddhist if r["tier"] == "A"]
    b_anti   = [r for r in buddhist if is_c_or_better(r["tier"])]
    b_cminus = [r for r in buddhist if is_cminus_or_better(r["tier"])]

    lines += [
        SEP,
        f"BUDDHIST HERITAGE — BIMODAL PATTERN  (N={len(buddhist)})",
        "",
        f"  Buddhist sites show a bimodal deviation distribution:",
        f"  tight harmonic hits (A+/A++ tier) alongside a disproportionate",
        f"  cluster of C-tier sites (close to inter-harmonic midpoints).",
        "",
        f"  A+/A++ sites  : {len(b_ap):>3}  ({100*len(b_ap)/len(buddhist):.1f}%)  "
        f"— on harmonic lines",
        f"  A sites       : {len(b_a):>3}  ({100*len(b_a)/len(buddhist):.1f}%)  "
        f"— inner harmonic band 7-33 km",
        f"  C-tier sites  : {len(b_anti):>3}  ({100*len(b_anti)/len(buddhist):.1f}%)  "
        f"— near midpoints (expected {100*n_anti/N:.1f}%)",
        f"  C-/C-- sites  : {len(b_cminus):>3}  ({100*len(b_cminus)/len(buddhist):.1f}%)  "
        f"— very close to midpoints (expected {P_NULL_CMINUS*100:.0f}%)",
        "",
    ]

    # Quartile breakdown
    devs_sorted = sorted(r["dev"] for r in records)
    q1 = devs_sorted[N // 4]
    q2 = devs_sorted[N // 2]
    q3 = devs_sorted[3 * N // 4]
    lines += [
        f"  Deviation quartiles (full corpus):  "
        f"Q1={q1:.4f}  Q2={q2:.4f}  Q3={q3:.4f}  max=0.0500",
        "",
        f"  {'Quartile':<26}  {'Buddhist':>9}  {'All sites':>9}  {'Enrichment'}",
        "  " + "-" * 60,
    ]
    for label, lo, hi in [
        ("Q1 tight  (0.000 – {:.4f})".format(q1), 0.0, q1),
        ("Q2        ({:.4f} – {:.4f})".format(q1, q2), q1, q2),
        ("Q3        ({:.4f} – {:.4f})".format(q2, q3), q2, q3),
        ("Q4 anti-node ({:.4f} – 0.050)".format(q3), q3, 1.0),
    ]:
        nb = sum(1 for r in buddhist if lo <= r["dev"] < hi)
        na = sum(1 for r in records  if lo <= r["dev"] < hi)
        enr = (nb / len(buddhist)) / (na / N) if na else 0
        lines.append(
            f"  {label:<26}  {nb:>4}/{len(buddhist):>3} = {100*nb/len(buddhist):>4.0f}%  "
            f"{na:>4}/{N:>4} = {100*na/N:>4.0f}%  {enr:.2f}x"
        )
    lines.append("")

    # Buddhist A+/A sites (tight hits)
    lines += [
        f"  BUDDHIST TIGHT HITS (A+ and A-tier) — on harmonic lines:",
        "",
        f"  {'Tier':<5}  {'Dev (beru)':>10}  {'km':>6}  {'Harmonic':>8}  {'Lon°E':>9}  Site",
        "  " + "-" * 80,
    ]
    for r in sorted(b_ap + b_a, key=lambda x: x["dev"]):
        lines.append(
            f"  {r['tier']:<5}  {r['dev']:>10.5f}  {r['km']:>6.1f}  "
            f"{r['harmonic']:>8.1f}  {r['lon']:>9.4f}  {r['name']}"
        )
    lines.append("")

    # Buddhist C-tier sites — with distance-from-midpoint column
    lines += [
        f"  BUDDHIST C-TIER SITES (dist_mid <= {TIER_C_MAX} beru) — near inter-harmonic midpoints:",
        f"  'Dist midpt' = distance from the perfect inter-harmonic midpoint (0 = exactly between harmonics)",
        "",
        f"  {'Tier':<5}  {'Dev (beru)':>10}  {'km':>6}  {'Dist midpt':>12}  {'Near harm':>9}  {'Lon°E':>9}  Site",
        "  " + "-" * 94,
    ]
    for r in sorted(b_anti, key=lambda x: -x["dev"]):
        lines.append(
            f"  {r['tier']:<5}  {r['dev']:>10.5f}  {r['km']:>6.1f}  "
            f"{r['dist_mid_km']:>9.2f} km  "
            f"{r['harmonic']:>9.1f}  {r['lon']:>9.4f}  {r['name']}"
        )
    lines.append("")

    # ── Spectrum symmetry: tightest harmonic vs tightest anti-node ────────────
    # Rank all sites by closeness to harmonic (dev) and closeness to midpoint (dist_mid)
    tightest_harmonic  = sorted(records, key=lambda x: x["dev"])
    tightest_antinode  = sorted(records, key=lambda x: x["dist_mid"])

    # Rank of Lumbini in each list
    lumbini_h_rank = next(i+1 for i, r in enumerate(tightest_harmonic) if "lumbini" in r["name"].lower())
    lumbini_a_rank = next(i+1 for i, r in enumerate(tightest_antinode) if "lumbini" in r["name"].lower())
    nara_h_rank    = next((i+1 for i, r in enumerate(tightest_harmonic) if "nara" in r["name"].lower()), None)
    nara_a_rank    = next((i+1 for i, r in enumerate(tightest_antinode) if "nara" in r["name"].lower()), None)

    lum = next(r for r in records if "lumbini" in r["name"].lower())
    nara = next((r for r in records if "nara" in r["name"].lower() and "monuments" in r["name"].lower()), None)

    lines += [
        SEP,
        "SPECTRUM SYMMETRY: TIGHTEST HARMONIC vs TIGHTEST ANTI-NODE",
        "",
        f"  The harmonic system has two reference extremes:",
        f"    Harmonic line  (dev = 0.000) — sites precisely on a 0.1-beru multiple",
        f"    Midpoint       (dev = 0.050) — sites exactly between two harmonics",
        "",
        f"  Buddhism occupies both extremes at comparable precision (~1 km):",
        "",
        f"  {'Site':<55}  {'Rank':>5}  {'Dist harmonic':>14}  {'Dist midpoint':>14}",
        "  " + "-" * 96,
        f"  {lum['name']:<55}  #{lumbini_h_rank:<4}  "
        f"{lum['dev']*BERU*111:>11.2f} km  {lum['dist_mid_km']:>11.2f} km  (tightest Buddhist harmonic hit)",
    ]
    if nara:
        lines.append(
            f"  {nara['name']:<55}  #{nara_a_rank:<4}  "
            f"{nara['dev']*BERU*111:>11.2f} km  {nara['dist_mid_km']:>11.2f} km  (tightest Buddhist anti-node)"
        )

    lines += [
        "",
        f"  For reference — tightest harmonic and anti-node across all UNESCO sites:",
        "",
        f"  {'Site':<55}  {'Rank':>5}  {'Dist harmonic':>14}  {'Dist midpoint':>14}",
        "  " + "-" * 96,
    ]
    for rank, r in enumerate(tightest_harmonic[:5], 1):
        tags = ", ".join(r["tags"]) if r["tags"] else "—"
        lines.append(
            f"  {r['name'][:54]:<55}  #{rank:<4}  "
            f"{r['km']:>11.2f} km  {r['dist_mid_km']:>11.2f} km  [{tags}]"
        )
    lines.append("  ...")
    for rank, r in enumerate(tightest_antinode[:5], 1):
        tags = ", ".join(r["tags"]) if r["tags"] else "—"
        lines.append(
            f"  {r['name'][:54]:<55}  #{rank:<4}  "
            f"{r['km']:>11.2f} km  {r['dist_mid_km']:>11.2f} km  [{tags}]  (anti-node rank)"
        )
    lines.append("")

    # ── Full C-tier listing (all religions) ───────────────────────────────────
    lines += [
        SEP,
        f"FULL C-TIER LISTING  (dist_mid <= {TIER_C_MAX} beru, N={n_anti})",
        f"Sorted by deviation descending (closest to midpoint first).",
        f"'Dist midpt' = km from the perfect inter-harmonic midpoint.",
        "",
        f"  {'Tier':<5}  {'Dev':>8}  {'km':>6}  {'Dist midpt':>10}  {'Near harm':>9}  {'Lon°E':>9}  {'Tags':<22}  Site",
        "  " + "-" * 106,
    ]
    for r in sorted(antinode, key=lambda x: -x["dev"]):
        tag_str = ", ".join(r["tags"]) if r["tags"] else "—"
        lines.append(
            f"  {r['tier']:<5}  {r['dev']:>8.5f}  {r['km']:>6.1f}  "
            f"{r['dist_mid_km']:>7.2f} km  "
            f"{r['harmonic']:>9.1f}  {r['lon']:>9.4f}  "
            f"{tag_str:<22}  {r['name'][:46]}"
        )
    lines.append("")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(
        f"Written -> {OUT}\n"
        f"  C-tier (dist_mid<=0.010): {n_anti}/{N} sites ({100*n_anti/N:.1f}%)\n"
        f"  C-/C-- (dist_mid<=0.002): {n_deep}/{N} sites ({100*n_deep/N:.1f}%)\n"
        f"  Buddhist C-tier: {len(b_anti)}/{len(buddhist)} = "
        f"{100*len(b_anti)/len(buddhist):.1f}% (global rate {100*n_anti/N:.1f}%)\n"
        f"  Buddhist C-/C--: {len(b_cminus)}/{len(buddhist)} = "
        f"{100*len(b_cminus)/len(buddhist):.1f}% (global rate {P_NULL_CMINUS*100:.0f}%)"
    )


if __name__ == "__main__":
    run()
