"""
generate_audit_interharmonic.py
================================
Produce supplementary/audit/interharmonic_audit.txt

Analysis of UNESCO Cultural/Mixed sites that fall near the midpoints
between 0.1-beru harmonics — the inter-harmonic zone of the harmonic grid.

A site at the midpoint between two adjacent harmonics has deviation = 0.05 beru
(the maximum possible distance, ~166 km). Sites with dev >= 0.04 beru are
within 6.6 km of a midpoint.

Each C-tier site is tagged with which keyword categories matched it:
dome, founding (any founding/sacred/landscape set), religion (any religion),
mound — so the distribution can be inspected across all keyword audits.

Run from repo root:
    python3 tools/generate_audit_interharmonic.py
"""

import sys
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter

sys.path.insert(0, str(Path(__file__).parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
import json as _json
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_A_MAX, TIER_B_MAX,
    P_NULL_AP, P_NULL_A, P_NULL_C, P_NULL_CMINUS,
    MIDPOINT, TIER_C_MAX, TIER_CMINUS, TIER_CMINUS2,
    deviation as beru_dev, tier_label, is_aplus, is_a_or_better,
    is_c_or_better, is_cminus_or_better,
    load_keywords, load_religion_sets,
)

def _load_mound_keywords():
    kw = _json.loads((Path(__file__).parent.parent / "keywords.json").read_text())
    m = kw.get("mound_evolution", {})
    return (
        m.get("mound_unambiguous", []) +
        m.get("mound_ambiguous", []) +
        m.get("stupa", []) +
        m.get("dome", [])
    )
from scipy.stats import binomtest, fisher_exact

OUT = Path(__file__).parent.parent / "supplementary" / "audit" / "interharmonic_audit.txt"
SEP = "─" * 96

# ── Keyword categories for tagging ───────────────────────────────────────────
# Each entry: (display_name, [keywords])
# Religion sets are aggregated under "religion" for the summary table,
# and broken out individually in the full listing.

def _build_categories():
    founding_kws = (
        load_keywords("founding_capital") +
        load_keywords("sacred_origin") +
        load_keywords("founding_monument") +
        load_keywords("founding_axis") +
        load_keywords("ancient_landscape")
    )
    return [
        ("dome",     load_keywords("dome_forms")),
        ("founding", founding_kws),
        ("mound",    _load_mound_keywords()),
        ("religion", [k for _, kws in load_religion_sets() for k in kws]),
    ]

def _build_detail_tags():
    """Detailed tag list for per-site listing (religion broken out by name)."""
    founding_kws = (
        load_keywords("founding_capital") +
        load_keywords("sacred_origin") +
        load_keywords("founding_monument") +
        load_keywords("founding_axis") +
        load_keywords("ancient_landscape")
    )
    tags = [
        ("dome",     load_keywords("dome_forms")),
        ("founding", founding_kws),
        ("mound",    _load_mound_keywords()),
    ]
    for rname, kws in load_religion_sets():
        tags.append((rname, kws))
    return tags


def sig(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"


def run():
    corpus    = load_corpus()
    all_sites = cultural_sites_with_coords(corpus)

    CATEGORIES   = _build_categories()
    DETAIL_TAGS  = _build_detail_tags()
    RELIGION_SETS = load_religion_sets()

    records = []
    for s in all_sites:
        if s.year is None:
            continue
        dev  = beru_dev(s.longitude)
        bv   = abs(s.longitude - GERIZIM) / BERU
        near = round(bv / 0.1) * 0.1
        ext  = s.extended_description if s.extended_description else ""
        text = (s.site + " " + s.short_description + " " + ext).lower()

        # High-level category tags
        cat_tags = [cname for cname, kws in CATEGORIES
                    if any(k in text for k in kws)]
        # Detailed tags (for full listing)
        detail_tags = [tname for tname, kws in DETAIL_TAGS
                       if any(k in text for k in kws)]

        dist_mid = 0.05 - dev
        records.append({
            "name":          s.site,
            "lon":           s.longitude,
            "dev":           dev,
            "dist_mid":      dist_mid,
            "dist_mid_km":   dist_mid * BERU * 111.0,
            "beru_val":      bv,
            "harmonic":      near,
            "km":            dev * BERU * 111.0,
            "tier":          tier_label(dev),
            "text":          text,
            "cat_tags":      cat_tags,
            "detail_tags":   detail_tags,
            "direction":     "E" if s.longitude >= GERIZIM else "W",
        })

    N              = len(records)
    interharmonic  = [r for r in records if is_c_or_better(r["tier"])]
    deep_inter     = [r for r in records if is_cminus_or_better(r["tier"])]
    n_inter        = len(interharmonic)
    n_deep         = len(deep_inter)

    ts = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")

    lines = [
        "UNESCO HARMONIC INTER-HARMONIC ZONE AUDIT",
        f"Generated  : {ts}",
        f"Script     : tools/generate_audit_interharmonic.py",
        f"Corpus     : {N} Cultural/Mixed UNESCO sites with year",
        f"Anchor     : Gerizim {GERIZIM}°E   |   Beru = {BERU}°",
        f"Harmonics  : 0.0, 0.1, 0.2 ... beru (step = 0.1 beru = {0.1*BERU*111:.1f} km)",
        f"Midpoints  : 0.05, 0.15, 0.25 ... beru (inter-harmonic midpoints)",
        f"Max dev    : 0.050 beru = {0.05*BERU*111:.1f} km  (perfect midpoint = C-- if dist_mid=0)",
        "",
        f"C-tier  (dist_mid <= {TIER_C_MAX} beru): within {TIER_C_MAX*BERU*111:.1f} km of a midpoint  [null rate = {P_NULL_C*100:.0f}%]",
        f"C- tier (dist_mid <= {TIER_CMINUS} beru): within {TIER_CMINUS*BERU*111:.1f} km of a midpoint  [null rate = {P_NULL_CMINUS*100:.0f}%]",
        f"C-- tier(dist_mid <= {TIER_CMINUS2} beru): within {TIER_CMINUS2*BERU*111:.2f} km of a midpoint",
        "",
        f"Keyword categories used for site tagging:",
        f"  dome     — dome_forms keyword set",
        f"  founding — founding_capital + sacred_origin + founding_monument + founding_axis + ancient_landscape",
        f"  mound    — mound_evolution (tumulus/barrow/stupa/dome forms)",
        f"  religion — any religion set (Buddhism, Christianity, Islam, Judaism, Hinduism, Zoroastrianism)",
        f"  (See keywords.json for full keyword lists; per-religion detail in religion_keyword_audit.txt)",
        "",
    ]

    # ── Tier distribution ─────────────────────────────────────────────────────
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
        f"  C-tier total (C + C- + C--):  {n_inter:>5}  {100*n_inter/N:>5.1f}%  (null rate = {P_NULL_C*100:.0f}%)",
        f"  C-/C-- total:                 {n_deep:>5}  {100*n_deep/N:>5.1f}%  (null rate = {P_NULL_CMINUS*100:.0f}%)",
        "",
    ]

    # ── Per-category summary ──────────────────────────────────────────────────
    lines += [
        SEP,
        "PER-CATEGORY INTER-HARMONIC SUMMARY",
        "",
        f"  How many C-tier sites fall into each keyword category?",
        f"  H0: each category has the same C-tier rate as the full corpus ({100*n_inter/N:.1f}%).",
        f"  Fisher exact (one-tailed, greater): category over-represented in C-tier.",
        "",
        f"  {'Category':<14}  {'N total':>8}  {'N C-tier':>9}  {'C%':>6}  {'OR':>6}  {'Fisher p':>9}  {''}",
        "  " + "-" * 72,
    ]
    for cname, kws in CATEGORIES:
        cat_sites = [r for r in records if any(k in r["text"] for k in kws)]
        n_cat = len(cat_sites)
        if n_cat < 5:
            continue
        n_cat_inter = sum(1 for r in cat_sites if is_c_or_better(r["tier"]))
        n_not = n_cat - n_cat_inter
        table = [[n_cat_inter, n_not],
                 [n_inter - n_cat_inter, (N - n_inter) - n_not]]
        OR, p = fisher_exact(table, alternative="greater")
        lines.append(
            f"  {cname:<14}  {n_cat:>8}  {n_cat_inter:>9}  {100*n_cat_inter/n_cat:>5.1f}%  "
            f"{OR:>6.2f}  {p:>9.4f}  {sig(p)}"
        )
    lines += [
        "  " + "-" * 72,
        f"  {'Full corpus':<14}  {N:>8}  {n_inter:>9}  {100*n_inter/N:>5.1f}%",
        "",
    ]

    # ── Spectrum symmetry ─────────────────────────────────────────────────────
    tightest_harmonic = sorted(records, key=lambda x: x["dev"])
    tightest_inter    = sorted(records, key=lambda x: x["dist_mid"])

    lines += [
        SEP,
        "SPECTRUM SYMMETRY: TIGHTEST HARMONIC vs TIGHTEST INTER-HARMONIC",
        "",
        f"  Top-5 sites closest to a harmonic and top-5 closest to a midpoint.",
        f"  Tags show which keyword categories matched each site.",
        "",
        f"  Tightest harmonic and inter-harmonic across all UNESCO sites:",
        "",
        f"  {'Site':<55}  {'Rank':>5}  {'Dist harmonic':>14}  {'Dist midpoint':>14}",
        "  " + "-" * 96,
    ]
    for rank, r in enumerate(tightest_harmonic[:5], 1):
        tags = ", ".join(r["detail_tags"]) if r["detail_tags"] else "—"
        lines.append(
            f"  {r['name'][:54]:<55}  #{rank:<4}  "
            f"{r['km']:>11.2f} km  {r['dist_mid_km']:>11.2f} km  [{tags}]"
        )
    lines.append("  ...")
    for rank, r in enumerate(tightest_inter[:5], 1):
        tags = ", ".join(r["detail_tags"]) if r["detail_tags"] else "—"
        lines.append(
            f"  {r['name'][:54]:<55}  #{rank:<4}  "
            f"{r['km']:>11.2f} km  {r['dist_mid_km']:>11.2f} km  [{tags}]  (inter-harmonic rank)"
        )
    lines.append("")

    # ── Full C-tier listing ───────────────────────────────────────────────────
    lines += [
        SEP,
        f"FULL C-TIER LISTING  (dist_mid <= {TIER_C_MAX} beru, N={n_inter})",
        f"Sorted by deviation descending (closest to midpoint first).",
        f"'Dist midpt' = km from the perfect inter-harmonic midpoint.",
        f"Tags show keyword categories: dome / founding / mound / religion name(s).",
        "",
        f"  {'Tier':<5}  {'Dev':>8}  {'km':>6}  {'Dist midpt':>10}  {'Near harm':>9}  {'Lon°E':>9}  {'Tags':<26}  Site",
        "  " + "-" * 110,
    ]
    for r in sorted(interharmonic, key=lambda x: -x["dev"]):
        tag_str = ", ".join(r["detail_tags"]) if r["detail_tags"] else "—"
        lines.append(
            f"  {r['tier']:<5}  {r['dev']:>8.5f}  {r['km']:>6.1f}  "
            f"{r['dist_mid_km']:>7.2f} km  "
            f"{r['harmonic']:>9.1f}  {r['lon']:>9.4f}  "
            f"{tag_str:<26}  {r['name'][:44]}"
        )
    lines.append("")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(
        f"Written -> {OUT}\n"
        f"  C-tier (dist_mid<=0.010): {n_inter}/{N} sites ({100*n_inter/N:.1f}%)\n"
        f"  C-/C-- (dist_mid<=0.002): {n_deep}/{N} sites ({100*n_deep/N:.1f}%)"
    )


if __name__ == "__main__":
    run()
