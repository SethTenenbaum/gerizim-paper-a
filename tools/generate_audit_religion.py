"""
generate_audit_religion.py
==========================
Produce supplementary/audit/religion_keyword_audit.txt

Full per-religion breakdown with all sites, tiers, and distances.
Includes the Judaism A-tier Fisher exact result discussed in the paper.

Run from repo root:
    python3 tools/generate_audit_religion.py
"""

import sys
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_A_MAX,
    deviation as beru_deviation, tier_label, is_aplus, is_a_or_better,
    load_religion_sets, P_NULL_AP, P_NULL_A,
)
from scipy.stats import binomtest, fisher_exact

OUT = Path(__file__).parent.parent / "supplementary" / "audit" / "religion_keyword_audit.txt"
SEP = "─" * 88


def run():
    corpus   = load_corpus()
    cultural = cultural_sites_with_coords(corpus)

    sites = []
    for s in cultural:
        if s.year is None:
            continue
        d    = beru_deviation(s.longitude)
        text = (s.short_description + " " + s.site).lower()
        sites.append({
            "name": s.site,
            "lon":  s.longitude,
            "yr":   s.year,
            "dev":  d,
            "tier": tier_label(d),
            "km":   d * BERU * 111.0,
            "text": text,
        })

    N_ALL    = len(sites)
    n_a_all  = sum(1 for s in sites if is_a_or_better(s["tier"]))

    RELIGION_SETS = load_religion_sets()
    any_relig_kws = [k for _, kws in RELIGION_SETS for k in kws]
    any_relig     = [s for s in sites if any(k in s["text"] for k in any_relig_kws)]

    def stats(subset):
        n   = len(subset)
        nap = sum(1 for s in subset if is_aplus(s["tier"]))
        na  = sum(1 for s in subset if is_a_or_better(s["tier"]))
        p_ap = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue if n >= 1 else 1.0
        p_a  = binomtest(na,  n, P_NULL_A,  alternative="greater").pvalue if n >= 1 else 1.0
        enr_ap = (nap / n) / P_NULL_AP if n else 0
        enr_a  = (na  / n) / P_NULL_A  if n else 0
        return n, nap, na, p_ap, p_a, enr_ap, enr_a

    def sig(p):
        if p < 0.001: return "***"
        if p < 0.01:  return "**"
        if p < 0.05:  return "*"
        if p < 0.10:  return "~"
        return "ns"

    ts = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")
    lines = [
        "UNESCO WORLD RELIGION KEYWORD AUDIT",
        f"Generated : {ts}",
        f"Script    : tools/generate_audit_religion.py",
        f"Full corpus (Cultural/Mixed with year): N = {N_ALL}",
        f"A-tier null: {P_NULL_A:.0%}   |   A+ null: {P_NULL_AP:.0%}",
        f"Anchor: Gerizim {GERIZIM}°E   |   Beru = {BERU}°",
        "",
    ]

    # ── Summary table ─────────────────────────────────────────────────────────
    lines += [
        SEP,
        "SUMMARY TABLE",
        f"  {'Religion':<20} {'N':>5}  {'A+':>4} {'A+%':>6} {'A+p':>8} {'':>4}  "
        f"{'A':>4} {'A%':>6} {'Ap':>8} {'':>4}",
        "  " + "-" * 82,
    ]

    religion_data = {}
    for religion, kws in RELIGION_SETS:
        matched = [s for s in sites if any(k in s["text"] for k in kws)]
        n, nap, na, p_ap, p_a, enr_ap, enr_a = stats(matched)
        religion_data[religion] = (matched, n, nap, na, p_ap, p_a, enr_ap, enr_a, kws)
        if n < 5:
            continue
        lines.append(
            f"  {religion:<20} {n:>5}  {nap:>4} {100*nap/n:>5.1f}% {p_ap:>8.4f} {sig(p_ap):>4}  "
            f"{na:>4} {100*na/n:>5.1f}% {p_a:>8.4f} {sig(p_a):>4}"
        )

    n_u, nap_u, na_u, p_ap_u, p_a_u, enr_ap_u, enr_a_u = stats(any_relig)
    lines += [
        "  " + "-" * 82,
        f"  {'UNION (any)':<20} {n_u:>5}  {nap_u:>4} {100*nap_u/n_u:>5.1f}% {p_ap_u:>8.4f} "
        f"{sig(p_ap_u):>4}  {na_u:>4} {100*na_u/n_u:>5.1f}% {p_a_u:>8.4f} {sig(p_a_u):>4}",
        "",
        f"  Full corpus: N={N_ALL}, A+ null=4%, A null=20%",
        "",
    ]

    # ── Judaism Fisher exact (A-tier) ─────────────────────────────────────────
    jud_data = religion_data.get("Judaism")
    if jud_data:
        matched_jud, n_jud, nap_jud, na_jud, *_ = jud_data
        non_jud_a = n_a_all - na_jud
        non_jud_n = N_ALL - n_jud
        table = [[na_jud, n_jud - na_jud], [non_jud_a, non_jud_n - non_jud_a]]
        fisher_or, fisher_p = fisher_exact(table, alternative="greater")
        lines += [
            SEP,
            "JUDAISM: SMALL-SAMPLE DETAIL (Fisher exact for A-tier)",
            f"  N = {n_jud}  |  A-tier = {na_jud}/{ n_jud} ({100*na_jud/n_jud:.1f}%)",
            f"  vs rest of corpus: {non_jud_a}/{non_jud_n} ({100*non_jud_a/non_jud_n:.1f}%)",
            f"  Fisher exact (one-sided): p = {fisher_p:.4f}  {sig(fisher_p)}   OR = {fisher_or:.2f}x",
            f"",
            f"  NOTE: 8 of 9 A-tier Judaism sites are in the outer band (7-33 km, not A+).",
            f"  Geographic concentration in the Levant and Central Europe is a strong",
            f"  confound at this sample size. Reported as observed association only.",
            "",
        ]

    # ── Per-religion site listings ─────────────────────────────────────────────
    tier_order = {"A++": 0, "A+": 1, "A": 2, "B": 3, "C": 4, "C-": 5, "C--": 6}

    for religion, kws in RELIGION_SETS:
        if religion not in religion_data:
            continue
        matched, n, nap, na, p_ap, p_a, enr_ap, enr_a, _ = religion_data[religion]
        if n < 1:
            continue

        matched_sorted = sorted(matched, key=lambda s: (tier_order.get(s["tier"], 9), s["dev"]))
        aplus_sites = [s for s in matched_sorted if is_aplus(s["tier"])]
        a_sites     = [s for s in matched_sorted if s["tier"] == "A"]
        other_sites = [s for s in matched_sorted if s["tier"] not in ("A++", "A+", "A")]

        lines += [
            SEP,
            f"{religion.upper()}  (N={n}, keywords: {kws})",
            f"  A+  : {nap}/{n} = {100*nap/n:.1f}%  enrich {enr_ap:.2f}x  p={p_ap:.4f} {sig(p_ap)}",
            f"  A   : {na}/{n} = {100*na/n:.1f}%  enrich {enr_a:.2f}x  p={p_a:.4f} {sig(p_a)}",
            "",
        ]

        if aplus_sites:
            lines.append(f"  A+ SITES ({len(aplus_sites)}):")
            for s in aplus_sites:
                lines.append(f"    {s['tier']:<4}  {s['km']:>6.1f} km  {s['lon']:>9.4f}°E  {s['name']}")
            lines.append("")

        if a_sites:
            lines.append(f"  A-TIER SITES ({len(a_sites)}):")
            for s in a_sites:
                lines.append(f"    {s['tier']:<4}  {s['km']:>6.1f} km  {s['lon']:>9.4f}°E  {s['name']}")
            lines.append("")

        if other_sites:
            lines.append(f"  B/C SITES ({len(other_sites)}):")
            for s in other_sites:
                lines.append(f"    {s['tier']:<4}  {s['km']:>6.1f} km  {s['lon']:>9.4f}°E  {s['name']}")
            lines.append("")

    # ── Union listing ─────────────────────────────────────────────────────────
    lines += [
        SEP,
        f"UNION — ALL RELIGION-TAGGED SITES (N={n_u})",
        f"  A+  : {nap_u}/{n_u} = {100*nap_u/n_u:.1f}%  enrich {enr_ap_u:.2f}x  p={p_ap_u:.4f} {sig(p_ap_u)}",
        f"  A   : {na_u}/{n_u} = {100*na_u/n_u:.1f}%  enrich {enr_a_u:.2f}x  p={p_a_u:.4f} {sig(p_a_u)}",
        "",
        f"  {'Site':<55} {'Tier':<5} {'km':>7}  Religion tags",
        "  " + "-" * 82,
    ]
    any_relig_sorted = sorted(any_relig,
                              key=lambda s: (tier_order.get(s["tier"], 9), s["dev"]))
    for s in any_relig_sorted:
        tags = [r for r, kws in RELIGION_SETS
                if any(k in s["text"] for k in kws)]
        lines.append(
            f"  {s['name'][:54]:<55} {s['tier']:<5} {s['km']:>7.1f}  {', '.join(tags)}"
        )
    lines.append("")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Written → {OUT}  ({n_u} religion-tagged sites, {nap_u} A+)")


if __name__ == "__main__":
    run()
