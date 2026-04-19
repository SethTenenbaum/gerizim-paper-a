"""
generate_audit_religion.py
==========================
Produce supplementary/audit/religion_keyword_audit.txt

Full per-religion breakdown with all sites, tiers, and distances.
Includes A-side (A+/A) and C-side (C/C-/C--) statistics and Fisher exact tests.

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
    TIER_C_MAX, TIER_CMINUS, TIER_CMINUS2, MIDPOINT,
    deviation as beru_deviation, tier_label,
    is_aplus, is_a_or_better, is_c_or_better, is_cminus_or_better,
    load_religion_sets,
    P_NULL_AP, P_NULL_A, P_NULL_C, P_NULL_CMINUS, P_NULL_CMINUS2,
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
        ext  = s.extended_description if s.extended_description else ""
        text = (s.site + " " + s.short_description + " " + ext).lower()
        sites.append({
            "name":   s.site,
            "lon":    s.longitude,
            "yr":     s.year,
            "dev":    d,
            "tier":   tier_label(d),
            "km":     d * BERU * 111.0,
            "text":   text,
            "states": s.states or "",
        })

    N_ALL    = len(sites)
    n_ap_all = sum(1 for s in sites if is_aplus(s["tier"]))
    n_a_all  = sum(1 for s in sites if is_a_or_better(s["tier"]))
    n_c_all  = sum(1 for s in sites if is_c_or_better(s["tier"]))
    n_cm_all = sum(1 for s in sites if is_cminus_or_better(s["tier"]))

    RELIGION_SETS = load_religion_sets()
    any_relig_kws = [k for _, kws in RELIGION_SETS for k in kws]
    any_relig     = [s for s in sites if any(k in s["text"] for k in any_relig_kws)]

    def stats(subset):
        n    = len(subset)
        nap  = sum(1 for s in subset if is_aplus(s["tier"]))
        na   = sum(1 for s in subset if is_a_or_better(s["tier"]))
        nc   = sum(1 for s in subset if is_c_or_better(s["tier"]))
        ncm  = sum(1 for s in subset if is_cminus_or_better(s["tier"]))
        enr_ap = (nap / n) / P_NULL_AP     if n else 0
        enr_a  = (na  / n) / P_NULL_A      if n else 0
        enr_c  = (nc  / n) / P_NULL_C      if n else 0
        enr_cm = (ncm / n) / P_NULL_CMINUS if n else 0
        # Fisher exact: sub-corpus vs rest of full corpus (one-sided, greater)
        nap_rest = n_ap_all - nap;  n_rest = N_ALL - n
        na_rest  = n_a_all  - na
        nc_rest  = n_c_all  - nc
        ncm_rest = n_cm_all - ncm
        def _f(k, n_sub, k_rest, n_rest_):
            if n_sub < 1: return 1.0, float('nan')
            tbl = [[k, n_sub - k], [k_rest, n_rest_ - k_rest]]
            or_, p_ = fisher_exact(tbl, alternative="greater")
            return p_, or_
        p_ap, or_ap = _f(nap, n, nap_rest, n_rest)
        p_a,  or_a  = _f(na,  n, na_rest,  n_rest)
        p_c,  or_c  = _f(nc,  n, nc_rest,  n_rest)
        p_cm, or_cm = _f(ncm, n, ncm_rest, n_rest)
        return n, nap, na, nc, ncm, p_ap, p_a, p_c, p_cm, enr_ap, enr_a, enr_c, enr_cm, or_ap, or_a, or_c, or_cm

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
        f"A-tier null: {P_NULL_A:.0%}  |  A+ null: {P_NULL_AP:.0%}  |  "
        f"C-tier null: {P_NULL_C:.0%}  |  C- null: {P_NULL_CMINUS:.0%}",
        f"Anchor: Gerizim {GERIZIM}°E   |   Beru = {BERU}°",
        f"Extended descriptions searched: YES",
        "",
    ]

    # ── A-side summary table ──────────────────────────────────────────────────
    lines += [
        SEP,
        "A-SIDE SUMMARY  (harmonic alignment)  — Fisher exact: sub-corpus vs rest of full corpus",
        f"  {'Religion':<20} {'N':>5}  {'A+':>4} {'A+%':>6} {'A+p(F)':>9} {'OR':>6} {'':>4}  "
        f"{'A':>4} {'A%':>6} {'Ap(F)':>9} {'OR':>6} {'':>4}",
        "  " + "-" * 96,
    ]

    religion_data = {}
    for religion, kws in RELIGION_SETS:
        matched = [s for s in sites if any(k in s["text"] for k in kws)]
        r = stats(matched)
        religion_data[religion] = (matched, r, kws)
        n, nap, na, nc, ncm, p_ap, p_a, p_c, p_cm, enr_ap, enr_a, enr_c, enr_cm, or_ap, or_a, *_ = r
        if n < 5:
            continue
        lines.append(
            f"  {religion:<20} {n:>5}  {nap:>4} {100*nap/n:>5.1f}% {p_ap:>9.4f} {or_ap:>5.2f}x {sig(p_ap):>4}  "
            f"{na:>4} {100*na/n:>5.1f}% {p_a:>9.4f} {or_a:>5.2f}x {sig(p_a):>4}"
        )

    r_u = stats(any_relig)
    n_u, nap_u, na_u, nc_u, ncm_u, p_ap_u, p_a_u, p_c_u, p_cm_u, enr_ap_u, enr_a_u, enr_c_u, enr_cm_u, or_ap_u, or_a_u, or_c_u, or_cm_u = r_u
    lines += [
        "  " + "-" * 96,
        f"  {'UNION (any)':<20} {n_u:>5}  {nap_u:>4} {100*nap_u/n_u:>5.1f}% {p_ap_u:>9.4f} {or_ap_u:>5.2f}x "
        f"{sig(p_ap_u):>4}  {na_u:>4} {100*na_u/n_u:>5.1f}% {p_a_u:>9.4f} {or_a_u:>5.2f}x {sig(p_a_u):>4}",
        "",
        f"  Full corpus: N={N_ALL}, A+ rate={100*n_ap_all/N_ALL:.1f}% ({n_ap_all}), A rate={100*n_a_all/N_ALL:.1f}% ({n_a_all})",
        f"  Fisher exact: one-sided (greater), sub-corpus vs rest of full corpus.",
        "",
    ]

    # ── C-side summary table ──────────────────────────────────────────────────
    lines += [
        SEP,
        "C-SIDE SUMMARY  (inter-harmonic midpoint proximity)  — Fisher exact vs rest of corpus",
        f"  Full-corpus C rate: {100*n_c_all/N_ALL:.1f}%  |  C-/C-- rate: {100*n_cm_all/N_ALL:.1f}%",
        "",
        f"  {'Religion':<20} {'N':>5}  {'C':>4} {'C%':>6} {'Cp(F)':>9} {'OR':>6} {'':>4}  "
        f"{'C-':>4} {'C-%':>6} {'C-p(F)':>9} {'OR':>6} {'':>4}",
        "  " + "-" * 96,
    ]

    for religion, kws in RELIGION_SETS:
        if religion not in religion_data:
            continue
        matched, r, _ = religion_data[religion]
        n, nap, na, nc, ncm, p_ap, p_a, p_c, p_cm, enr_ap, enr_a, enr_c, enr_cm, or_ap, or_a, or_c, or_cm = r
        if n < 5:
            continue
        lines.append(
            f"  {religion:<20} {n:>5}  {nc:>4} {100*nc/n:>5.1f}% {p_c:>9.4f} {or_c:>5.2f}x {sig(p_c):>4}  "
            f"{ncm:>4} {100*ncm/n:>5.1f}% {p_cm:>9.4f} {or_cm:>5.2f}x {sig(p_cm):>4}"
        )

    lines += [
        "  " + "-" * 96,
        f"  {'UNION (any)':<20} {n_u:>5}  {nc_u:>4} {100*nc_u/n_u:>5.1f}% {p_c_u:>9.4f} {or_c_u:>5.2f}x "
        f"{sig(p_c_u):>4}  {ncm_u:>4} {100*ncm_u/n_u:>5.1f}% {p_cm_u:>9.4f} {or_cm_u:>5.2f}x {sig(p_cm_u):>4}",
        "",
        f"  Full corpus: N={N_ALL}, C-tier count={n_c_all} ({100*n_c_all/N_ALL:.1f}%), C-/C-- count={n_cm_all} ({100*n_cm_all/N_ALL:.1f}%)",
        "",
    ]

    tier_order = {"A++": 0, "A+": 1, "A": 2, "B": 3, "C": 4, "C-": 5, "C--": 6}

    # ── Cross-religion overlap ────────────────────────────────────────────────
    # For each site tagged by 2+ religions, show which ones and its tier
    overlap_sites = []
    for s in sites:
        tags = [r for r, kws in RELIGION_SETS if any(k in s["text"] for k in kws)]
        if len(tags) >= 2:
            overlap_sites.append((s, tags))
    overlap_sites.sort(key=lambda x: (tier_order.get(x[0]["tier"], 9), x[0]["dev"]))

    # Pairwise overlap matrix
    relig_names = [r for r, _ in RELIGION_SETS]
    lines += [
        SEP,
        f"CROSS-RELIGION OVERLAP  ({len(overlap_sites)} sites tagged by 2+ religions)",
        "",
        "  Pairwise overlap counts (sites matching both religions):",
        "  " + " ".join(f"{r[:6]:>7}" for r in relig_names),
    ]
    for r1 in relig_names:
        kws1 = next(kws for r, kws in RELIGION_SETS if r == r1)
        row = f"  {r1:<14}"
        for r2 in relig_names:
            if r2 == r1:
                row += f"  {'—':>5}"
            else:
                kws2 = next(kws for r, kws in RELIGION_SETS if r == r2)
                n_both = sum(1 for s in sites
                             if any(k in s["text"] for k in kws1)
                             and any(k in s["text"] for k in kws2))
                row += f"  {n_both:>5}"
        lines.append(row)

    lines += [
        "",
        f"  {'Site':<55} {'Tier':<5} {'km':>7}  Tags",
        "  " + "-" * 88,
    ]
    for s, tags in overlap_sites:
        country = f"  ({s['states']})" if s.get("states") else ""
        lines.append(
            f"  {s['name'][:54]:<55} {s['tier']:<5} {s['km']:>7.1f}  {', '.join(tags)}{country}"
        )
    lines.append("")

    # ── Judaism Fisher exact (A+ and A-tier) — note results from main table ──
    jud_data = religion_data.get("Judaism")
    if jud_data:
        matched_jud, r_jud, _ = jud_data
        n_jud, nap_jud, na_jud, nc_jud, ncm_jud, p_ap_jud, p_a_jud, p_c_jud, p_cm_jud, \
            enr_ap_jud, enr_a_jud, enr_c_jud, enr_cm_jud, or_ap_jud, or_a_jud, *_ = r_jud
        non_jud_ap = n_ap_all - nap_jud
        non_jud_a  = n_a_all - na_jud
        non_jud_n  = N_ALL - n_jud
        non_jud_c  = n_c_all - nc_jud
        lines += [
            SEP,
            "JUDAISM: SMALL-SAMPLE DETAIL  (N=50; all p-values are Fisher exact, same as main table)",
            f"  N = {n_jud}  (full corpus N={N_ALL})",
            f"",
            f"  A+ tier:  {nap_jud}/{n_jud} ({100*nap_jud/n_jud:.1f}%)  vs rest: {non_jud_ap}/{non_jud_n} "
            f"({100*non_jud_ap/non_jud_n:.1f}%)  Fisher p={p_ap_jud:.4f} {sig(p_ap_jud)}  OR={or_ap_jud:.2f}x",
            f"  A tier:   {na_jud}/{n_jud} ({100*na_jud/n_jud:.1f}%)  vs rest: {non_jud_a}/{non_jud_n} "
            f"({100*non_jud_a/non_jud_n:.1f}%)  Fisher p={p_a_jud:.4f} {sig(p_a_jud)}  OR={or_a_jud:.2f}x",
            f"  C tier:   {nc_jud}/{n_jud} ({100*nc_jud/n_jud:.1f}%)  vs rest: {non_jud_c}/{non_jud_n} "
            f"({100*non_jud_c/non_jud_n:.1f}%)  Fisher p={p_c_jud:.4f} {sig(p_c_jud)}  OR=see main table",
            "",
            f"  NOTE: Judaism A-tier sites are concentrated in the outer A-band (7–33 km).",
            f"  Geographic clustering in the Levant and Central Europe is a strong",
            f"  confound at this sample size. Reported as observed association only.",
            "",
        ]

    # ── Per-religion site listings ─────────────────────────────────────────────

    for religion, kws in RELIGION_SETS:
        if religion not in religion_data:
            continue
        matched, r, _ = religion_data[religion]
        n, nap, na, nc, ncm, p_ap, p_a, p_c, p_cm, enr_ap, enr_a, enr_c, enr_cm, or_ap, or_a, or_c, or_cm = r
        if n < 1:
            continue

        matched_sorted = sorted(matched, key=lambda s: (tier_order.get(s["tier"], 9), s["dev"]))
        aplus_sites = [s for s in matched_sorted if is_aplus(s["tier"])]
        a_sites     = [s for s in matched_sorted if s["tier"] == "A"]
        b_sites     = [s for s in matched_sorted if s["tier"] == "B"]
        c_sites     = [s for s in matched_sorted if s["tier"] == "C"]
        cm_sites    = [s for s in matched_sorted if s["tier"] in ("C-", "C--")]

        # Per-keyword tier breakdown (sorted by total hits descending)
        kw_tier_rows = []
        for k in sorted(kws, key=lambda k: -sum(1 for s in matched if k in s["text"])):
            kw_sites = [s for s in matched if k in s["text"]]
            if not kw_sites:
                continue
            kw_ap  = sum(1 for s in kw_sites if is_aplus(s["tier"]))
            kw_a   = sum(1 for s in kw_sites if s["tier"] == "A")
            kw_b   = sum(1 for s in kw_sites if s["tier"] == "B")
            kw_c   = sum(1 for s in kw_sites if s["tier"] == "C")
            kw_cm  = sum(1 for s in kw_sites if s["tier"] in ("C-", "C--"))
            kw_tier_rows.append(
                f"    {k:<20}  total={len(kw_sites):>3}  A+={kw_ap:>2}  A={kw_a:>3}  B={kw_b:>3}  C={kw_c:>3}  C-={kw_cm:>2}"
            )

        lines += [
            SEP,
            f"{religion.upper()}  (N={n})",
            f"  Keywords : {', '.join(kws)}",
            f"  A+  : {nap}/{n} = {100*nap/n:.1f}%  enrich {enr_ap:.2f}x  Fisher p={p_ap:.4f} {sig(p_ap)}  OR={or_ap:.2f}x",
            f"  A   : {na}/{n} = {100*na/n:.1f}%  enrich {enr_a:.2f}x  Fisher p={p_a:.4f} {sig(p_a)}  OR={or_a:.2f}x",
            f"  C   : {nc}/{n} = {100*nc/n:.1f}%  enrich {enr_c:.2f}x  Fisher p={p_c:.4f} {sig(p_c)}  OR={or_c:.2f}x",
            f"  C-  : {ncm}/{n} = {100*ncm/n:.1f}%  enrich {enr_cm:.2f}x  Fisher p={p_cm:.4f} {sig(p_cm)}  OR={or_cm:.2f}x",
            "",
            f"  Keyword tier breakdown (sorted by total hits):",
            f"    {'keyword':<20}  {'total':>7}  {'A+':>4}  {'A':>5}  {'B':>5}  {'C':>5}  {'C-':>4}",
            "    " + "-" * 60,
        ] + kw_tier_rows + [""]

        def kw_hit(s):
            """Comma-separated keywords that matched this site."""
            return ", ".join(k for k in kws if k in s["text"])

        def site_line(s):
            country = f"  ({s['states']})" if s.get("states") else ""
            return (f"    {s['tier']:<4}  {s['km']:>6.1f} km  {s['lon']:>9.4f}°E"
                    f"  [{kw_hit(s)}]  {s['name']}{country}")

        def site_line_c(s):
            dist_mid = (MIDPOINT - s["dev"]) * BERU * 111.0
            country = f"  ({s['states']})" if s.get("states") else ""
            return (f"    {s['tier']:<4}  {s['km']:>6.1f} km  {s['lon']:>9.4f}°E"
                    f"  dist_mid={dist_mid:>5.1f} km  [{kw_hit(s)}]  {s['name']}{country}")

        if aplus_sites:
            lines.append(f"  A+ SITES ({len(aplus_sites)}):")
            for s in aplus_sites:
                lines.append(site_line(s))
            lines.append("")

        if a_sites:
            lines.append(f"  A-TIER SITES ({len(a_sites)}):")
            for s in a_sites:
                lines.append(site_line(s))
            lines.append("")

        if b_sites:
            lines.append(f"  B SITES ({len(b_sites)}):")
            for s in b_sites:
                lines.append(site_line(s))
            lines.append("")

        if c_sites:
            lines.append(f"  C SITES ({len(c_sites)}) — near inter-harmonic midpoint:")
            for s in c_sites:
                lines.append(site_line_c(s))
            lines.append("")

        if cm_sites:
            lines.append(f"  C-/C-- SITES ({len(cm_sites)}) — very close to midpoint:")
            for s in cm_sites:
                lines.append(site_line_c(s))
            lines.append("")

    # ── Union listing ─────────────────────────────────────────────────────────
    lines += [
        SEP,
        f"UNION — ALL RELIGION-TAGGED SITES (N={n_u})",
        f"  A+  : {nap_u}/{n_u} = {100*nap_u/n_u:.1f}%  enrich {enr_ap_u:.2f}x  p={p_ap_u:.4f} {sig(p_ap_u)}",
        f"  A   : {na_u}/{n_u} = {100*na_u/n_u:.1f}%  enrich {enr_a_u:.2f}x  p={p_a_u:.4f} {sig(p_a_u)}",
        f"  C   : {nc_u}/{n_u} = {100*nc_u/n_u:.1f}%  enrich {enr_c_u:.2f}x  p={p_c_u:.4f} {sig(p_c_u)}",
        f"  C-  : {ncm_u}/{n_u} = {100*ncm_u/n_u:.1f}%  enrich {enr_cm_u:.2f}x  p={p_cm_u:.4f} {sig(p_cm_u)}",
        "",
        f"  {'Site':<55} {'Tier':<5} {'km':>7}  Religion (keyword)",
        "  " + "-" * 90,
    ]
    any_relig_sorted = sorted(any_relig,
                              key=lambda s: (tier_order.get(s["tier"], 9), s["dev"]))
    for s in any_relig_sorted:
        tag_parts = []
        for r, rkws in RELIGION_SETS:
            hit_kws = [k for k in rkws if k in s["text"]]
            if hit_kws:
                tag_parts.append(f"{r}({', '.join(hit_kws)})")
        lines.append(
            f"  {s['name'][:54]:<55} {s['tier']:<5} {s['km']:>7.1f}  {'; '.join(tag_parts)}"
        )
    lines.append("")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines), encoding="utf-8")
    n_c_relig = sum(1 for s in any_relig if is_c_or_better(s["tier"]))
    print(f"Written → {OUT}  ({n_u} religion-tagged sites, {nap_u} A+, {n_c_relig} C-tier)")


if __name__ == "__main__":
    run()
