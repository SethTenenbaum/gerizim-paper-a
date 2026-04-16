"""
generate_audit_aplus_sites.py
=============================
Produce supplementary/audit/aplus_sites_audit.txt

Simple flat listing of every A-tier UNESCO Cultural/Mixed site
under two anchors: Gerizim and Jerusalem.

For each anchor, lists all A++ / A+ / A sites sorted by deviation
(tightest first). Includes a head-to-head summary and an overlap
section showing which sites appear in both.

Run from repo root:
    python3 tools/generate_audit_aplus_sites.py
"""

import sys
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APP, TIER_APLUS, TIER_A_MAX,
    P_NULL_AP, P_NULL_A,
    tier_label, is_aplus, is_a_or_better,
)
from scipy.stats import binomtest

OUT = Path(__file__).parent.parent / "supplementary" / "audit" / "aplus_sites_audit.txt"
SEP = "─" * 96

JERUSALEM = 35.2317   # from config anchor_comparison


def site_deviation(lon, anchor):
    arc  = abs(lon - anchor)
    bv   = arc / BERU
    near = round(bv / 0.1) * 0.1
    return abs(bv - near), near


def classify_sites(all_sites, anchor):
    """Return list of dicts for every A-tier site under the given anchor."""
    out = []
    for s in all_sites:
        dev, harmonic = site_deviation(s.longitude, anchor)
        t = tier_label(dev)
        if not is_a_or_better(t):
            continue
        direction = "E" if s.longitude >= anchor else "W"
        out.append({
            "name":      s.site,
            "lon":       s.longitude,
            "dev":       dev,
            "harmonic":  harmonic,
            "km":        dev * BERU * 111.0,
            "tier":      t,
            "direction": direction,
        })
    return sorted(out, key=lambda x: x["dev"])


def sig(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"


def run():
    corpus    = load_corpus()
    all_sites = cultural_sites_with_coords(corpus)
    N_ALL     = len(all_sites)

    # Gerizim is external — no exclusion needed; full corpus N=1011.
    # Jerusalem is an inscribed member and self-excludes when used as anchor
    # (matches Sweep B in the main paper, N=1010).
    jer_corpus = [s for s in all_sites if abs(s.longitude - JERUSALEM) >= 0.001]
    N_JER      = len(jer_corpus)   # should be N_ALL - 1

    ger  = classify_sites(all_sites, GERIZIM)
    jer  = classify_sites(jer_corpus, JERUSALEM)

    def counts(sites):
        app = sum(1 for s in sites if s["tier"] == "A++")
        ap  = sum(1 for s in sites if s["tier"] in ("A++", "A+"))
        a   = len(sites)
        return app, ap, a

    g_app, g_ap, g_a = counts(ger)
    j_app, j_ap, j_a = counts(jer)

    # Binomial p-values
    # Gerizim: full corpus N=1011; Jerusalem: self-excluded corpus N=1010
    p_ger_ap = binomtest(g_ap, N_ALL, P_NULL_AP, alternative="greater").pvalue
    p_jer_ap = binomtest(j_ap, N_JER, P_NULL_AP, alternative="greater").pvalue
    p_ger_a  = binomtest(g_a,  N_ALL, P_NULL_A,  alternative="greater").pvalue
    p_jer_a  = binomtest(j_a,  N_JER, P_NULL_A,  alternative="greater").pvalue

    # Overlap: sites in A+ for both anchors
    ger_ap_names = {s["name"] for s in ger if s["tier"] in ("A++", "A+")}
    jer_ap_names = {s["name"] for s in jer if s["tier"] in ("A++", "A+")}
    both_ap      = ger_ap_names & jer_ap_names
    ger_only_ap  = ger_ap_names - jer_ap_names
    jer_only_ap  = jer_ap_names - ger_ap_names

    ts = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")

    lines = [
        "UNESCO A-TIER SITES: GERIZIM vs JERUSALEM ANCHOR COMPARISON",
        f"Generated : {ts}",
        f"Script    : tools/generate_audit_aplus_sites.py",
        f"Corpus    : {N_ALL} Cultural/Mixed UNESCO sites with coordinates",
        f"Beru      : {BERU}°  |  A+ threshold: <= {TIER_APLUS} beru "
        f"({TIER_APLUS*BERU*111:.1f} km)  |  A threshold: <= {TIER_A_MAX} beru "
        f"({TIER_A_MAX*BERU*111:.0f} km)",
        f"Gerizim   : {GERIZIM}°E  (external; N={N_ALL}, no exclusion)",
        f"Jerusalem : {JERUSALEM}°E  (self-excluded as anchor; N={N_JER})  "
        f"[{abs(GERIZIM - JERUSALEM):.4f}° = {abs(GERIZIM - JERUSALEM)*111:.1f} km separation]",
        "",
    ]

    # ── Summary table ─────────────────────────────────────────────────────────
    lines += [
        SEP,
        "SUMMARY",
        "",
        f"  {'Anchor':<14}  {'N corpus':>9}  {'A++':>5}  {'A+':>5}  {'A (all)':>8}  "
        f"{'A+ %':>6}  {'A+ p':>9}  {'':>4}  {'A %':>6}  {'A p':>9}  {'':>4}",
        "  " + "-" * 86,
        f"  {'Gerizim':<14}  {N_ALL:>9}  {g_app:>5}  {g_ap:>5}  {g_a:>8}  "
        f"{100*g_ap/N_ALL:>5.1f}%  {p_ger_ap:>9.4f}  {sig(p_ger_ap):>4}  "
        f"{100*g_a/N_ALL:>5.1f}%  {p_ger_a:>9.4f}  {sig(p_ger_a):>4}",
        f"  {'Jerusalem':<14}  {N_JER:>9}  {j_app:>5}  {j_ap:>5}  {j_a:>8}  "
        f"{100*j_ap/N_JER:>5.1f}%  {p_jer_ap:>9.4f}  {sig(p_jer_ap):>4}  "
        f"{100*j_a/N_JER:>5.1f}%  {p_jer_a:>9.4f}  {sig(p_jer_a):>4}",
        "",
        f"  A+ null rate: {P_NULL_AP:.0%}   A null rate: {P_NULL_A:.0%}",
        f"  A+ overlap (both anchors): {len(both_ap)} sites",
        f"  A+ unique to Gerizim: {len(ger_only_ap)} sites",
        f"  A+ unique to Jerusalem: {len(jer_only_ap)} sites",
        "",
    ]

    # ── Gerizim site listing ──────────────────────────────────────────────────
    def site_rows(sites, label):
        rows = []
        rows.append(SEP)
        app_c = sum(1 for s in sites if s["tier"] == "A++")
        ap_c  = sum(1 for s in sites if s["tier"] in ("A++", "A+"))
        a_c   = len(sites)
        rows.append(f"{label.upper()} — ALL A-TIER SITES  "
                    f"(A++={app_c}, A+={ap_c}, A={a_c - ap_c}, total={a_c})")
        rows.append(
            f"  {'Tier':<5}  {'Dev (beru)':>10}  {'km':>6}  "
            f"{'Dir':>3}  {'Harmonic':>8}  {'Lon°E':>9}  Site"
        )
        rows.append("  " + "-" * 88)
        prev_tier = None
        for s in sites:
            if s["tier"] != prev_tier:
                if prev_tier is not None:
                    rows.append("")
                rows.append(f"  --- {s['tier']} ---")
                prev_tier = s["tier"]
            rows.append(
                f"  {s['tier']:<5}  {s['dev']:>10.6f}  {s['km']:>6.1f}  "
                f"{s['direction']:>3}  {s['harmonic']:>8.1f}  {s['lon']:>9.4f}  {s['name']}"
            )
        rows.append("")
        return rows

    lines += site_rows(ger, f"Gerizim ({GERIZIM}°E)")
    lines += site_rows(jer, f"Jerusalem ({JERUSALEM}°E)")

    # ── Overlap section ───────────────────────────────────────────────────────
    lines += [
        SEP,
        f"A+ OVERLAP: SITES IN A+ BAND FOR BOTH ANCHORS  (N={len(both_ap)})",
        "",
        f"  These sites fall within {TIER_APLUS*BERU*111:.1f} km of a 0.1-beru harmonic",
        f"  under both Gerizim and Jerusalem.  They do not discriminate between anchors.",
        "",
        f"  {'Site':<55}  {'Ger tier':>8}  {'Ger km':>7}  {'Jer tier':>8}  {'Jer km':>7}",
        "  " + "-" * 90,
    ]
    ger_by_name = {s["name"]: s for s in ger}
    jer_by_name = {s["name"]: s for s in jer}
    for name in sorted(both_ap):
        g = ger_by_name[name]
        j = jer_by_name[name]
        lines.append(
            f"  {name[:54]:<55}  {g['tier']:>8}  {g['km']:>7.1f}  {j['tier']:>8}  {j['km']:>7.1f}"
        )
    lines.append("")

    lines += [
        SEP,
        f"A+ UNIQUE TO GERIZIM  (N={len(ger_only_ap)}) — not A+ under Jerusalem",
        "",
        f"  {'Site':<55}  {'Ger tier':>8}  {'Ger km':>7}  {'Jer tier':>8}  {'Jer km':>7}",
        "  " + "-" * 90,
    ]
    ger_only_sorted = sorted(
        [s for s in ger if s["name"] in ger_only_ap],
        key=lambda x: x["dev"]
    )
    for s in ger_only_sorted:
        j_info = jer_by_name.get(s["name"])
        j_tier = j_info["tier"] if j_info else "—"
        j_km   = j_info["km"]   if j_info else 0.0
        j_km_s = f"{j_km:>7.1f}" if j_info else f"{'—':>7}"
        lines.append(
            f"  {s['name'][:54]:<55}  {s['tier']:>8}  {s['km']:>7.1f}  {j_tier:>8}  {j_km_s}"
        )
    lines.append("")

    lines += [
        SEP,
        f"A+ UNIQUE TO JERUSALEM  (N={len(jer_only_ap)}) — not A+ under Gerizim",
        "",
        f"  {'Site':<55}  {'Jer tier':>8}  {'Jer km':>7}  {'Ger tier':>8}  {'Ger km':>7}",
        "  " + "-" * 90,
    ]
    jer_only_sorted = sorted(
        [s for s in jer if s["name"] in jer_only_ap],
        key=lambda x: x["dev"]
    )
    for s in jer_only_sorted:
        g_info = ger_by_name.get(s["name"])
        g_tier = g_info["tier"] if g_info else "—"
        g_km   = g_info["km"]   if g_info else 0.0
        g_km_s = f"{g_km:>7.1f}" if g_info else f"{'—':>7}"
        lines.append(
            f"  {s['name'][:54]:<55}  {s['tier']:>8}  {s['km']:>7.1f}  {g_tier:>8}  {g_km_s}"
        )
    lines.append("")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(
        f"Written -> {OUT}\n"
        f"  Gerizim : A++={g_app}  A+={g_ap}  A={g_a}\n"
        f"  Jerusalem: A++={j_app}  A+={j_ap}  A={j_a}\n"
        f"  A+ overlap: {len(both_ap)}  |  Gerizim-unique: {len(ger_only_ap)}  "
        f"|  Jerusalem-unique: {len(jer_only_ap)}"
    )


if __name__ == "__main__":
    run()
