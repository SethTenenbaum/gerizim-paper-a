"""
generate_audit_corridor.py
==========================
Produce supplementary/audit/corridor_audit.txt

Full corridor precision audit: Gerizim vs Jerusalem anchor comparison.
Lists all 17 harmonics, their A+ sites, occupancy stats, and the
anchor-comparison table showing why Gerizim outperforms Jerusalem.

Run from repo root:
    python3 tools/generate_audit_corridor.py
"""

import json
import math
import sys
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import (
    GERIZIM, BERU, TIER_APP, TIER_APLUS, TIER_A_MAX, TIER_B_MAX,
    P_NULL_AP,
    deviation as beru_deviation, tier_label, is_aplus, is_a_or_better,
    load_notable_anchors,
)
from scipy.stats import binomtest, fisher_exact, mannwhitneyu

OUT = Path(__file__).parent.parent / "supplementary" / "audit" / "corridor_audit.txt"
SEP = "─" * 100

_ROOT = Path(__file__).parent.parent
with open(_ROOT / "config.json") as _f:
    _CONFIG = json.load(_f)

_corr = _CONFIG["corridor"]
LUMBINI_LON    = _corr["lumbini_lon"]
CORRIDOR_MAX   = _corr["corridor_max_beru"]
SURVEY_MAX     = _corr["survey_max_beru"]

CORRIDOR_HARMONICS = [round(n * 0.1, 1) for n in range(0, int(CORRIDOR_MAX * 10) + 1)]
SURVEY_HARMONICS   = [round(n * 0.1, 1) for n in range(0, int(SURVEY_MAX  * 10) + 1)]
N_HARMONICS = len(CORRIDOR_HARMONICS)


def best_site_at_harmonic(harmonic, sites, anchor=GERIZIM):
    best = None
    for s in sites:
        arc  = abs(s.longitude - anchor)
        bv   = arc / BERU
        near = round(bv * 10) / 10
        if near != harmonic:
            continue
        dev = abs(bv - near)
        if dev <= TIER_APLUS:
            if best is None or dev < best["dev"]:
                direction = "E" if s.longitude >= anchor else "W"
                best = {
                    "site": s.site, "lon": s.longitude, "dev": dev,
                    "km": dev * BERU * 111, "tier": tier_label(dev),
                    "direction": direction,
                }
    return best


def all_sites_at_harmonic(harmonic, sites, anchor=GERIZIM, max_tier="A+"):
    out = []
    for s in sites:
        arc  = abs(s.longitude - anchor)
        bv   = arc / BERU
        near = round(bv * 10) / 10
        if near != harmonic:
            continue
        dev  = abs(bv - near)
        tier = tier_label(dev)
        if is_aplus(tier):
            direction = "E" if s.longitude >= anchor else "W"
            out.append({
                "site": s.site, "lon": s.longitude, "dev": dev,
                "km": dev * BERU * 111, "tier": tier, "direction": direction,
            })
    return sorted(out, key=lambda x: x["dev"])


def corridor_occupancy(anchor_lon, sites, n_harmonics=17):
    harmonics = [round(n * 0.1, 1) for n in range(0, n_harmonics)]
    hits = 0
    app  = 0
    for n in harmonics:
        best_dev = None
        for s in sites:
            arc  = abs(s.longitude - anchor_lon)
            bv   = arc / BERU
            near = round(bv * 10) / 10
            if near != n:
                continue
            dev = abs(bv - near)
            if dev <= TIER_APLUS:
                if best_dev is None or dev < best_dev:
                    best_dev = dev
        if best_dev is not None:
            hits += 1
            if best_dev <= TIER_APP:
                app += 1
    return hits, app


def sig(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"


def run():
    corpus    = load_corpus()
    all_sites = cultural_sites_with_coords(corpus)

    ts = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")

    lines = [
        "GERIZIM CORRIDOR PRECISION AUDIT",
        f"Generated : {ts}",
        f"Script    : tools/generate_audit_corridor.py",
        f"Anchor    : Gerizim {GERIZIM}°E   |   Beru = {BERU}°",
        f"Corridor  : {GERIZIM}°E (0.0 beru) to Lumbini {LUMBINI_LON}°E "
        f"({abs(LUMBINI_LON - GERIZIM)/BERU:.4f} beru ≈ 1.6 harmonic)",
        f"Harmonics : {N_HARMONICS} (0.0 to {CORRIDOR_MAX:.1f} beru, step 0.1)",
        f"A+ thresh : ≤ {TIER_APLUS} beru (≤ {TIER_APLUS*BERU*111:.1f} km)",
        f"Corpus    : {len(all_sites)} UNESCO Cultural/Mixed sites with coordinates",
        "",
    ]

    # ── Per-harmonic table ────────────────────────────────────────────────────
    lines += [
        SEP,
        "PER-HARMONIC SITE TABLE",
        f"  {'Harm':>5}  {'Harmonic°E':>10}  {'Dir':>3}  {'Best A+ site':<52}  {'Dev':>8}  {'km':>5}  Tier",
        "  " + "-" * 90,
    ]

    corridor_hits   = []
    corridor_misses = []

    for n in CORRIDOR_HARMONICS:
        harmonic_e_lon = GERIZIM + n * BERU
        harmonic_w_lon = GERIZIM - n * BERU
        best      = best_site_at_harmonic(n, all_sites)
        all_ap    = all_sites_at_harmonic(n, all_sites)

        if best:
            corridor_hits.append({**best, "harmonic": n,
                                   "harmonic_lon": harmonic_e_lon, "all_ap": all_ap})
            marker    = " A++" if best["tier"] == "A++" else " A+"
            dir_str   = best["direction"]
            lines.append(
                f"  {n:>5.1f}  {harmonic_e_lon:>10.3f}  {dir_str:>3}  "
                f"{best['site']:<52}  {best['dev']:>8.6f}  {best['km']:>5.1f}  {best['tier']}{marker}"
            )
        else:
            corridor_misses.append(n)
            lines.append(
                f"  {n:>5.1f}  {harmonic_e_lon:>10.3f}    E  "
                f"{'--- no A+ site ---':<52}"
            )

    n_hits = len(corridor_hits)
    n_app  = sum(1 for h in corridor_hits if h["tier"] == "A++")

    lines += [
        "",
        f"  Corridor occupancy: {n_hits}/{N_HARMONICS} harmonics have >= 1 A+ site",
        f"  A++: {n_app}   A+: {n_hits - n_app}   Misses: {len(corridor_misses)}",
        "",
    ]

    # ── Full A+ catalogue per harmonic ────────────────────────────────────────
    lines += [
        SEP,
        "FULL A+ SITE CATALOGUE PER CORRIDOR HARMONIC",
        "",
    ]
    for h in corridor_hits:
        harmonic_lon_e = GERIZIM + h["harmonic"] * BERU
        harmonic_lon_w = GERIZIM - h["harmonic"] * BERU
        lines.append(
            f"  {h['harmonic']:.1f} beru  "
            f"(E: {harmonic_lon_e:.3f}°E  |  W: {harmonic_lon_w:.3f}°E)"
        )
        for s in h["all_ap"]:
            marker = " [A++]" if s["tier"] == "A++" else " [A+]"
            lines.append(
                f"      {s['direction']}  {s['site']:<55}  "
                f"lon={s['lon']:>9.4f}  dev={s['dev']:.6f}  {s['km']:>5.1f} km{marker}"
            )
        lines.append("")

    # ── Statistical tests ─────────────────────────────────────────────────────
    bt      = binomtest(n_hits, N_HARMONICS, P_NULL_AP, alternative="greater")
    p_binom = bt.pvalue

    # Fisher combined probability method
    log_p_sum = 0.0
    chi_obs   = 0.0
    for h in corridor_hits:
        p_i = h["dev"] / TIER_B_MAX
        p_i = max(p_i, 1e-15)
        log_p_sum += math.log(p_i)
        chi_obs   += -2 * math.log(p_i)
    from scipy.stats import chi2
    df_fisher = 2 * n_hits
    p_fisher  = chi2.sf(chi_obs, df_fisher)

    # A++ enrichment vs non-corridor
    corridor_set      = set(CORRIDOR_HARMONICS)
    survey_non_corr   = [n for n in SURVEY_HARMONICS if n not in corridor_set]
    n_ncorr_app       = 0
    n_ncorr_total     = len(survey_non_corr)
    for n in survey_non_corr:
        b = best_site_at_harmonic(n, all_sites)
        if b and b["tier"] == "A++":
            n_ncorr_app += 1

    ft_table = [
        [n_app,       N_HARMONICS - n_app],
        [n_ncorr_app, n_ncorr_total - n_ncorr_app],
    ]
    ft_OR, ft_p = fisher_exact(ft_table, alternative="greater")

    # Mann-Whitney precision comparison
    corr_devs   = [h["dev"] for h in corridor_hits]
    all_ap_devs = [beru_deviation(s.longitude) for s in all_sites
                   if is_aplus(tier_label(beru_deviation(s.longitude)))]
    corr_mean   = sum(corr_devs) / len(corr_devs) if corr_devs else 0.0
    global_mean = sum(all_ap_devs) / len(all_ap_devs) if all_ap_devs else 0.0
    mw_stat, mw_p = mannwhitneyu(corr_devs, all_ap_devs, alternative="less")

    lines += [
        SEP,
        "STATISTICAL TESTS",
        "",
        f"  TEST 1: Binomial occupancy",
        f"    H0: each harmonic independently occupied with P0 = {P_NULL_AP:.2f}",
        f"    Observed: {n_hits}/{N_HARMONICS} = {100*n_hits/N_HARMONICS:.0f}%  "
        f"(expected: {N_HARMONICS*P_NULL_AP:.1f})",
        f"    Binomial p (one-tailed): {p_binom:.2e}  {sig(p_binom)}",
        "",
        f"  TEST 2: Fisher combined probability (precision within corridor)",
        f"    For each occupied harmonic: P_i = dev_i / {TIER_B_MAX} (CDF under uniform null)",
        f"    X2 = {chi_obs:.2f}  (df = {df_fisher})",
        f"    Fisher combined p: {p_fisher:.2e}  {sig(p_fisher)}",
        "",
        f"  TEST 3: A++ enrichment in corridor vs non-corridor",
        f"    Corridor    ({N_HARMONICS:>2} harmonics): {n_app} A++ sites",
        f"    Non-corridor ({n_ncorr_total:>2} harmonics): {n_ncorr_app} A++ sites",
        f"    Fisher exact (one-tailed): OR = {ft_OR:.2f},  p = {ft_p:.4f}  {sig(ft_p)}",
        "",
        f"  TEST 4: Corridor precision vs global A+ mean deviation",
        f"    Corridor A+ mean dev: {corr_mean:.6f} beru  ({corr_mean*BERU*111:.2f} km)",
        f"    Global A+ mean dev:   {global_mean:.6f} beru  ({global_mean*BERU*111:.2f} km)",
        f"    Precision ratio:      {global_mean/max(corr_mean,1e-10):.2f}x",
        f"    Mann-Whitney U (corridor < global): p = {mw_p:.4f}  {sig(mw_p)}",
        "",
    ]

    # ── Anchor comparison ─────────────────────────────────────────────────────
    lines += [
        SEP,
        "ANCHOR COMPARISON: CORRIDOR OCCUPANCY (17 harmonics)",
        f"  {'Anchor':<22}  {'Lon°E':>8}  {'Hits/17':>8}  {'A++':>5}  {'Occupancy%':>11}  {'Bar (17 cells)'}",
        "  " + "-" * 85,
    ]

    anchors = _corr["anchor_comparison"]
    for a in anchors:
        hits, app_n = corridor_occupancy(a["lon"], all_sites, n_harmonics=17)
        pct = 100 * hits / 17
        bar = "#" * hits + "." * (17 - hits)
        lines.append(
            f"  {a['label']:<22}  {a['lon']:>8.4f}  {hits:>6}/17  {app_n:>5}  "
            f"{pct:>10.0f}%  [{bar}]  ({a['note']})"
        )

    # Add notable anchors for context
    lines += [
        "",
        "  Extended notable-anchor comparison (Gerizim vs 17-harmonic corridor):",
        f"  {'Anchor':<30}  {'Lon°E':>8}  {'Hits/17':>8}  {'A++':>5}",
        "  " + "-" * 60,
    ]
    notable = load_notable_anchors()
    for label, lon in notable.items():
        hits, app_n = corridor_occupancy(lon, all_sites, n_harmonics=17)
        lines.append(f"  {label:<30}  {lon:>8.4f}  {hits:>6}/17  {app_n:>5}")

    lines.append("")

    # ── Summary ───────────────────────────────────────────────────────────────
    lines += [
        SEP,
        "SUMMARY",
        "",
        f"  {n_hits}/{N_HARMONICS} harmonics in the Gerizim-Lumbini corridor contain at least",
        f"  one UNESCO World Heritage site within the A+ band (< {TIER_APLUS*BERU*111:.1f} km of harmonic).",
        f"",
        f"  Binomial p = {p_binom:.2e}  {sig(p_binom)}",
        f"  Fisher combined p = {p_fisher:.2e}  {sig(p_fisher)}",
        f"  A++ enrichment vs non-corridor: OR = {ft_OR:.2f}, p = {ft_p:.4f}  {sig(ft_p)}",
        f"  Corridor A+ sites are {global_mean/max(corr_mean,1e-10):.2f}x more precise than the global A+ mean.",
        "",
    ]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Written -> {OUT}  ({n_hits}/{N_HARMONICS} corridor harmonics occupied)")


if __name__ == "__main__":
    run()
