"""
reporting.py — Output formatting for analysis results.

All print/display logic lives here so that analysis scripts
contain only data processing and statistical calculations.
"""

from typing import List, Dict, Any, Optional, Tuple

from lib.beru import is_aplus, is_a_or_better, CONFIG
from lib.stats import (
    BinomialResult, ChiSquareResult, FisherResult,
    CochranArmitageResult, SpearmanResult,
    significance_label, bonferroni,
)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _W(width: int = 90) -> str:
    return "=" * width

def _sep(width: int = 90) -> str:
    return "─" * width


# ---------------------------------------------------------------------------
# Section headers
# ---------------------------------------------------------------------------

def print_header(title: str, width: int = 90):
    """Print a minor section header."""
    print()
    print(_W(width))
    print(f"  {title}")
    print(_W(width))


def print_subheader(title: str):
    """Print a sub-section header (no equals bars)."""
    print()
    print(f"  {title}")
    print("  " + _sep(len(title) + 2))


def print_major_header(title: str, width: int = 90):
    """Print a major section header with description lines below the bar."""
    print()
    print(_W(width))
    print(f"  {title}")
    print(_W(width))


# ---------------------------------------------------------------------------
# Configuration summary
# ---------------------------------------------------------------------------

def print_config_summary(anchor_name: str = "gerizim"):
    """Print the configuration being used, sourced from config.json."""
    anchor = CONFIG["anchors"][anchor_name]
    tiers = CONFIG["tiers"]
    print(f"Anchor: {anchor['label']} ({anchor['longitude']}°E)")
    print(f"Beru unit: {CONFIG['units']['beru']['degrees']}° of arc")
    print(f"Tier thresholds:")
    for name, tier in tiers.items():
        print(f"  {name}: ≤ {tier['max_deviation_beru']} beru "
              f"(~{tier['approx_km']} km)")


# ---------------------------------------------------------------------------
# Population summary
# ---------------------------------------------------------------------------

def print_population_summary(sites: List[Dict], excluded: List[Dict],
                             rejected: List[Dict], total: int):
    """Print how the test population was constructed."""
    print_subheader("Population")
    print(f"UNESCO corpus entries: {total}")
    print(f"Sites matching keywords: {len(sites) + len(excluded) + len(rejected)}")
    print(f"Context-rejected (ambiguous keyword, no architectural context): {len(rejected)}")
    print(f"Manually excluded: {len(excluded)}")
    print(f"Final test population: {len(sites)}")

    if any("match_source" in s for s in sites):
        n_xml = sum(1 for s in sites if s.get("match_source") == "xml")
        n_ext = sum(1 for s in sites if s.get("match_source") == "extended")
        print(f"Matched via XML text: {n_xml}")
        print(f"Matched via extended OUV text only: {n_ext}")


# ---------------------------------------------------------------------------
# Site table
# ---------------------------------------------------------------------------

def print_site_table(sites: List[Dict]):
    """Print one row per site, sorted by deviation."""
    print_subheader("Sites (sorted by deviation from nearest harmonic)")
    fmt = "{name:<50s}  {lon:>8.4f}°  {dev:>7.5f} beru  {dev_km:>6.1f} km  {tier}  [{kw}]"
    for s in sorted(sites, key=lambda x: x["dev"]):
        tier_mark = s["tier"]
        if is_aplus(tier_mark):
            tier_mark += " <--"
        print(fmt.format(
            name=s["name"][:50],
            lon=s["lon"],
            dev=s["dev"],
            dev_km=s["dev_km"],
            tier=tier_mark,
            kw=", ".join(s.get("keywords", []))[:30],
        ))


# ---------------------------------------------------------------------------
# Enrichment row table
#
# Used by: origin_sites_test, founding_enrichment_test, meta_keyword_test,
#          tumulus_dome_evolution_test, …
# ---------------------------------------------------------------------------

_ENR_HEADER_FMT = (
    "  {label:<{lw}}  {n:>5}  {nap:>4}  {rate:>6}  {enr:>7}  {p:>10}  {sig_h:<4}  {extra}"
)
_ENR_ROW_FMT = (
    "  {label:<{lw}}  {n:>5}  {nap:>4}  {rate:>5.1f}%  {enr:>6.2f}×  {p:>10.4f}  {sig_v:<4}  {extra}"
)

def print_enrichment_header(label_width: int = 30, extra_header: str = ""):
    """Print the column header for an enrichment table."""
    print(_ENR_HEADER_FMT.format(
        label="Group", lw=label_width,
        n="N", nap="A+", rate=" A+%", enr="Enrich",
        p="  p(A+)", sig_h="Sig", extra=extra_header,
    ))
    print("  " + _sep(label_width + 50))


def print_enrichment_row(label: str, n: int, n_aplus: int, p: float,
                         null_rate: float = 0.04, label_width: int = 30,
                         extra: str = ""):
    """Print one enrichment row (N, A+, rate, enrichment, p, sig)."""
    rate = n_aplus / n if n else 0.0
    enr  = rate / null_rate if null_rate else 0.0
    sig  = significance_label(p)
    print(_ENR_ROW_FMT.format(
        label=label, lw=label_width,
        n=n, nap=n_aplus, rate=100 * rate,
        enr=enr, p=p, sig_v=sig, extra=extra,
    ))


def print_enrichment_table(rows: List[Dict], null_rate: float = 0.04,
                           label_width: int = 30, extra_header: str = ""):
    """
    Print a full enrichment table from a list of row dicts.

    Each dict must have: label, n, n_aplus, p.
    Optional keys: extra (str appended to the row).
    """
    print_enrichment_header(label_width, extra_header)
    for r in rows:
        print_enrichment_row(
            label=r["label"], n=r["n"], n_aplus=r["n_aplus"],
            p=r["p"], null_rate=null_rate,
            label_width=label_width, extra=r.get("extra", ""),
        )


# ---------------------------------------------------------------------------
# Cohort table
# ---------------------------------------------------------------------------

def print_cohort_table(rows: List[Dict], null_rate: float = 0.04):
    """
    Print a cohort-level enrichment table.

    Each dict must have: label, n, n_aplus, p.
    """
    print_enrichment_table(rows, null_rate=null_rate, label_width=32)


# ---------------------------------------------------------------------------
# Category / thematic breakdown table
#
# Used by: founding_sites_analysis, founding_enrichment_test
# ---------------------------------------------------------------------------

def print_category_table(rows: List[Dict], n_corpus: int, n_aplus: int,
                         null_rate: float = 0.04, label_width: int = 25):
    """
    Print per-category enrichment (corpus rate, A+ rate, enrichment, Fisher p).

    Each dict must have: label, n_corpus, n_aplus_cat, n_nap_cat, p_fisher.
    """
    print(
        f"  {'Category':<{label_width}}  {'Corpus':>8}  {'Corpus%':>8}"
        f"  {'A+':>5}  {'A+%':>7}  {'Enrich':>7}  {'p(Fisher)':>10}  Sig"
    )
    print("  " + _sep(label_width + 65))
    for r in rows:
        n_c   = r["n_corpus"]
        n_ap  = r["n_aplus_cat"]
        rate_c = n_c / n_corpus if n_corpus else 0.0
        rate_ap = n_ap / n_aplus if n_aplus else 0.0
        enr    = rate_ap / rate_c if rate_c > 0 else float("inf")
        sig    = significance_label(r["p_fisher"])
        print(
            f"  {r['label']:<{label_width}}"
            f"  {n_c:>5}/{n_corpus}"
            f"  {100*rate_c:>6.1f}%"
            f"  {n_ap:>3}/{n_aplus}"
            f"  {100*rate_ap:>5.1f}%"
            f"  {enr:>6.2f}×"
            f"  {r['p_fisher']:>10.6f}"
            f"  {sig}"
        )


# ---------------------------------------------------------------------------
# Thematic site listing (founding_sites_analysis style)
# ---------------------------------------------------------------------------

def print_thematic_site_rows(sites: List[Dict], matched_kw_key: str = "matched_kw"):
    """
    Print one row per site: dev, km, year, name, matched keyword.

    Each dict must have: dev, km, yr, name.
    Optional: matched_kw_key (defaults to s['matched_kw']).
    """
    print(
        f"  {'Dev(beru)':<11}  {'km':>5}  {'Yr':<5}  {'Site':<50}  Matched keyword"
    )
    print("  " + _sep(100))
    for s in sorted(sites, key=lambda x: x["yr"]):
        kw = s.get(matched_kw_key, "—")
        print(
            f"  {s['dev']:<11.6f}  {s['km']:>5.1f}  {s['yr']:<5}"
            f"  {s['name'][:50]:<50}  \"{kw}\""
        )


def print_thematic_breakdown(by_cat: Dict[str, List[Dict]],
                             cat_labels: Dict[str, str],
                             kw_lists: Dict[str, List[str]],
                             priority: List[str] = None):
    """
    Print a full thematic breakdown: one section per category.

    Parameters
    ----------
    by_cat     : {cat_code: [site_dict, ...]}
    cat_labels : {cat_code: human-readable label}
    kw_lists   : {cat_code: [keyword, ...]}   — used to find the matched kw
    priority   : ordered list of cat codes (default FSMXL?)
    """
    if priority is None:
        priority = ["F", "S", "M", "X", "L", "?"]
    for cat in priority:
        items = by_cat.get(cat, [])
        if not items:
            continue
        print(f"\n  ── {cat_labels.get(cat, cat)}  (N={len(items)}) ──")
        kws = kw_lists.get(cat, [])
        print(
            f"  {'Dev(beru)':<11}  {'km':>5}  {'Yr':<5}  {'Site':<50}  Matched keyword"
        )
        print("  " + _sep(100))
        for s in sorted(items, key=lambda x: x["yr"]):
            matched_kw = next((kw for kw in kws if kw in s.get("text", "")), "—")
            print(
                f"  {s['dev']:<11.6f}  {s['km']:>5.1f}  {s['yr']:<5}"
                f"  {s['name'][:50]:<50}  \"{matched_kw}\""
            )


def print_unclassified_sites(sites: List[Dict]):
    """Print the unclassified A+ sites block for transparency."""
    if not sites:
        return
    print()
    print(_W(95))
    print("  UNCLASSIFIED A+ SITES (no config.json keyword matched)")
    print("  Shown for transparency — these are NOT excluded from any test.")
    print("  Extend the keyword lists in config.json if any should be classified.")
    print(_W(95))
    print(f"  {'Dev(beru)':<11}  {'km':>5}  {'Yr':<5}  Site")
    print("  " + _sep(70))
    for s in sorted(sites, key=lambda x: x["yr"]):
        print(f"  {s['dev']:<11.6f}  {s['km']:>5.1f}  {s['yr']:<5}  {s['name']}")


# ---------------------------------------------------------------------------
# Top-N precision table
# ---------------------------------------------------------------------------

def print_top_sites(sites: List[Dict], n: int = 10, title: str = None):
    """Print the N most precise A+ hits (smallest deviation)."""
    title = title or f"TOP {n} MOST PRECISE HITS (smallest deviation from harmonic)"
    print_header(title, width=95)
    print(f"\n  {'Dev (beru)':<12}  {'km':>5}  {'Yr':<5}  {'Cats':<6}  Site")
    print("  " + _sep(85))
    for s in sorted(sites, key=lambda x: x["dev"])[:n]:
        cats_str = "".join(sorted(s.get("cats", set()))) if s.get("cats") else "?"
        print(
            f"  {s['dev']:<12.6f}  {s['km']:>5.1f}  {s['yr']:<5}"
            f"  {cats_str:<6}  {s['name'][:60]}"
        )


# ---------------------------------------------------------------------------
# Statistical results
# ---------------------------------------------------------------------------

def print_binomial_results(label: str, result: BinomialResult,
                           bonferroni_k: int = None):
    """Print one binomial test result."""
    sig = significance_label(result.p_value)
    line = (
        f"  {label}: {result.observed}/{result.total} = "
        f"{result.observed_rate:.1%}  "
        f"(expected {result.expected:.1f}, "
        f"enrichment {result.ratio:.2f}x, "
        f"p = {result.p_value:.4f} {sig})"
    )
    if bonferroni_k:
        adj = bonferroni(result.p_value, bonferroni_k)
        line += f"  [Bonferroni k={bonferroni_k}: p_adj = {adj:.4f}]"
    print(line)


def print_chi_square_result(result: ChiSquareResult):
    """Print chi-square uniformity test result."""
    sig = significance_label(result.p_value)
    print(f"  Chi-square uniformity ({result.n_bins} bins): "
          f"statistic = {result.statistic:.2f}, "
          f"p = {result.p_value:.4f} {sig}")
    print(f"  Bin counts: {result.observed_bins}")


def print_fisher_result(label: str, result: FisherResult):
    """Print a Fisher exact test result."""
    sig = significance_label(result.p_value)
    print(f"  {label}: OR = {result.odds_ratio:.2f}, "
          f"p = {result.p_value:.6f} {sig}")


def print_fisher_inline(label: str, or_val: float, p: float):
    """Print a Fisher exact test result from raw values (no FisherResult object)."""
    sig = significance_label(p)
    print(f"  {label}: OR = {or_val:.2f},  p = {p:.4f}  {sig}")


def print_binomial_inline(label: str, k: int, n: int, null_rate: float, p: float):
    """Print a binomial test result from raw values."""
    sig = significance_label(p)
    rate = 100 * k / n if n else 0.0
    print(
        f"  {label}: {k}/{n} = {rate:.1f}%  "
        f"vs null {100*null_rate:.1f}%  →  p = {p:.6f}  {sig}"
    )


def print_cochran_armitage_result(result: CochranArmitageResult):
    """Print a Cochran-Armitage trend test result."""
    sig = significance_label(result.p_value)
    print(f"  Cochran-Armitage Z = {result.z_statistic:.4f}, "
          f"p = {result.p_value:.6f} {sig}")
    print(f"  Direction: {result.direction}")


def print_spearman_result(label: str, result: SpearmanResult):
    """Print a Spearman correlation result."""
    sig = significance_label(result.p_value)
    print(f"  {label}: rho = {result.rho:.4f}, "
          f"p = {result.p_value:.6f} {sig}")


# ---------------------------------------------------------------------------
# Tier counts
# ---------------------------------------------------------------------------

def print_tier_counts(sites: List[Dict]):
    """Print how many sites fall in each tier."""
    from collections import Counter
    tiers = Counter(s["tier"] for s in sites)
    print_subheader("Tier distribution")
    for tier in ["A++", "A+", "A", "B", "C"]:
        count = tiers.get(tier, 0)
        pct = 100 * count / len(sites) if sites else 0
        print(f"  {tier}: {count} ({pct:.1f}%)")


# ---------------------------------------------------------------------------
# Hit lists
# ---------------------------------------------------------------------------

def print_hit_list(label: str, sites: List[Dict]):
    """Print a list of notable sites with their UNESCO descriptions."""
    print_subheader(label)
    for s in sorted(sites, key=lambda x: x["dev"]):
        print(f"  {s['name'][:60]}")
        print(f"    deviation: {s['dev']:.5f} beru ({s['dev_km']:.1f} km), "
              f"lon: {s['lon']:.4f}°")
        if s.get("sentences"):
            print(f"    UNESCO: \"{s['sentences'][0][:160]}\"")
        print()


# ---------------------------------------------------------------------------
# Anchor sweep table
# ---------------------------------------------------------------------------

def print_anchor_sweep_result(anchors: Dict[str, Dict]):
    """
    Print the anchor sweep comparison table.

    Parameters
    ----------
    anchors : dict
        {label: {count_aplus, pctile_aplus, count_a, pctile_a}}
    """
    print_subheader("Anchor sweep comparison")
    print(f"  {'Anchor':<26s}  {'A+ count':>8}  {'Pctile':>7}  "
          f"{'A count':>8}  {'Pctile':>7}")
    for label, data in anchors.items():
        print(f"  {label:<26s}  {data['count_aplus']:>8}  "
              f"{data['pctile_aplus']:>6.0f}th  "
              f"{data['count_a']:>8}  {data['pctile_a']:>6.0f}th")


# ---------------------------------------------------------------------------
# Keyword composition
# ---------------------------------------------------------------------------

def print_keyword_composition(sites: List[Dict]):
    """Print which keywords matched how many sites."""
    from collections import Counter
    counts = Counter()
    for s in sites:
        for kw in s.get("keywords", []):
            counts[kw] += 1
    print_subheader("Keyword composition")
    for kw, count in counts.most_common():
        print(f"  {kw}: {count} sites")


# ---------------------------------------------------------------------------
# Keyword enrichment table (meta_keyword_test style)
# ---------------------------------------------------------------------------

def print_keyword_enrichment_table(results: List[Dict], n_show: int = 30,
                                   title: str = "Top keywords by χ²-uniform p"):
    """
    Print a ranked keyword enrichment table.

    Each dict must have: word, n, n_a, rate_a, p_a, chi_p.
    """
    print(f"\n  {title}:")
    print(f"  {'Word':<25s}  {'N':>5}  {'Tier-A':>6}  {'Rate-A':>7}  {'p_A':>8}  {'χ²p':>8}")
    print(f"  {'─'*70}")
    for r in results[:n_show]:
        sig = significance_label(r.get("chi_p", r.get("p_a", 1.0)))
        print(
            f"  {r['word']:<25s}  {r['n']:>5d}  {r['n_a']:>6d}"
            f"  {100*r['rate_a']:>5.1f}%  {r['p_a']:>8.4f}  {r.get('chi_p', float('nan')):>8.4f}"
            f"  {sig}"
        )


# ---------------------------------------------------------------------------
# Summary verdict
# ---------------------------------------------------------------------------

def print_verdict(n: int, result_aplus: BinomialResult, result_a: BinomialResult,
                  chi_result: ChiSquareResult, bonferroni_k: int = 5):
    """Print a concise summary of all results."""
    print_subheader("Summary")
    print(f"Population: {n} sites")
    print()
    print("Primary test (Tier-A+):")
    print_binomial_results("Tier-A+", result_aplus, bonferroni_k)
    print()
    print("Secondary test (Tier-A):")
    print_binomial_results("Tier-A", result_a, bonferroni_k)
    print()
    print("Uniformity test:")
    print_chi_square_result(chi_result)


# ---------------------------------------------------------------------------
# LaTeX macros
# ---------------------------------------------------------------------------

def print_latex_macros(macros: Dict[str, Any]):
    """
    Print LaTeX \\newcommand macros for inclusion in the manuscript.

    Parameters
    ----------
    macros : dict
        {macro_name: value}. Values are formatted automatically:
        int → plain, float → .4f, str → as-is.
    """
    print_header("LATEX MACROS")
    for name, value in macros.items():
        if isinstance(value, int):
            formatted = str(value)
        elif isinstance(value, float):
            if abs(value) < 0.0001 and value != 0:
                formatted = f"{value:.1e}"
            elif abs(value) >= 100:
                formatted = f"{value:.1f}"
            else:
                formatted = f"{value:.4f}"
        else:
            formatted = str(value)
        print(f"  \\newcommand{{\\{name}}}{{{formatted}}}")
