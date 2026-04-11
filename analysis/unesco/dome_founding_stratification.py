"""
dome_founding_stratification.py
================================
Cross-tabulate the Test-2 dome population (spherical_monument_raw_sweep)
against the founding/sacred classifier (lib/founding_filter.py) and
compare Tier-A+ rates for:
  - dome sites that are also F/S (founding capital or sacred origin)
  - dome sites that are not F/S

Tests whether dome enrichment is separable from prestige-site concentration.
"""

import sys
import re
from pathlib import Path
from scipy.stats import binomtest, fisher_exact

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus
from lib.beru import GERIZIM, BERU, TIER_APLUS, P_NULL_AP, deviation as beru_dev, is_aplus, tier_label
from lib.dome_filter import FORM_KEYWORDS, FORM_KEYWORD_RES, validate_keyword_match
from lib.founding_filter import classify_site as founding_classify

def get_dome_sites():
    """Return all dome-population sites (raw keyword sweep, Test 2)."""
    corpus = load_corpus()
    sites = []
    for site_obj in corpus:
        if site_obj.category == "Natural" or not site_obj.has_coords:
            continue
        full_text = site_obj.full_text
        all_text = (site_obj.short_description or "") + " " + (site_obj.extended_description or "")
        raw_matched = [k for k in FORM_KEYWORDS if FORM_KEYWORD_RES[k].search(full_text)]
        if not raw_matched:
            continue
        validated_keywords = [k for k in raw_matched
                               if validate_keyword_match(all_text, k)[0]]
        if not validated_keywords:
            continue
        dev = beru_dev(site_obj.longitude)
        ap  = is_aplus(tier_label(dev))
        founding_cats = set(founding_classify(site_obj).keys())
        is_fs = bool(founding_cats & {"F", "S"})
        sites.append({
            "name": site_obj.site,
            "lon":  site_obj.longitude,
            "dev":  dev,
            "ap":   ap,
            "founding_cats": founding_cats,
            "is_fs": is_fs,
        })
    return sites

def binomial_row(label, n, k, p0=P_NULL_AP):
    if n == 0:
        return f"  {label}: n=0"
    result = binomtest(k, n, p0, alternative="greater")
    rate = 100 * k / n
    enrich = rate / (100 * p0)
    return (f"  {label}: n={n}, A+={k}, rate={rate:.1f}%, "
            f"enrich={enrich:.2f}x, p={result.pvalue:.4f}")

def main():
    sites = get_dome_sites()
    total = len(sites)

    fs     = [s for s in sites if s["is_fs"]]
    non_fs = [s for s in sites if not s["is_fs"]]

    n_fs,     n_non_fs     = len(fs),     len(non_fs)
    fs_ap,    non_fs_ap    = sum(s["ap"] for s in fs), sum(s["ap"] for s in non_fs)
    fs_rate   = 100 * fs_ap    / n_fs     if n_fs     else 0.0
    nfs_rate  = 100 * non_fs_ap / n_non_fs if n_non_fs else 0.0
    fs_enrich = (fs_rate  / (100 * P_NULL_AP)) if n_fs     else 0.0
    nfs_enrich= (nfs_rate / (100 * P_NULL_AP)) if n_non_fs else 0.0

    bt_fs  = binomtest(fs_ap,     n_fs,     P_NULL_AP, alternative="greater")
    bt_nfs = binomtest(non_fs_ap, n_non_fs, P_NULL_AP, alternative="greater")

    table = [[fs_ap,     n_fs     - fs_ap],
             [non_fs_ap, n_non_fs - non_fs_ap]]
    odds, p_fish = fisher_exact(table, alternative="greater")

    # ── Human-readable summary ──────────────────────────────────────────────
    print(f"\nDome × Founding-category stratification")
    print(f"  Total dome population (Test 2, raw sweep): n={total}")
    print()
    print(binomial_row("F/S dome sites     (founding capital or sacred origin)", n_fs,     fs_ap))
    print(binomial_row("non-F/S dome sites (all other categories)",              n_non_fs, non_fs_ap))
    print()
    print(f"  Fisher exact (F/S vs non-F/S, A+ rate): OR={odds:.2f}, p={p_fish:.4f}")
    print()
    print("  Non-F/S dome A+ sites:")
    for s in non_fs:
        if s["ap"]:
            cats = s["founding_cats"] if s["founding_cats"] else {"?"}
            print(f"    {s['name'][:55]:<55}  lon={s['lon']:.3f}  dev={s['dev']:.4f}  cats={cats}")

    # ── LaTeX macros ────────────────────────────────────────────────────────
    print(f"  \\newcommand{{\\domeStratNfs}}{{{n_fs}}}                % dome sites that are F/S")
    print(f"  \\newcommand{{\\domeStratNnfs}}{{{n_non_fs}}}              % dome sites that are not F/S")
    print(f"  \\newcommand{{\\domeStratFsAp}}{{{fs_ap}}}               % A+ count, F/S dome stratum")
    print(f"  \\newcommand{{\\domeStratNfsAp}}{{{non_fs_ap}}}            % A+ count, non-F/S dome stratum")
    print(f"  \\newcommand{{\\domeStratFsRate}}{{{fs_rate:.1f}}}         % A+ rate (%), F/S stratum")
    print(f"  \\newcommand{{\\domeStratNfsRate}}{{{nfs_rate:.1f}}}       % A+ rate (%), non-F/S stratum")
    print(f"  \\newcommand{{\\domeStratFsEnrich}}{{{fs_enrich:.2f}}}     % enrichment, F/S stratum")
    print(f"  \\newcommand{{\\domeStratNfsEnrich}}{{{nfs_enrich:.2f}}}   % enrichment, non-F/S stratum")
    print(f"  \\newcommand{{\\domeStratFsP}}{{{bt_fs.pvalue:.4f}}}       % binomial p, F/S stratum")
    print(f"  \\newcommand{{\\domeStratNfsP}}{{{bt_nfs.pvalue:.4f}}}     % binomial p, non-F/S stratum")
    print(f"  \\newcommand{{\\domeStratFisherOR}}{{{odds:.2f}}}          % Fisher OR, F/S vs non-F/S")
    print(f"  \\newcommand{{\\domeStratFisherP}}{{{p_fish:.4f}}}         % Fisher p, F/S vs non-F/S")

if __name__ == "__main__":
    main()
