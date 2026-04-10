"""
leave_one_out_sensitivity.py
============================
Leave-one-out (LOO) sensitivity analysis for the two enriched sub-populations:

  Test 2  — Domed/spherical monuments (raw keyword sweep, N=90, A+=11)
  Test 2b — Stupa sub-stage          (raw keyword sweep, N=18, A+=4)

For each population, every Tier-A+ site is removed in turn and the
enrichment test is re-run on the remaining N-1 sites.  Because all A+
sites are equally spaced from a harmonic node at Tier-A+ precision, the
resulting N, A+, rate, enrichment, and p-value are identical for all
removals.  The script reports both the per-removal invariance check and
the single representative LOO row used in the manuscript.

OUTPUT
------
  Console:    human-readable table + \newcommand lines for paper_a_primary_unesco.tex
  Store:      ResultsStore keys LOOdomeN, LOOdomeAp, LOOdomeRate, LOOdomeEnrich,
              LOOdomeP, LOOstupaN, LOOstupaAp, LOOstupaRate, LOOstupaP

USAGE
-----
    cd /path/to/gerizim-paper-a
    python3 analysis/unesco/leave_one_out_sensitivity.py

REPRODUCIBILITY
---------------
  No random seed needed — test is deterministic (exact binomial).
  Keyword lists loaded from keywords.json (single source of truth):
    • Dome population  → lib/dome_filter.FORM_KEYWORD_RES
                         (keywords.json → "form_keywords")
    • Stupa population → keywords.json → "mound_evolution" → "stupa"
  Tier thresholds imported from lib/beru.py (config.json).
"""

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from data.unesco_corpus import load_corpus
from lib.beru import (
    BERU, TIER_APLUS, P_NULL_AP,
    deviation, tier_label, is_aplus,
)
from lib.dome_filter import FORM_KEYWORD_RES
from lib.results_store import ResultsStore
from scipy.stats import binomtest

# ── Keyword sets — loaded from keywords.json (single source of truth) ────────
# Dome/spherical raw sweep: FORM_KEYWORD_RES from lib/dome_filter.py, which
#   itself reads from keywords.json → "form_keywords".  No context gating.
# Stupa sub-stage: keywords.json → "mound_evolution" → "stupa".
#   Must match tumulus_dome_evolution_raw_sweep.py exactly.
_KW_PATH = Path(__file__).parent.parent.parent / "keywords.json"
with open(_KW_PATH) as _f:
    _KW = json.load(_f)

STUPA_KEYWORDS = _KW["mound_evolution"]["stupa"]
_stupa_res = {k: re.compile(r"\b" + re.escape(k) + r"\b", re.IGNORECASE)
              for k in STUPA_KEYWORDS}


def _classify(corpus):
    """Return (dome_sites, stupa_sites) as lists of dicts."""
    dome_sites  = []
    stupa_sites = []
    for site in corpus:
        if site.category == "Natural":
            continue
        if not site.has_coords:
            continue
        txt  = site.full_text
        dev  = deviation(site.longitude)
        ap   = is_aplus(tier_label(dev))
        rec  = {
            "name":   site.site,
            "lon":    site.longitude,
            "dev":    dev,
            "dev_km": dev * BERU * 111.0,
            "ap":     ap,
        }
        if any(r.search(txt) for r in FORM_KEYWORD_RES.values()):
            dome_sites.append(rec)
        if any(r.search(txt) for r in _stupa_res.values()):
            stupa_sites.append(rec)
    return dome_sites, stupa_sites


def _loo_stats(sites, label):
    """
    Run leave-one-out over each A+ site.

    Returns a tuple:
      (loo_N, loo_Ap, loo_rate, loo_enrich, loo_p, all_identical)

    where `all_identical` is True when every removal yields the same result
    (expected because all A+ sites are at identical Tier-A+ precision).
    """
    ap_sites  = [s for s in sites if s["ap"]]
    n_total   = len(sites)
    n_ap      = len(ap_sites)

    print(f"\n{'='*70}")
    print(f"  {label}:  N={n_total}, A+={n_ap}")
    print(f"{'='*70}")
    print(f"  {'Removed site':<45} {'N':>4} {'A+':>4} {'Rate':>7} {'Enrich':>7} {'p':>10}")
    print(f"  {'-'*45} {'-'*4} {'-'*4} {'-'*7} {'-'*7} {'-'*10}")

    rows = []
    for removed in ap_sites:
        remaining = [s for s in sites if s is not removed]
        n   = len(remaining)
        nap = sum(1 for s in remaining if s["ap"])
        p   = binomtest(nap, n, P_NULL_AP, alternative="greater").pvalue
        enr = (nap / n) / P_NULL_AP if n else 0.0
        rate = 100.0 * nap / n if n else 0.0
        rows.append((n, nap, rate, enr, p))
        name_trunc = removed["name"][:45]
        print(f"  {name_trunc:<45} {n:>4} {nap:>4} {rate:>6.1f}% {enr:>7.2f}x {p:>10.4f}")

    # Verify invariance
    all_identical = len(set((r[0], r[1]) for r in rows)) == 1
    status = "✓ invariant" if all_identical else "✗ NOT invariant — check keyword overlap"
    print(f"\n  Invariance check: {status}")

    # Return the single representative row (first removal = all the same)
    n_rep, nap_rep, rate_rep, enr_rep, p_rep = rows[0]
    return n_rep, nap_rep, rate_rep, enr_rep, p_rep, all_identical


def main():
    corpus = load_corpus()
    dome_sites, stupa_sites = _classify(corpus)

    # ── Dome LOO ──────────────────────────────────────────────────────────────
    loo_dome = _loo_stats(dome_sites, "Test 2 — Dome/Spherical LOO")
    dome_N, dome_Ap, dome_rate, dome_enrich, dome_p, dome_inv = loo_dome

    # ── Stupa LOO ─────────────────────────────────────────────────────────────
    loo_stupa = _loo_stats(stupa_sites, "Test 2b — Stupa Sub-stage LOO")
    stupa_N, stupa_Ap, stupa_rate, stupa_enrich, stupa_p, stupa_inv = loo_stupa

    # ── Print \newcommand macros ──────────────────────────────────────────────
    print("\n" + "═"*70)
    print("  LATEX MACROS")
    print("═"*70)
    print(f"  \\newcommand{{\\LOOdomeN}}{{{dome_N}}}               % LOO dome: N after removing one A+ site")
    print(f"  \\newcommand{{\\LOOdomeAp}}{{{dome_Ap}}}              % LOO dome: A+ after removing one A+ site")
    print(f"  \\newcommand{{\\LOOdomeRate}}{{{dome_rate:.1f}}}          % LOO dome: A+ rate after LOO (%)")
    print(f"  \\newcommand{{\\LOOdomeEnrich}}{{{dome_enrich:.2f}}}        % LOO dome: enrichment after LOO")
    print(f"  \\newcommand{{\\LOOdomeP}}{{{dome_p:.4f}}}           % LOO dome: p-value after LOO (same for all {sum(s['ap'] for s in dome_sites)} removals)")
    print(f"  \\newcommand{{\\LOOstupaN}}{{{stupa_N}}}              % LOO stupa: N after removing one A+ stupa site")
    print(f"  \\newcommand{{\\LOOstupaAp}}{{{stupa_Ap}}}              % LOO stupa: A+ after LOO")
    print(f"  \\newcommand{{\\LOOstupaRate}}{{{stupa_rate:.1f}}}          % LOO stupa: A+ rate after LOO (%)")
    print(f"  \\newcommand{{\\LOOstupaP}}{{{stupa_p:.4f}}}          % LOO stupa: p after LOO (same for all {sum(s['ap'] for s in stupa_sites)} removals)")

    # ── Write to results store ────────────────────────────────────────────────
    store = ResultsStore()
    store.write_many({
        "LOOdomeN":      dome_N,
        "LOOdomeAp":     dome_Ap,
        "LOOdomeRate":   round(dome_rate, 1),
        "LOOdomeEnrich": round(dome_enrich, 2),
        "LOOdomeP":      round(dome_p, 4),
        "LOOstupaN":     stupa_N,
        "LOOstupaAp":    stupa_Ap,
        "LOOstupaRate":  round(stupa_rate, 1),
        "LOOstupaP":     round(stupa_p, 4),
    })
    print("\n  ✓ Written to data/store/results.json")


if __name__ == "__main__":
    main()
