r"""
keyword_sensitivity_stupa.py
=============================
Keyword-stratified sensitivity check for Test 2 (domed monuments).

Partitions the dome population (N=\NcircTotal, raw keyword sweep) by the
keyword that triggered inclusion:

  Stage     Keywords
  -------   -----------------------------------------
  Mound     tumulus, barrow, kofun
  Stupa     stupa, stupas
  Dome      dome, domed, domes, tholos, spherical

(A site matched by multiple stages is counted in its most specific stage:
 stupa > dome > mound, to prevent double-counting.)

For each stage and the combined population, computes A++/A+/A counts,
rates, Fisher OR and p vs. full corpus null, and emits LaTeX macros:

  \kwSensNMound, \kwSensNStupa, \kwSensNDome
  \kwSensAppNMound  \kwSensApNMound  \kwSensANMound
  \kwSensAppNStupa  \kwSensApNStupa  \kwSensANStupa
  \kwSensAppNDome   \kwSensApNDome   \kwSensANDome
  \kwSensAppNCombined ...

  \kwSensAppORMound  \kwSensAppPMound  (etc.)

The key "leave-stupa-out" macros are:
  \\kwSensNoStupaN          -- n for dome+mound without stupa
  \\kwSensNoStupaAppN       -- A++ count
  \\kwSensNoStupaAppRate    -- A++ rate (%)
  \\kwSensNoStupaAppOR      -- Fisher OR vs null
  \\kwSensNoStupaAppP       -- Fisher p
  \\kwSensNoStupaApN        -- A+ count
  \\kwSensNoStupaApRate     -- A+ rate (%)
  \\kwSensNoStupaApOR
  \\kwSensNoStupaApP

Run from repo root:
    python3 analysis/unesco/keyword_sensitivity_stupa.py
"""

import sys
import re
from pathlib import Path
from scipy.stats import fisher_exact, binomtest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from data.unesco_corpus import load_corpus
from lib.beru import (
    GERIZIM, TIER_APP, TIER_APLUS, TIER_A_MAX,
    P_NULL_APP, P_NULL_AP, P_NULL_A,
    deviation as beru_dev, tier_label,
)
from lib.dome_filter import FORM_KEYWORDS, FORM_KEYWORD_RES, validate_keyword_match
from lib.results_store import ResultsStore

# ── Keyword stage definitions ─────────────────────────────────────────────────
MOUND_KW   = {"tumulus", "barrow", "kofun"}
STUPA_KW   = {"stupa", "stupas"}
DOME_KW    = {"dome", "domed", "domes", "tholos", "spherical"}

# Full corpus size (for Fisher contingency table)
def _full_corpus_n():
    corpus = load_corpus()
    return sum(
        1 for s in corpus
        if s.category != "Natural" and s.has_coords
    )

def _classify_site(site_obj):
    """Return (stage, dev) or None if site is not in dome population."""
    full_text = site_obj.full_text
    all_text = (site_obj.short_description or "") + " " + (site_obj.extended_description or "")
    raw_matched = {k for k in FORM_KEYWORDS if FORM_KEYWORD_RES[k].search(full_text)}
    if not raw_matched:
        return None
    validated = {k for k in raw_matched if validate_keyword_match(all_text, k)[0]}
    if not validated:
        return None
    # Assign to most specific stage (stupa > dome > mound)
    if validated & STUPA_KW:
        stage = "stupa"
    elif validated & DOME_KW:
        stage = "dome"
    elif validated & MOUND_KW:
        stage = "mound"
    else:
        stage = "dome"  # fallback
    dev = beru_dev(site_obj.longitude)
    return stage, dev

def _fisher_or_p(k, n, p_null, N_total):
    """One-sided Fisher exact (greater) vs. full-corpus background."""
    # Contingency: [[k, n-k], [N_total*p_null, N_total*(1-p_null)]]
    # Use integer approximation for background cell
    bg_hit   = round(N_total * p_null)
    bg_miss  = N_total - bg_hit
    table = [[k, n - k], [bg_hit, bg_miss]]
    _, p = fisher_exact(table, alternative="greater")
    if n == 0:
        return float("nan"), float("nan")
    obs_rate = k / n
    if bg_hit == 0:
        OR = float("inf")
    else:
        OR = (k / max(n - k, 1)) / (bg_hit / max(bg_miss, 1))
    return round(OR, 2), p

def _fmt_p(p):
    if p < 0.001:
        return "< 0.001"
    return f"{p:.4f}"

def main():
    corpus = load_corpus()
    N_total = _full_corpus_n()

    sites = []
    for s in corpus:
        if s.category == "Natural" or not s.has_coords:
            continue
        result = _classify_site(s)
        if result is None:
            continue
        stage, dev = result
        t = tier_label(dev)
        sites.append({"stage": stage, "dev": dev, "tier": t})

    stages = ["mound", "stupa", "dome"]
    stage_sites = {st: [s for s in sites if s["stage"] == st] for st in stages}
    combined    = sites
    no_stupa    = [s for s in sites if s["stage"] != "stupa"]

    store = ResultsStore()

    def _emit_stage(label, subset, prefix):
        n    = len(subset)
        n_app = sum(1 for s in subset if s["dev"] <= TIER_APP)
        n_ap  = sum(1 for s in subset if s["dev"] <= TIER_APLUS)
        n_a   = sum(1 for s in subset if s["dev"] <= TIER_A_MAX)

        data = {
            f"{prefix}N":       n,
            f"{prefix}AppN":    n_app,
            f"{prefix}ApN":     n_ap,
            f"{prefix}AN":      n_a,
        }

        if n > 0:
            data[f"{prefix}AppRate"] = round(100 * n_app / n, 1)
            data[f"{prefix}ApRate"]  = round(100 * n_ap  / n, 1)
            data[f"{prefix}ARate"]   = round(100 * n_a   / n, 1)
            OR_app, p_app = _fisher_or_p(n_app, n, P_NULL_APP, N_total)
            OR_ap,  p_ap  = _fisher_or_p(n_ap,  n, P_NULL_AP,  N_total)
            OR_a,   p_a   = _fisher_or_p(n_a,   n, P_NULL_A,   N_total)
            data[f"{prefix}AppOR"] = OR_app
            data[f"{prefix}AppP"]  = _fmt_p(p_app)
            data[f"{prefix}ApOR"]  = OR_ap
            data[f"{prefix}ApP"]   = _fmt_p(p_ap)
            data[f"{prefix}AOR"]   = OR_a
            data[f"{prefix}AP"]    = _fmt_p(p_a)

        store.write_many(data)

        app_rate = f"{100*n_app/n:.1f}%" if n > 0 else "n/a"
        ap_rate  = f"{100*n_ap/n:.1f}%"  if n > 0 else "n/a"
        print(f"\n{label} (n={n})")
        print(f"  A++: {n_app} ({app_rate}) OR={data.get(f'{prefix}AppOR','--')} p={data.get(f'{prefix}AppP','--')}")
        print(f"  A+:  {n_ap}  ({ap_rate})  OR={data.get(f'{prefix}ApOR','--')}  p={data.get(f'{prefix}ApP','--')}")
        print(f"  A:   {n_a}   ({f'{100*n_a/n:.1f}%' if n>0 else 'n/a'})")

    for st, label, prefix in [
        ("mound", "Mound",    "kwSensMound"),
        ("stupa", "Stupa",    "kwSensStupa"),
        ("dome",  "Dome",     "kwSensDome"),
    ]:
        _emit_stage(label, stage_sites[st], prefix)

    _emit_stage("Combined", combined, "kwSensCombined")
    _emit_stage("No-stupa (dome + mound)", no_stupa, "kwSensNoStupa")

    store.write_many({})  # flush / no-op; write_many already persists

    # ── Print \newcommand lines for the macro pipeline ────────────────────────
    for key, val in sorted(store.all().items()):
        if key.startswith("kwSens"):
            print(f"\\newcommand{{\\{key}}}{{{val}}}")

    ns  = store.read("kwSensNoStupaN", "?")
    na = store.read("kwSensNoStupaAppN", "?")
    nr = store.read("kwSensNoStupaAppRate", "?")
    np_ = store.read("kwSensNoStupaAppP", "?")
    print(f"\nMacros written. Leave-stupa-out: n={ns}, A++={na} ({nr}%), p={np_}")

if __name__ == "__main__":
    main()
