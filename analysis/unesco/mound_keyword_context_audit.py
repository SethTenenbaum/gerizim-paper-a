"""
mound_keyword_context_audit.py
================================
Context audit for the hemispherical mound-evolution keyword sweep.

PURPOSE
-------
Determines the false-positive rate for the mound/stupa/barrow/kofun keyword
set used in tumulus_dome_evolution_*.py (Test 2b).

For each Cultural/Mixed site whose full text contains any mound-evolution
keyword, this script:

  1. Shows the matching keyword(s).
  2. For the "mound" keyword (the one ambiguous term in the set), determines
     whether each matching sentence contains an archaeological context term
     (accepted) or lacks one (potentially false positive).
  3. Classifies each site as:
       ACCEPTED_UNAMB  — matched only unambiguous keywords (tumulus/tumuli/
                         barrow/barrows/kofun/stupa/stupas/dagoba/chorten/tholos)
       ACCEPTED_AMB    — matched "mound" and the sentence contains an
                         archaeological context term
       REJECTED_AMB    — matched "mound" and no matching sentence has an
                         archaeological context term (potential FP)
       ACCEPTED_DOME   — matched dome/domed/domes/spherical and the sentence
                         passes the dome_filter context check
       REJECTED_DOME   — matched a dome ambiguous keyword but failed context
  4. Computes:
       N_raw             — total keyword matches (no validation)
       N_rejected        — sites rejected after context validation
       FP_rate           — N_rejected / N_raw (upper-bound false-positive rate)
       Signal comparison — Tier-A+ counts for raw vs. validated populations

RELATIONSHIP TO OTHER SCRIPTS
------------------------------
  spherical_monument_raw_sweep.py — analogous script for the dome-only corpus
  tumulus_dome_evolution_raw_sweep.py  — raw sweep without context audit detail
  tumulus_dome_evolution_test.py       — context-validated (exploratory) version

USAGE
-----
    cd /path/to/gerizim-paper-a
    python3 analysis/unesco/mound_keyword_context_audit.py
    python3 analysis/unesco/mound_keyword_context_audit.py --verbose
    python3 analysis/unesco/mound_keyword_context_audit.py --latex

    --verbose  : print every matching sentence for each site
    --latex    : print LaTeX \\newcommand macros for the manuscript

OUTPUT
------
    Console report showing per-site classification, population-level FP rate,
    and statistical comparison between raw and context-validated populations.
"""

import re
import sys
import json
import numpy as np
from pathlib import Path
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_A_MAX,
    P_NULL_AP, P_NULL_A,
    deviation as _beru_dev, tier_label, is_aplus, is_a_or_better,
)
from lib.stats import significance_label as sig
from lib.dome_filter import (
    UNAMBIGUOUS_KEYWORDS  as DOME_UNAMBIGUOUS,
    AMBIGUOUS_KEYWORDS    as DOME_AMBIGUOUS,
    validate_keyword_match as dome_validate,
)

VERBOSE  = "--verbose" in sys.argv
LATEX    = "--latex"   in sys.argv
# When run without --latex (e.g. during reproduce_all_macros.sh --macros-only
# which greps for \newcommand lines), macros are always emitted.
# --latex additionally suppresses the large per-site tables and replaces them
# with a compact macro-only output.
MACROS_ONLY = LATEX

# ── Keyword sets ──────────────────────────────────────────────────────────────
_KW_PATH = Path(__file__).parent.parent.parent / "keywords.json"
with open(_KW_PATH) as _f:
    _KW = json.load(_f)

_evo = _KW["mound_evolution"]

MOUND_UNAMB      = _evo["mound_unambiguous"]   # tumulus, tumuli, barrow, barrows, kofun
MOUND_AMB        = _evo["mound_ambiguous"]      # mound
STUPA_KWS        = _evo["stupa"]                # stupa, stupas, dagoba, chorten
DOME_KWS         = _evo["dome"]                 # tholos, dome, domed, domes, spherical

ALL_KEYWORDS = MOUND_UNAMB + MOUND_AMB + STUPA_KWS + DOME_KWS

KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    for kw in ALL_KEYWORDS
}

# "mound" archaeological context patterns (from keywords.json)
MOUND_POS_CONTEXT = _evo["mound_positive_context"]
MOUND_POS_RES     = [re.compile(r"\b" + re.escape(p) + r"\b", re.IGNORECASE)
                     for p in MOUND_POS_CONTEXT]

def _mound_sentence_ok(sentence: str) -> bool:
    """Accept a 'mound' sentence if it contains an archaeological context term."""
    return any(pat.search(sentence) for pat in MOUND_POS_RES)

# ── Validation logic ──────────────────────────────────────────────────────────
VERDICT_UNAMB        = "ACCEPTED (unambiguous keyword)"
VERDICT_MOUND_ARCH   = "ACCEPTED (mound + arch. context)"
VERDICT_MOUND_FP     = "REJECTED (mound — no arch. context)"
VERDICT_DOME_VALID   = "ACCEPTED (dome/spherical — context validated)"
VERDICT_DOME_REJECT  = "REJECTED (dome/spherical — context failed)"


def audit_site(site_obj):
    """
    Classify a site's keyword matches.

    Returns a dict with:
      raw_kws      : list of all keywords that matched
      verdicts     : dict kw -> verdict string
      accepted_kws : list of kws that passed
      rejected_kws : list of kws that failed
      accepted      : True if ANY keyword passed validation
      mound_sentences_ok   : list of validated 'mound' sentences
      mound_sentences_bad  : list of rejected 'mound' sentences
    """
    full = site_obj.full_text  # already lowercased

    raw_kws = [kw for kw in ALL_KEYWORDS if KEYWORD_RES[kw].search(full)]
    if not raw_kws:
        return None

    all_text = (site_obj.short_description or "") + " " + (site_obj.extended_description or "")
    sentences = re.split(r"(?<=[.!?])\s+", all_text)

    verdicts             = {}
    accepted_kws         = []
    rejected_kws         = []
    mound_sentences_ok   = []
    mound_sentences_bad  = []

    for kw in raw_kws:
        kw_re = KEYWORD_RES[kw]

        if kw in MOUND_UNAMB or kw in STUPA_KWS:
            # Always unambiguous in an archaeological context
            verdicts[kw] = VERDICT_UNAMB
            accepted_kws.append(kw)

        elif kw in MOUND_AMB:
            # "mound" — needs an archaeological context sentence
            ok_sents  = [s.strip() for s in sentences if kw_re.search(s) and _mound_sentence_ok(s)]
            bad_sents = [s.strip() for s in sentences if kw_re.search(s) and not _mound_sentence_ok(s)]
            if ok_sents:
                verdicts[kw] = VERDICT_MOUND_ARCH
                accepted_kws.append(kw)
                mound_sentences_ok.extend(ok_sents)
            else:
                verdicts[kw] = VERDICT_MOUND_FP
                rejected_kws.append(kw)
                mound_sentences_bad.extend(bad_sents)

        elif kw in DOME_AMBIGUOUS:
            is_valid, _, note = dome_validate(all_text, kw)
            if is_valid:
                verdicts[kw] = VERDICT_DOME_VALID
                accepted_kws.append(kw)
            else:
                verdicts[kw] = VERDICT_DOME_REJECT
                rejected_kws.append(kw)

        else:
            # tholos — unambiguous dome-family
            verdicts[kw] = VERDICT_UNAMB
            accepted_kws.append(kw)

    return {
        "raw_kws":             raw_kws,
        "verdicts":            verdicts,
        "accepted_kws":        accepted_kws,
        "rejected_kws":        rejected_kws,
        "accepted":            bool(accepted_kws),
        "mound_sentences_ok":  mound_sentences_ok,
        "mound_sentences_bad": mound_sentences_bad,
    }

# ── Load corpus and audit ─────────────────────────────────────────────────────
corpus = load_corpus()

raw_sites       = []   # ALL sites matching any keyword
accepted_sites  = []   # sites that pass context validation
rejected_sites  = []   # sites that fail all keyword validations

for site_obj in corpus:
    if site_obj.category == "Natural":
        continue
    if not site_obj.has_coords:
        continue

    result = audit_site(site_obj)
    if result is None:
        continue

    dev   = _beru_dev(site_obj.longitude)
    tier  = tier_label(dev)

    entry = {
        "name":   site_obj.site,
        "lon":    site_obj.longitude,
        "dev":    dev,
        "dev_km": dev * BERU * 111.0,
        "tier":   tier,
        "ap":     is_aplus(tier),
        **result,
    }

    raw_sites.append(entry)

    if result["accepted"]:
        accepted_sites.append(entry)
    else:
        rejected_sites.append(entry)

# ── Counts for FP-rate calculation ───────────────────────────────────────────
# We want the FP rate for the "mound" keyword specifically: of all sites
# where "mound" is the *deciding* keyword (i.e. removing "mound" would
# exclude the site because no other raw keyword also passes validation),
# how many are false positives?
#
# Algorithm:
#   For each site that matched "mound" and has NO co-occurring unambiguous
#   keyword (tumulus/kofun/barrow/stupa/tholos), check whether any OTHER
#   raw keyword (dome/domed/domes/spherical) independently passes its own
#   context validation.  If not, this site is included/excluded solely
#   because of "mound" — that's the meaningful denominator.
#
# This is Group 3 in the debug analysis:
#   N = 22 (sites included or rejected solely because of mound)
#   Rejected = 10  →  45.5% FP rate

def _mound_is_deciding(entry) -> bool:
    """Return True if removing 'mound' would change this site's inclusion."""
    if "mound" not in entry["raw_kws"]:
        return False
    # If any unambiguous co-keyword exists, the site is accepted regardless
    if any(kw in entry["raw_kws"] for kw in MOUND_UNAMB + STUPA_KWS + ["tholos"]):
        return False
    # Check whether any ambiguous dome co-keyword independently passes
    # (uses the dome_filter result already stored in verdicts)
    for kw, verdict in entry["verdicts"].items():
        if kw == "mound":
            continue
        if kw in DOME_AMBIGUOUS and verdict == VERDICT_DOME_VALID:
            return False   # dome kw would include the site anyway
    return True

# Sites where mound is the deciding keyword
mound_deciding        = [e for e in raw_sites if _mound_is_deciding(e)]
mound_deciding_acc    = [e for e in mound_deciding if e["accepted"]]
mound_deciding_rej    = [e for e in mound_deciding if not e["accepted"]]

# Keep old names for compatibility with the FP-rate summary block
mound_only_raw      = mound_deciding
mound_only_rejected = mound_deciding_rej

N_raw       = len(raw_sites)
N_accepted  = len(accepted_sites)
N_rejected  = len(rejected_sites)

# Tier-A+ in raw vs validated populations
n_ap_raw      = sum(1 for e in raw_sites      if e["ap"])
n_ap_accepted = sum(1 for e in accepted_sites if e["ap"])
n_a_raw       = sum(1 for e in raw_sites      if is_a_or_better(e["tier"]))
n_a_accepted  = sum(1 for e in accepted_sites if is_a_or_better(e["tier"]))

bt_ap_raw  = binomtest(n_ap_raw,      N_raw,      P_NULL_AP, alternative="greater")
bt_ap_val  = binomtest(n_ap_accepted, N_accepted, P_NULL_AP, alternative="greater")
bt_a_raw   = binomtest(n_a_raw,       N_raw,      P_NULL_A,  alternative="greater")
bt_a_val   = binomtest(n_a_accepted,  N_accepted, P_NULL_A,  alternative="greater")

enr_ap_raw = (n_ap_raw  / N_raw)      / P_NULL_AP if N_raw      else 0
enr_ap_val = (n_ap_accepted / N_accepted) / P_NULL_AP if N_accepted else 0

fp_rate_overall = N_rejected / N_raw if N_raw else 0

# Per-keyword FP breakdown
kw_stats = {}
for kw in ALL_KEYWORDS:
    sites_with_kw = [e for e in raw_sites if kw in e["raw_kws"]]
    sites_acc     = [e for e in sites_with_kw if kw in e["accepted_kws"]]
    sites_rej     = [e for e in sites_with_kw if kw in e["rejected_kws"]]
    if not sites_with_kw:
        continue
    kw_stats[kw] = {
        "n_raw": len(sites_with_kw),
        "n_acc": len(sites_acc),
        "n_rej": len(sites_rej),
        "fp":    len(sites_rej) / len(sites_with_kw),
    }

# ── Print ─────────────────────────────────────────────────────────────────────
SEP  = "=" * 110
SEP2 = "─" * 110

print()
print(SEP)
print("  UNESCO WHC — MOUND/STUPA/BARROW/KOFUN KEYWORD CONTEXT AUDIT")
print(f"  Keywords: {ALL_KEYWORDS}")
print(f"  Anchor: Gerizim {GERIZIM}°E  |  BERU = {BERU}°")
print(f"  Ambiguous: {MOUND_AMB + DOME_AMBIGUOUS}  "
      f"(require sentence-level context validation)")
print(f"  Unambiguous: {MOUND_UNAMB + STUPA_KWS}  "
      f"(accepted without sentence check)")
print(SEP)

# ── Per-keyword FP summary ────────────────────────────────────────────────────
print()
print("  PER-KEYWORD CLASSIFICATION SUMMARY")
print()
print(f"  {'Keyword':<14}  {'N_raw':>6}  {'Accepted':>8}  {'Rejected':>8}  "
      f"{'FP rate':>8}  Type")
print(f"  {'─'*70}")
for kw in ALL_KEYWORDS:
    if kw not in kw_stats:
        continue
    ks = kw_stats[kw]
    kw_type = ("UNAMB" if kw in MOUND_UNAMB + STUPA_KWS
                else ("AMB-mound" if kw in MOUND_AMB
                      else ("AMB-dome" if kw in DOME_AMBIGUOUS
                            else "UNAMB-dome")))
    fp_str = f"{100*ks['fp']:>6.1f}%" if ks["n_rej"] > 0 else "   0.0%"
    print(f"  {kw:<14}  {ks['n_raw']:>6}  {ks['n_acc']:>8}  {ks['n_rej']:>8}  "
          f"{fp_str:>8}  {kw_type}")

# ── Per-site table ────────────────────────────────────────────────────────────
print()
print(SEP)
print("  ALL SITES — KEYWORD AUDIT TABLE")
print(SEP)
print()
print(SEP2)
print(f"  {'Site':<50}  {'Lon':>8}  {'Dev':>7}  {'km':>6}  T    Status")
print(SEP2)
for e in sorted(raw_sites, key=lambda x: x["dev"]):
    mark    = " ◀◀ A+" if e["ap"] else (" ◀ A" if e["tier"] == "A" else "")
    status  = "OK" if e["accepted"] else "FP?"
    kw_str  = ", ".join(e["raw_kws"])[:38]
    verdict = "; ".join(f"{k}:{v[:14]}" for k, v in e["verdicts"].items())
    print(f"  {e['name']:<50}  {e['lon']:>8.4f}  {e['dev']:>7.5f}  "
          f"{e['dev_km']:>6.1f}  {e['tier']}{mark:<7}  [{status}]  kw=[{kw_str}]")
    if VERBOSE:
        print(f"    Verdicts: {verdict}")
        for sent in e["mound_sentences_ok"]:
            print(f"      ✓ {sent[:120]}")
        for sent in e["mound_sentences_bad"]:
            print(f"      ✗ {sent[:120]}")
        print()

# ── Context-rejected sites ────────────────────────────────────────────────────
print()
print(SEP)
print(f"  CONTEXT-REJECTED SITES  (potential false positives)  N = {N_rejected}")
print(SEP)
print()
if rejected_sites:
    for e in sorted(rejected_sites, key=lambda x: x["dev"]):
        print(f"  [{e['tier']:>3}] {e['dev_km']:>5.1f} km  {e['name']}")
        for kw, verdict in e["verdicts"].items():
            print(f"         kw='{kw}'  →  {verdict}")
        if e["mound_sentences_bad"]:
            for sent in e["mound_sentences_bad"][:3]:
                print(f"           ✗  {sent[:110]}")
else:
    print("  None — all keyword matches passed context validation.")

# ── Statistical comparison ────────────────────────────────────────────────────
print()
print(SEP)
print("  FALSE-POSITIVE RATE SUMMARY")
print(SEP)
print()
print(f"  N_raw (any keyword, no validation) : {N_raw}")
print(f"  N_accepted (pass context validation): {N_accepted}")
print(f"  N_rejected (fail all validations)  : {N_rejected}")
print(f"  Overall FP rate (site-level)       : {100*fp_rate_overall:.1f}%  "
      f"({N_rejected}/{N_raw} sites rejected)")
print()
print(f"  Note: 'rejected' means the site matched ONLY ambiguous keywords")
print(f"  and all matching sentences failed context validation.  A site")
print(f"  with any unambiguous keyword (tumulus/kofun/stupa/…) is always accepted.")
print()

# FP rate for "mound" specifically (sites where mound is the deciding keyword)
if mound_only_raw:
    n_mo_raw = len(mound_only_raw)
    n_mo_rej = len(mound_only_rejected)
    print(f"  'mound' keyword — sites where it is the deciding inclusion factor")
    print(f"  (no unambiguous co-keyword; any dome/spherical co-kws also failed):")
    print(f"    N_deciding: {n_mo_raw}")
    print(f"    N_rejected: {n_mo_rej}")
    if n_mo_raw > 0:
        print(f"    FP rate   : {100*n_mo_rej/n_mo_raw:.1f}%  ({n_mo_rej}/{n_mo_raw} sites)")
    print()

print(SEP)
print("  STATISTICAL COMPARISON — RAW SWEEP vs CONTEXT-VALIDATED")
print(SEP)
print()
print(f"  {'Metric':<50}  {'Raw':>12}  {'Validated':>12}")
print(f"  {'─'*78}")
print(f"  {'N (population)':<50}  {N_raw:>12}  {N_accepted:>12}")
print(f"  {'n Tier-A+ hits':<50}  {n_ap_raw:>12}  {n_ap_accepted:>12}")
print(f"  {'Tier-A+ rate':<50}  {100*n_ap_raw/N_raw:>11.1f}%  "
      f"{100*n_ap_accepted/N_accepted:>11.1f}%")
print(f"  {'Enrichment vs 4% null (A+)':<50}  {enr_ap_raw:>11.2f}×  {enr_ap_val:>11.2f}×")
print(f"  {'Binomial p (A+)':<50}  {bt_ap_raw.pvalue:>12.4f}  {bt_ap_val.pvalue:>12.4f}")
print(f"  {'Significance (A+)':<50}  {sig(bt_ap_raw.pvalue):>12}  {sig(bt_ap_val.pvalue):>12}")
print(f"  {'n Tier-A hits':<50}  {n_a_raw:>12}  {n_a_accepted:>12}")
print(f"  {'Binomial p (A)':<50}  {bt_a_raw.pvalue:>12.4f}  {bt_a_val.pvalue:>12.4f}")
print()
if N_rejected > 0:
    print(f"  Interpretation:")
    print(f"    Context-filtering removes {N_rejected} site(s) ({100*fp_rate_overall:.1f}% of raw matches).")
    print(f"    If the signal survives or strengthens in the validated pop, context-")
    print(f"    filtering is doing real work and is not generating the result.")
    print(f"    If p_validated < p_raw: ✓ filtering removes noise (FPs dilute the signal).")
    print(f"    If p_validated > p_raw: ✗ filtering is suspiciously discarding aligned sites.")
else:
    print(f"  No sites were rejected — the raw and validated populations are identical.")
    print(f"  The false positive rate for these keywords is effectively 0% in this corpus.")
print()

# ── Accepted sites sorted by tier ────────────────────────────────────────────
print(SEP)
print(f"  ACCEPTED SITES  (pass context validation)  N = {N_accepted}")
print(SEP)
print()
print(SEP2)
print(f"  {'Site':<50}  {'Lon':>8}  {'km':>6}  T    Accepted keywords")
print(SEP2)
for e in sorted(accepted_sites, key=lambda x: x["dev"]):
    mark   = " ◀◀ A+" if e["ap"] else (" ◀ A" if e["tier"] == "A" else "")
    kw_str = ", ".join(e["accepted_kws"])
    print(f"  {e['name']:<50}  {e['lon']:>8.4f}  {e['dev_km']:>6.1f}  "
          f"{e['tier']}{mark:<7}  [{kw_str}]")

# ── LaTeX macros (always emitted — greppable by reproduce_all_macros.sh) ────
fp_rate_pct = round(100 * fp_rate_overall, 1)
mound_fp    = round(100 * len(mound_only_rejected) / len(mound_only_raw), 1) \
              if mound_only_raw else 0.0

if LATEX:
    print()
    print(SEP)
    print("  LATEX MACROS  (mound keyword context audit)")
    print(SEP)
    print()

print(f"  % Mound-evolution keyword audit results")
print(f"  \\newcommand{{\\NmoundRaw}}{{{N_raw}}}           % raw-sweep population (no context filter)")
print(f"  \\newcommand{{\\NmoundAccepted}}{{{N_accepted}}}       % context-validated population")
print(f"  \\newcommand{{\\NmoundRejected}}{{{N_rejected}}}        % sites failing all context validation")
print(f"  \\newcommand{{\\moundFPRate}}{{{fp_rate_pct}}}      % site-level FP rate (%) for mound-evo keywords")
print(f"  \\newcommand{{\\moundKwFPRate}}{{{mound_fp}}}     % FP rate (%) for 'mound' keyword specifically")
print(f"  \\newcommand{{\\NmoundDeciding}}{{{len(mound_only_raw)}}}       % sites where 'mound' is the deciding inclusion factor")
print(f"  \\newcommand{{\\NmoundDecidingRej}}{{{len(mound_only_rejected)}}}    % of those, rejected by context validation")
print(f"  \\newcommand{{\\NmoundRawAp}}{{{n_ap_raw}}}        % Tier-A+ in raw sweep")
print(f"  \\newcommand{{\\NmoundAccAp}}{{{n_ap_accepted}}}       % Tier-A+ in validated population")
print(f"  \\newcommand{{\\pmoundRawAp}}{{{bt_ap_raw.pvalue:.4f}}}   % p-value A+ (raw sweep)")
print(f"  \\newcommand{{\\pmoundAccAp}}{{{bt_ap_val.pvalue:.4f}}}   % p-value A+ (validated)")
print(f"  \\newcommand{{\\moundEnrichRaw}}{{{enr_ap_raw:.2f}}}  % enrichment A+ (raw sweep)")
print(f"  \\newcommand{{\\moundEnrichAcc}}{{{enr_ap_val:.2f}}}  % enrichment A+ (validated)")

print()
print(f"  Audit complete.  Run with --verbose to see matching sentences per site.")
print(f"                   Run with --latex for compact macro-only output.")
print()
