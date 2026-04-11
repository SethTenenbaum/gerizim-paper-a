"""
tumulus_dome_evolution_test.py
================================
Context-validated version of the hemispherical mound-evolution corpus test.

RELATIONSHIP TO RAW SWEEP
--------------------------
  tumulus_dome_evolution_raw_sweep.py  — RAW, no context validation.  This is
      the primary confirmatory version (Test 2b / Sensitivity), maximally
      resistant to accusations of post-hoc filtering.

  THIS FILE (tumulus_dome_evolution_test.py)  — CONTEXT-VALIDATED, Exploratory.
      Applies the same sentence-level disambiguation used by
      spherical_monument_test.py:
        • Mound keywords (tumulus/tumuli/barrow/kofun/mound) are treated as
          UNAMBIGUOUS in an archaeological context — accepted without validation.
        • Stupa/dagoba/chorten are unambiguous — accepted without validation.
        • Dome/domed/domes/spherical are AMBIGUOUS — accepted only when the
          matching sentence contains a positive architectural-context term and
          no negative-context term, via lib/dome_filter.validate_keyword_match.

      If context-filtering removes a site from the validated corpus that was
      present in the raw sweep, it means the keyword matched in a non-architectural
      sentence (e.g. a natural mound, a geological dome).  The validated corpus
      should therefore show a STRONGER signal than the raw sweep, not a weaker one.

EXPLORATORY STATUS
------------------
This version is labelled Exploratory because the specific mound keyword set was
refined after the initial corpus inspection (post-hoc). The RAW sweep (which
includes everything that matches any mound/stupa/dome keyword) is the
pre-specified confirmatory version.

USAGE
-----
    cd /path/to/gerizim-paper-a
    python3 analysis/unesco/tumulus_dome_evolution_test.py
    python3 analysis/unesco/tumulus_dome_evolution_test.py --audit

OUTPUT
------
    Console table with per-stage statistics and binomial tests for the
    context-validated population, plus LaTeX \\newcommand macros and
    ResultsStore write for pEvoAp_validated.
"""

import re
import sys
import numpy as np
from pathlib import Path
from scipy.stats import binomtest, fisher_exact, chisquare

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_A_MAX,
    P_NULL_AP, P_NULL_A,
    deviation as _beru_dev, tier_label, is_aplus, is_a_or_better,
)
from lib.stats import significance_label as sig
from lib.results_store import ResultsStore
from lib.dome_filter import (
    UNAMBIGUOUS_KEYWORDS as DOME_UNAMBIGUOUS,
    AMBIGUOUS_KEYWORDS   as DOME_AMBIGUOUS,
    validate_keyword_match,
)

# ── Load keyword sets from keywords.json ─────────────────────────────────────
import json
_KW_PATH = Path(__file__).parent.parent.parent / "keywords.json"
with open(_KW_PATH) as _f:
    _KW = json.load(_f)

_evo = _KW["mound_evolution"]
# Mound family — unambiguous in archaeological context (accepted without
# sentence-level context validation, consistent with stupa/tholos treatment).
MOUND_KEYWORDS_UNAMB = _evo["mound_unambiguous"]

# "mound" is slightly ambiguous (natural mounds exist) so we require the
# sentence to contain an archaeological context term.
MOUND_KEYWORDS_AMB   = _evo["mound_ambiguous"]

MOUND_KEYWORDS_ALL   = MOUND_KEYWORDS_UNAMB + MOUND_KEYWORDS_AMB

# Archaeological-context patterns for the word "mound" (must appear in the
# same sentence for the match to count).
MOUND_POSITIVE_CONTEXT_RES = [
    re.compile(r"\b" + p + r"\b", re.IGNORECASE)
    for p in _evo["mound_positive_context"]
]

def _validate_mound_sentence(sentence: str) -> bool:
    """Accept 'mound' if the sentence contains an archaeological context term."""
    return any(pat.search(sentence) for pat in MOUND_POSITIVE_CONTEXT_RES)

# Stupa family — already unambiguous via dome_filter
STUPA_KEYWORDS = _evo["stupa"]

# Dome family — dome_filter handles these with full POSITIVE/NEGATIVE context
DOME_KEYWORDS  = _evo["dome"]

ALL_KEYWORDS = MOUND_KEYWORDS_ALL + STUPA_KEYWORDS + DOME_KEYWORDS

# Compile word-boundary regexes for all keywords
KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    for kw in ALL_KEYWORDS
}

# ── Stage classification ──────────────────────────────────────────────────────
def classify_stage(matched_kws: list) -> set:
    stages = set()
    for kw in matched_kws:
        if kw in MOUND_KEYWORDS_ALL:
            stages.add("mound")
        elif kw in STUPA_KEYWORDS:
            stages.add("stupa")
        elif kw in DOME_KEYWORDS:
            stages.add("dome")
    return stages

# ── Date extraction (identical to raw sweep) ──────────────────────────────────
BCE_RE            = re.compile(r"(\d{3,4})\s*(?:BCE|BC|B\.C\.E?\.?)",  re.IGNORECASE)
CE_RE             = re.compile(r"(\d{1,4})\s*(?:CE|AD|A\.D\.)",          re.IGNORECASE)
CENTURY_BCE_RE    = re.compile(r"(\d{1,2})(?:st|nd|rd|th)\s*(?:century|c\.?)\s*(?:BCE|BC)", re.IGNORECASE)
CENTURY_CE_RE     = re.compile(r"(\d{1,2})(?:st|nd|rd|th)\s*(?:century|c\.?)\s*(?:CE|AD)",  re.IGNORECASE)
MILLENNIUM_BCE_RE = re.compile(r"(\d{1,2})(?:st|nd|rd|th)\s*millennium\s*(?:BCE|BC)",        re.IGNORECASE)

def _ord(n: str) -> str:
    return {"1": "st", "2": "nd", "3": "rd"}.get(n, "th")

def extract_earliest_date(text: str):
    dates = []
    for m in BCE_RE.finditer(text):
        dates.append((-int(m.group(1)), f"{m.group(1)} BCE"))
    for m in MILLENNIUM_BCE_RE.finditer(text):
        n = m.group(1)
        dates.append((-int(n) * 1000, f"{n}{_ord(n)} millennium BCE"))
    for m in CENTURY_BCE_RE.finditer(text):
        n = m.group(1)
        dates.append((-int(n) * 100, f"{n}{_ord(n)} c. BCE"))
    for m in CE_RE.finditer(text):
        yr = int(m.group(1))
        if yr < 50:
            continue
        dates.append((yr, f"{yr} CE"))
    for m in CENTURY_CE_RE.finditer(text):
        n = m.group(1)
        dates.append((int(n) * 100, f"{n}{_ord(n)} c. CE"))
    if dates:
        dates.sort(key=lambda x: x[0])
        return dates[0]
    return None, None

# ── Context-validated corpus loader ───────────────────────────────────────────
def load_validated_corpus():
    """
    Returns (selected, raw_rejected, total_xml).

    selected      — list of site dicts that PASS context validation
    raw_rejected  — list of dicts for sites whose keywords were all
                    rejected by context rules (with notes)
    total_xml     — total corpus size before any filtering
    """
    corpus = load_corpus()
    total_xml    = len(corpus)
    selected     = []
    raw_rejected = []

    for site_obj in corpus:
        if site_obj.category == "Natural":
            continue
        if not site_obj.has_coords:
            continue

        full = site_obj.full_text   # lowercased

        # Stage 1: raw keyword match
        raw_matched = [kw for kw in ALL_KEYWORDS if KEYWORD_RES[kw].search(full)]
        if not raw_matched:
            continue

        # Stage 2: context validation per keyword
        all_text = (site_obj.short_description or "")
        if site_obj.extended_description:
            all_text += " " + site_obj.extended_description

        sentences = re.split(r"(?<=[.!?])\s+", all_text)

        validated_kws    = []
        validation_notes = {}

        for kw in raw_matched:
            kw_re = KEYWORD_RES[kw]

            if kw in MOUND_KEYWORDS_UNAMB:
                # Unambiguous archaeological terms — accept immediately
                validated_kws.append(kw)
                validation_notes[kw] = "unambiguous"

            elif kw in MOUND_KEYWORDS_AMB:
                # "mound" — require archaeological context sentence
                valid_sents = [
                    s.strip() for s in sentences
                    if kw_re.search(s) and _validate_mound_sentence(s)
                ]
                if valid_sents:
                    validated_kws.append(kw)
                    validation_notes[kw] = "mound-context-validated"
                else:
                    validation_notes[kw] = "mound-rejected-no-arch-context"

            elif kw in STUPA_KEYWORDS or kw in DOME_UNAMBIGUOUS:
                # Unambiguous dome-family terms (stupa, dagoba, chorten, tholos)
                validated_kws.append(kw)
                validation_notes[kw] = "unambiguous"

            else:
                # Ambiguous dome-family terms (dome, domed, domes, spherical)
                is_valid, _, note = validate_keyword_match(all_text, kw)
                validation_notes[kw] = note
                if is_valid:
                    validated_kws.append(kw)

        if not validated_kws:
            raw_rejected.append({
                "name":   site_obj.site,
                "lon":    site_obj.longitude,
                "raw":    raw_matched,
                "notes":  validation_notes,
            })
            continue

        stages    = classify_stage(validated_kws)
        dev       = _beru_dev(site_obj.longitude)
        tier      = tier_label(dev)
        date_yr, date_label = extract_earliest_date(all_text)

        selected.append({
            "name":       site_obj.site,
            "lon":        site_obj.longitude,
            "dev":        dev,
            "dev_km":     dev * BERU * 111.0,
            "tier":       tier,
            "ap":         is_aplus(tier),
            "stages":     stages,
            "keywords":   validated_kws,
            "raw_kws":    raw_matched,
            "notes":      validation_notes,
            "date_yr":    date_yr,
            "date_label": date_label,
        })

    return selected, raw_rejected, total_xml

# ── --audit mode ──────────────────────────────────────────────────────────────
if "--audit" in sys.argv:
    sel, rejected, total = load_validated_corpus()
    print(f"\nAUDIT — tumulus_dome_evolution_test.py (context-validated)")
    print(f"Total corpus : {total}  |  Included: {len(sel)}  |  Rejected: {len(rejected)}")
    print()
    print("INCLUDED:")
    print("─" * 100)
    for s in sorted(sel, key=lambda x: x["name"]):
        print(f"  {s['name']}")
        print(f"    VALIDATED kws : {s['keywords']}  stages: {sorted(s['stages'])}")
        print(f"    notes: {s['notes']}")
    print()
    print("CONTEXT-REJECTED:")
    print("─" * 100)
    for r in rejected:
        print(f"  {r['name']}")
        print(f"    RAW: {r['raw']}   notes: {r['notes']}")
    sys.exit(0)

# ── Main analysis ─────────────────────────────────────────────────────────────
selected, raw_rejected, total_xml = load_validated_corpus()

dome_only  = [e for e in selected if "dome" in e["stages"] or "stupa" in e["stages"]
              and "mound" not in e["stages"]]
mound_only = [e for e in selected if "mound" in e["stages"]
              and "dome" not in e["stages"] and "stupa" not in e["stages"]]
overlap    = [e for e in selected if "mound" in e["stages"]
              and ("dome" in e["stages"] or "stupa" in e["stages"])]

# ── Stats ─────────────────────────────────────────────────────────────────────
N    = len(selected)
n_ap = sum(1 for e in selected if e["ap"])
n_a  = sum(1 for e in selected if is_a_or_better(e["tier"]))
n_b  = sum(1 for e in selected if e["tier"] == "B")
n_c  = sum(1 for e in selected if e["tier"] == "C")

bt_ap = binomtest(n_ap, N, P_NULL_AP, alternative="greater")
bt_a  = binomtest(n_a,  N, P_NULL_A,  alternative="greater")
enr_ap = (n_ap / N) / P_NULL_AP if N else 0
enr_a  = (n_a  / N) / P_NULL_A  if N else 0

obs_bins = [0] * 5
for e in selected:
    obs_bins[min(int(e["dev"] / 0.010), 4)] += 1
_, chi_p = chisquare(obs_bins, f_exp=[N / 5.0] * 5)

# Anchor sweep
sweep  = np.arange(34.0, 37.001, 0.001)
cnt_Ap = np.array([
    sum(1 for e in selected
        if abs(abs(e["lon"] - a) / BERU - round(abs(e["lon"] - a) / BERU * 10) / 10) <= TIER_APLUS)
    for a in sweep
])
pctile = float(np.mean(cnt_Ap <= n_ap)) * 100

# Stage sub-stats
def stage_stats(entries):
    ns  = len(entries)
    nap = sum(1 for e in entries if e["ap"])
    na  = sum(1 for e in entries if is_a_or_better(e["tier"]))
    p   = binomtest(nap, ns, P_NULL_AP, alternative="greater").pvalue if ns > 0 else 1.0
    enr = (nap / ns) / P_NULL_AP if ns > 0 else 0
    return ns, nap, na, p, enr

n_mound_stage    = sum(1 for e in selected if "mound" in e["stages"])
n_stupa_stage    = sum(1 for e in selected if "stupa" in e["stages"])
n_dome_stage     = sum(1 for e in selected if "dome"  in e["stages"])
n_mound_stage_ap = sum(1 for e in selected if "mound" in e["stages"] and e["ap"])
n_stupa_stage_ap = sum(1 for e in selected if "stupa" in e["stages"] and e["ap"])
n_dome_stage_ap  = sum(1 for e in selected if "dome"  in e["stages"] and e["ap"])
n_mo_ap          = sum(1 for e in mound_only if e["ap"])

p_mound = binomtest(n_mound_stage_ap, n_mound_stage, P_NULL_AP, "greater").pvalue if n_mound_stage else 1.0
p_stupa = binomtest(n_stupa_stage_ap, n_stupa_stage, P_NULL_AP, "greater").pvalue if n_stupa_stage else 1.0
p_dome  = binomtest(n_dome_stage_ap,  n_dome_stage,  P_NULL_AP, "greater").pvalue if n_dome_stage  else 1.0

# ── Print ─────────────────────────────────────────────────────────────────────
SEP  = "=" * 100
SEP2 = "─" * 100

print()
print(SEP)
print("  UNESCO WHC — HEMISPHERICAL MOUND EVOLUTION  (CONTEXT-VALIDATED, Exploratory)")
print("  Exploratory robustness check of Test 2b (raw sweep primary)")
print(f"  Keywords: mound={MOUND_KEYWORDS_ALL}  stupa={STUPA_KEYWORDS}  dome={DOME_KEYWORDS}")
print(f"  Anchor: Gerizim {GERIZIM}°E  |  BERU = {BERU}°  |  "
      f"Tier-A+ ≤ {TIER_APLUS} beru (≤ {TIER_APLUS*BERU*111:.0f} km)")
print(f"  Context-rejected: {len(raw_rejected)} sites  |  Included: {N}")
print(f"  H₀: longitudes uniform w.r.t. 0.1-beru harmonics  |  H₁: rate > geometric null")
print(SEP)

# Per-site table
print()
print(SEP2)
print(f"  {'Site':<52}  {'Lon':>8}  {'Dev':>7}  {'km':>6}  T    Stages       kw (validated)")
print(SEP2)
for e in sorted(selected, key=lambda x: x["dev"]):
    mark   = " ◀◀ A+" if e["ap"] else (" ◀ A" if e["tier"] == "A" else "")
    stages = "+".join(sorted(e["stages"]))
    kw_str = ", ".join(e["keywords"])[:38]
    print(f"  {e['name']:<52}  {e['lon']:>8.4f}  {e['dev']:>7.5f}  "
          f"{e['dev_km']:>6.1f}  {e['tier']}{mark:<7}  {stages:<12}  [{kw_str}]")

# Combined stats
print()
print(SEP)
print("  STATISTICAL RESULTS — CONTEXT-VALIDATED HEMISPHERICAL POPULATION")
print(SEP)
print()
print(f"  {'Metric':<45}  {'Obs':>5}  {'Exp(H₀)':>9}  {'Enrich':>7}  {'p':>8}  Sig")
print(f"  {'-'*85}")
print(f"  {'Tier-A+  (≤0.002 beru, ≤6.7 km)  ◀ PRIMARY':<45}  {n_ap:>5}  "
      f"{N*P_NULL_AP:>9.2f}  {enr_ap:>6.2f}×  {bt_ap.pvalue:>8.4f}  {sig(bt_ap.pvalue)}")
print(f"  {'Tier-A   (≤0.010 beru, ≤33 km)':<45}  {n_a:>5}  "
      f"{N*P_NULL_A:>9.2f}  {enr_a:>6.2f}×  {bt_a.pvalue:>8.4f}  {sig(bt_a.pvalue)}")
print(f"  {'Tier-B   (≤0.050 beru)':<45}  {n_b:>5}")
print(f"  {'Tier-C   (>0.050 beru)':<45}  {n_c:>5}")
print(f"  {'χ²-uniform (5 bins, df=4)':<45}  {'':>5}  {'':>9}  {'':>7}  {chi_p:>8.4f}  {sig(chi_p)}")
print(f"  {'Anchor sweep percentile (34–37°E)':<45}  {n_ap:>5}  {'':>9}  {'':>7}  "
      f"{'':>8}  {pctile:.0f}th pctile")

# By stage
print()
print(SEP)
print("  BY EVOLUTIONARY STAGE  (context-validated)")
print(SEP)
print()
print(f"  {'Stage':<35}  {'N':>5}  {'A+':>4}  {'A':>4}  {'Enrich':>7}  {'p (A+)':>8}  Sig")
print(f"  {'-'*75}")
stage_stat_results = {}   # keyed by short name: "mound", "stupa", "dome"
for stage_key, stage_name, stage_kws in [
    ("mound", "Mound (tumulus/barrow/kofun/mound)", MOUND_KEYWORDS_ALL),
    ("stupa", "Stupa (stupa/dagoba/chorten)",       STUPA_KEYWORDS),
    ("dome",  "Dome  (tholos/dome/domed/spherical)", DOME_KEYWORDS),
]:
    entries = [e for e in selected if any(kw in e["keywords"] for kw in stage_kws)]
    ns, nap_s, na_s, p_s, enr_s = stage_stats(entries)
    stage_stat_results[stage_key] = {"n": ns, "nap": nap_s, "rate": nap_s/ns if ns else 0, "p": p_s, "enr": enr_s}
    print(f"  {stage_name:<35}  {ns:>5}  {nap_s:>4}  {na_s:>4}  {enr_s:>6.2f}×  {p_s:>8.4f}  {sig(p_s)}")

# Context-rejected sites
if raw_rejected:
    print()
    print(SEP)
    print(f"  CONTEXT-REJECTED SITES  ({len(raw_rejected)} — keyword matched but failed validation)")
    print(SEP)
    for r in raw_rejected:
        print(f"    ✗  {r['name'][:70]}")
        print(f"       raw kw: {r['raw']}   notes: {r['notes']}")

# A+ hits
print()
print(SEP)
print(f"  ALL TIER-A+ HITS  ({n_ap})")
print(SEP)
print()
for e in sorted([e for e in selected if e["ap"]], key=lambda x: x["dev"]):
    stages = "+".join(sorted(e["stages"]))
    kws    = ", ".join(e["keywords"])
    date_s = e["date_label"] or "?"
    print(f"  {e['dev_km']:>5.1f} km  {e['name']:<52}  [{stages:<11}]  kw={kws}  ~{date_s}")

# Chronological breakdown
print()
print(SEP)
print("  CHRONOLOGICAL BREAKDOWN")
print(SEP)
dated = [(e, e["date_yr"]) for e in selected if e["date_yr"] is not None]
eras  = [
    ("Deep antiquity (pre-1000 BCE)",   lambda y: y < -1000),
    ("Classical (1000 BCE – 1 CE)",     lambda y: -1000 <= y < 1),
    ("Late antiquity (1–500 CE)",       lambda y: 1 <= y < 500),
    ("Medieval (500–1500 CE)",          lambda y: 500 <= y < 1500),
    ("Early modern (1500–1800 CE)",     lambda y: 1500 <= y < 1800),
    ("Modern (post-1800 CE)",           lambda y: y >= 1800),
]
print()
print(f"  {'Era':<35}  {'N':>4}  {'A+':>4}  {'Rate':>6}  Stages")
print(f"  {'─'*70}")
for era_name, era_fn in eras:
    ee     = [e for e, y in dated if era_fn(y)]
    n_e    = len(ee)
    n_e_ap = sum(1 for e in ee if e["ap"])
    stg    = ", ".join(sorted({s for e in ee for s in e["stages"]}))
    rate   = f"{100*n_e_ap/n_e:.1f}%" if n_e > 0 else "—"
    print(f"  {era_name:<35}  {n_e:>4}  {n_e_ap:>4}  {rate:>6}  {stg}")

# ── LaTeX macros ──────────────────────────────────────────────────────────────
print()
print(SEP)
print("  LATEX MACROS  (context-validated — Exploratory)")
print(SEP)
print()
print(f"  \\newcommand{{\\NevoValidTotal}}{{{N}}}              % total dome-evolution corpus (context-validated)")
print(f"  \\newcommand{{\\NevoValidAp}}{{{n_ap}}}               % A+ sites (context-validated)")
print(f"  \\newcommand{{\\evoValidApRate}}{{{100*n_ap/N:.1f}}}             % A+ rate, context-validated (%)")
print(f"  \\newcommand{{\\evoValidEnrichAp}}{{{enr_ap:.2f}}}            % enrichment ratio A+ (context-validated)")
print(f"  \\newcommand{{\\pEvoApValidated}}{{{bt_ap.pvalue:.4f}}}          % p-value A+ binomial, context-validated (Exploratory)")
print(f"  \\newcommand{{\\pEvoAValidated}}{{{bt_a.pvalue:.4f}}}           % p-value A  binomial, context-validated (Exploratory)")
print(f"  \\newcommand{{\\NevoValidRejected}}{{{len(raw_rejected)}}}           % sites rejected by context validation")

# NOTE: Stage-level evo* macros (evoDomeApRate, pEvoDome, etc.) are emitted by
# tumulus_dome_evolution_raw_sweep.py (authoritative). Do not emit them here.
print()

# ── Write to results store ────────────────────────────────────────────────────
ResultsStore().write_many({
    "pEvoAp_validated": bt_ap.pvalue,     # binomial p, A+ (evolution, context-validated) — Exploratory 2bx
    "pEvoA_validated":  bt_a.pvalue,      # binomial p, A  (evolution, context-validated)
    "NevoValidTotal":   N,                # corpus size, context-validated
    "NevoValidAp":      n_ap,             # A+ hits, context-validated
    "NevoValidRejected": len(raw_rejected), # sites removed by context filter
    **{f"evo{sfx}ApRate": round(100 * stage_stat_results[k]["rate"], 1)
       for k, sfx in [("dome","Dome"),("mound","Mound"),("stupa","Stupa")]
       if k in stage_stat_results},
    **{f"evo{sfx}N": stage_stat_results[k]["n"]
       for k, sfx in [("dome","Dome"),("mound","Mound"),("stupa","Stupa")]
       if k in stage_stat_results},
    **{f"pEvo{sfx}": stage_stat_results[k]["p"]
       for k, sfx in [("dome","Dome"),("mound","Mound"),("stupa","Stupa")]
       if k in stage_stat_results},
})
print(f"  ✓ Results written to data/store/results.json")
print(f"    pEvoAp_validated = {bt_ap.pvalue:.6f}  {sig(bt_ap.pvalue)}")
