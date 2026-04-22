"""
tumulus_dome_evolution_raw_sweep.py
=====================================
Raw positive keyword sweep for the hemispherical mound-evolution corpus.

Like tumulus_dome_evolution_test.py (Test 2b) but performs NO context
validation and NO disambiguation.  Any Cultural/Mixed site whose full text
(name + short_description + extended OUV text) contains one or more of the
mound/stupa/dome FORM_KEYWORDS is included outright.

PURPOSE
-------
Provides a version of the evolution test that is maximally resistant to
accusations of exploratory bias from the context-filtering step.

  - Raw sweep is the primary confirmatory version (Test 2b).
  - Context-validated version (tumulus_dome_evolution_test.py) is retained as
    an Exploratory cross-check: it should yield a similar or stronger signal
    because context-filtering removes false-positive keyword hits.

USAGE
-----
    cd /path/to/gerizim-analysis
    python3 analysis/unesco/tumulus_dome_evolution_raw_sweep.py
    python3 analysis/unesco/tumulus_dome_evolution_raw_sweep.py --audit

OUTPUT
------
    Console table with per-stage statistics and binomial tests for the raw
    (over-inclusive) population.
"""

import re
import sys
import numpy as np
from pathlib import Path
from scipy.stats import binomtest, fisher_exact, chisquare

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus
from lib.beru import (
    GERIZIM, BERU, TIER_APP, TIER_APLUS, TIER_A_MAX, HARMONIC_STEP,
    P_NULL_APP, P_NULL_AP, P_NULL_A,
    TIER_APLUS_LABEL, TIER_A_LABEL, TIER_B_LABEL, TIER_C_LABEL,
    deviation as _beru_dev, tier_label, is_aplusplus, is_aplus, is_a_or_better,
)
from lib.stats import significance_label as sig
from lib.results_store import ResultsStore

# ── Load keyword sets from keywords.json ─────────────────────────────────────
import json
_KW_PATH = Path(__file__).parent.parent.parent / "keywords.json"
with open(_KW_PATH) as _f:
    _KW = json.load(_f)

_evo = _KW["mound_evolution"]
# Mound family — always unambiguous in the mound-evolution context
MOUND_KEYWORDS = _evo["mound_unambiguous"] + _evo["mound_ambiguous"]

# Stupa family — unambiguous monumental spherical architecture
STUPA_KEYWORDS = _evo["stupa"]

# Dome family — includes words that are ambiguous in isolation; we skip
# context-checking here intentionally (raw sweep).
DOME_KEYWORDS = _evo["dome"]

ALL_KEYWORDS = MOUND_KEYWORDS + STUPA_KEYWORDS + DOME_KEYWORDS

# Compile word-boundary regexes
KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    for kw in ALL_KEYWORDS
}

# ── Stage classification ──────────────────────────────────────────────────────
def classify_stage(matched_kws: list) -> set:
    stages = set()
    for kw in matched_kws:
        if kw in MOUND_KEYWORDS:
            stages.add("mound")
        elif kw in STUPA_KEYWORDS:
            stages.add("stupa")
        elif kw in DOME_KEYWORDS:
            stages.add("dome")
    return stages

# ── Date extraction (same regexes as validated version) ──────────────────────
BCE_RE           = re.compile(r"(\d{3,4})\s*(?:BCE|BC|B\.C\.E?\.?)", re.IGNORECASE)
CE_RE            = re.compile(r"(\d{1,4})\s*(?:CE|AD|A\.D\.)",        re.IGNORECASE)
CENTURY_BCE_RE   = re.compile(r"(\d{1,2})(?:st|nd|rd|th)\s*(?:century|c\.?)\s*(?:BCE|BC)", re.IGNORECASE)
CENTURY_CE_RE    = re.compile(r"(\d{1,2})(?:st|nd|rd|th)\s*(?:century|c\.?)\s*(?:CE|AD)",  re.IGNORECASE)
MILLENNIUM_BCE_RE = re.compile(r"(\d{1,2})(?:st|nd|rd|th)\s*millennium\s*(?:BCE|BC)", re.IGNORECASE)

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

# ── Load corpus — raw sweep, no validation ────────────────────────────────────
corpus = load_corpus()

# ── Full-corpus background counts (for Fisher exact) ─────────────────────────
_all_cultural = [s for s in corpus if s.category != "Natural" and s.has_coords]
N_CORPUS     = len(_all_cultural)
K_APP_CORPUS = sum(1 for s in _all_cultural if is_aplusplus(tier_label(_beru_dev(s.longitude))))
K_AP_CORPUS  = sum(1 for s in _all_cultural if is_aplus(tier_label(_beru_dev(s.longitude))))
K_A_CORPUS   = sum(1 for s in _all_cultural if is_a_or_better(tier_label(_beru_dev(s.longitude))))

selected   = []   # all raw matches
dome_only  = []   # dome/stupa stage only (no mound kw)
mound_only = []   # mound stage only (no dome/stupa kw)
overlap    = []   # both mound and dome/stupa kw

for site_obj in corpus:
    if site_obj.category == "Natural":
        continue
    if not site_obj.has_coords:
        continue

    full = site_obj.full_text  # already lowercased

    matched_kws = [kw for kw in ALL_KEYWORDS if KEYWORD_RES[kw].search(full)]
    if not matched_kws:
        continue

    stages = classify_stage(matched_kws)
    dev    = _beru_dev(site_obj.longitude)
    tier   = tier_label(dev)
    date_yr, date_label = extract_earliest_date(
        (site_obj.short_description or "") + " " + (site_obj.extended_description or "")
    )

    entry = {
        "name":       site_obj.site,
        "lon":        site_obj.longitude,
        "dev":        dev,
        "dev_km":     dev * BERU * 111.0,
        "tier":       tier,
        "ap":         is_aplus(tier),
        "stages":     stages,
        "keywords":   matched_kws,
        "date_yr":    date_yr,
        "date_label": date_label,
    }

    selected.append(entry)
    has_mound       = "mound" in stages
    has_dome_stupa  = "dome" in stages or "stupa" in stages

    if has_mound and has_dome_stupa:
        overlap.append(entry)
    elif has_mound:
        mound_only.append(entry)
    else:
        dome_only.append(entry)

# ── Stats ─────────────────────────────────────────────────────────────────────
N     = len(selected)
n_app = sum(1 for e in selected if is_aplusplus(e["tier"]))
n_ap  = sum(1 for e in selected if e["ap"])
n_a   = sum(1 for e in selected if is_a_or_better(e["tier"]))
n_b   = sum(1 for e in selected if e["tier"] == "B")
n_c   = sum(1 for e in selected if e["tier"] == "C")

bt_app = binomtest(n_app, N, P_NULL_APP, alternative="greater")
bt_ap = binomtest(n_ap, N, P_NULL_AP, alternative="greater")
bt_a  = binomtest(n_a,  N, P_NULL_A,  alternative="greater")

enr_app = (n_app / N) / P_NULL_APP if N else 0
enr_ap = (n_ap / N) / P_NULL_AP if N else 0
enr_a  = (n_a  / N) / P_NULL_A  if N else 0

obs_bins = [0] * 5
for e in selected:
    obs_bins[min(int(e["dev"] / TIER_A_MAX), 4)] += 1
_, chi_p = chisquare(obs_bins, f_exp=[N / 5.0] * 5)

# Anchor sweep
sweep  = np.arange(34.0, 37.001, 0.001)
cnt_Ap = np.array([
    sum(1 for e in selected
        if abs(abs(e["lon"] - a) / BERU - round(abs(e["lon"] - a) / BERU / HARMONIC_STEP) * HARMONIC_STEP) <= TIER_APLUS)
    for a in sweep
])
pctile = float(np.mean(cnt_Ap <= n_ap)) * 100

# Stage sub-stats
def stage_stats(entries, stage_name):
    ns  = len(entries)
    nap = sum(1 for e in entries if e["ap"])
    na  = sum(1 for e in entries if is_a_or_better(e["tier"]))
    p   = binomtest(nap, ns, P_NULL_AP, alternative="greater").pvalue if ns > 0 else 1.0
    enr = (nap / ns) / P_NULL_AP if ns > 0 else 0
    return ns, nap, na, p, enr

# ── Print ─────────────────────────────────────────────────────────────────────
SEP  = "=" * 100
SEP2 = "─" * 100

print()
print(SEP)
print("  UNESCO WHC — HEMISPHERICAL MOUND EVOLUTION  (RAW KEYWORD SWEEP, no context validation)")
print("  Test 2b primary confirmatory version")
print(f"  Keywords: mound={MOUND_KEYWORDS}  stupa={STUPA_KEYWORDS}  dome={DOME_KEYWORDS}")
print(f"  Anchor: Gerizim {GERIZIM}°E  |  BERU = {BERU}°  |  "
      f"Tier-A+ ≤ {TIER_APLUS} beru (≤ {TIER_APLUS*BERU*111:.0f} km)")
print(f"  H₀: longitudes uniform w.r.t. 0.1-beru harmonics  |  H₁: rate > geometric null")
print(f"  Population: N={N}  ({len(dome_only)} dome/stupa-only, "
      f"{len(mound_only)} mound-only, {len(overlap)} overlap)")
print(SEP)

# ── Per-site table ────────────────────────────────────────────────────────────
print()
print(SEP2)
print(f"  {'Site':<52}  {'Lon':>8}  {'Dev':>7}  {'km':>6}  T    Stages       Keywords")
print(SEP2)
for e in sorted(selected, key=lambda x: x["dev"]):
    mark   = " ◀◀ A+" if e["ap"] else (" ◀ A" if e["tier"] == "A" else "")
    stages = "+".join(sorted(e["stages"]))
    kw_str = ", ".join(e["keywords"])[:38]
    print(f"  {e['name']:<52}  {e['lon']:>8.4f}  {e['dev']:>7.5f}  "
          f"{e['dev_km']:>6.1f}  {e['tier']}{mark:<7}  {stages:<12}  [{kw_str}]")

# ── Combined stats ────────────────────────────────────────────────────────────
print()
print(SEP)
print("  STATISTICAL RESULTS — COMBINED HEMISPHERICAL POPULATION (RAW)")
print(SEP)
print()
print(f"  {'Metric':<45}  {'Obs':>5}  {'Exp(H₀)':>9}  {'Enrich':>7}  {'p':>8}  Sig")
print(f"  {'-'*85}")
print(f"  {f'Tier-A+  ({TIER_APLUS_LABEL})  ◀ PRIMARY':<45}  {n_ap:>5}  "
      f"{N*P_NULL_AP:>9.2f}  {enr_ap:>6.2f}×  {bt_ap.pvalue:>8.4f}  {sig(bt_ap.pvalue)}")
print(f"  {f'Tier-A   ({TIER_A_LABEL})':<45}  {n_a:>5}  "
      f"{N*P_NULL_A:>9.2f}  {enr_a:>6.2f}×  {bt_a.pvalue:>8.4f}  {sig(bt_a.pvalue)}")
print(f"  {f'Tier-B   ({TIER_B_LABEL})':<45}  {n_b:>5}")
print(f"  {f'Tier-C   ({TIER_C_LABEL})':<45}  {n_c:>5}")
print(f"  {'χ²-uniform (5 bins, df=4)':<45}  {'':>5}  {'':>9}  {'':>7}  {chi_p:>8.4f}  {sig(chi_p)}")
print(f"  {'Anchor sweep percentile (34–37°E)':<45}  {n_ap:>5}  {'':>9}  {'':>7}  "
      f"{'':>8}  {pctile:.0f}th pctile")

# ── By evolutionary stage ─────────────────────────────────────────────────────
print()
print(SEP)
print("  BY EVOLUTIONARY STAGE")
print(SEP)
print()
print(f"  {'Stage':<35}  {'N':>5}  {'A+':>4}  {'A':>4}  {'Enrich':>7}  {'p (A+)':>8}  Sig")
print(f"  {'-'*75}")

for stage_name, stage_kws in [
    ("Mound (tumulus/tumuli/barrow/kofun)", MOUND_KEYWORDS),
    ("Stupa (stupa/stupas/dagoba/chorten)", STUPA_KEYWORDS),
    ("Dome  (tholos/dome/domed/spherical)", DOME_KEYWORDS),
]:
    entries = [e for e in selected if any(kw in e["keywords"] for kw in stage_kws)]
    ns, nap, na, p_s, enr = stage_stats(entries, stage_name)
    rate_str = f"{100*nap/ns:.1f}%" if ns > 0 else "—"
    print(f"  {stage_name:<35}  {ns:>5}  {nap:>4}  {na:>4}  {enr:>6.2f}×  {p_s:>8.4f}  {sig(p_s)}")

# ── Incremental contribution of mound keywords ────────────────────────────────
print()
print(SEP)
print("  INCREMENTAL CONTRIBUTION OF MOUND KEYWORDS")
print(SEP)
print()
N_do = len(dome_only)
N_mo = len(mound_only)
N_ov = len(overlap)
n_do_ap = sum(1 for e in dome_only  if e["ap"])
n_mo_ap = sum(1 for e in mound_only if e["ap"])
n_ov_ap = sum(1 for e in overlap    if e["ap"])

print(f"  Dome/stupa-only (equivalent to Test 2 raw):  N={N_do:>4},  A+={n_do_ap}")
print(f"  Mound-only (unique to evolution test):       N={N_mo:>4},  A+={n_mo_ap}")
print(f"  Overlap (both mound + dome/stupa kw):        N={N_ov:>4},  A+={n_ov_ap}")
print(f"  Combined total:                              N={N:>4},  A+={n_ap}")

if N_mo > 0:
    p_mo = binomtest(n_mo_ap, N_mo, P_NULL_AP, alternative="greater").pvalue
    enr_mo = (n_mo_ap / N_mo) / P_NULL_AP
    print(f"\n  Mound-only binomial: {n_mo_ap}/{N_mo} = {100*n_mo_ap/N_mo:.1f}%,  "
          f"{enr_mo:.2f}× vs null,  p = {p_mo:.4f}  {sig(p_mo)}")

if N_mo > 0 and N_do > 0:
    t2_non_ap = N_do - n_do_ap
    mo_non_ap = N_mo - n_mo_ap
    _, p_fish = fisher_exact([[n_mo_ap, mo_non_ap], [n_do_ap, t2_non_ap]], alternative="two-sided")
    print(f"  Fisher exact (mound-only vs dome-only A+ rate): p = {p_fish:.4f}  {sig(p_fish)}")

# ── Tier-A+ hits ──────────────────────────────────────────────────────────────
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

# ── Chronological breakdown ───────────────────────────────────────────────────
print()
print(SEP)
print("  CHRONOLOGICAL BREAKDOWN")
print(SEP)
dated  = [(e, e["date_yr"]) for e in selected if e["date_yr"] is not None]
eras   = [
    ("Deep antiquity (pre-1000 BCE)",  lambda y: y < -1000),
    ("Classical (1000 BCE – 1 CE)",    lambda y: -1000 <= y < 1),
    ("Late antiquity (1–500 CE)",      lambda y: 1 <= y < 500),
    ("Medieval (500–1500 CE)",         lambda y: 500 <= y < 1500),
    ("Early modern (1500–1800 CE)",    lambda y: 1500 <= y < 1800),
    ("Modern (post-1800 CE)",          lambda y: y >= 1800),
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

# ── Comparison vs context-validated version ───────────────────────────────────
print()
print(SEP)
print("  COMPARISON: RAW SWEEP (Test 2b primary) vs CONTEXT-VALIDATED (Exploratory)")
print(SEP)
print()
print(f"  {'':35}  {'Raw sweep':>12}  {'Validated':>12}")
print(f"  {'─'*62}")
print(f"  {'N (population)':<35}  {N:>12}  {'104':>12}")
print(f"  {'n A+ hits':<35}  {n_ap:>12}  {'13':>12}")
print(f"  {'A+ rate':<35}  {100*n_ap/N:>11.1f}%  {'12.5%':>12}")
print(f"  {'Binomial p (A+)':<35}  {bt_ap.pvalue:>12.4f}  {'0.0003':>12}")
print(f"  {f'Enrichment vs {P_NULL_AP:.0%} null':<35}  {enr_ap:>11.2f}×  {'3.12×':>12}")
print(f"  {'Anchor pctile (A+)':<35}  {pctile:>11.0f}th  ")

# ── Audit ─────────────────────────────────────────────────────────────────────
if "--audit" in sys.argv:
    print()
    print(SEP)
    print(f"  FULL AUDIT — ALL {N} SITES")
    print(SEP)
    print()
    for e in sorted(selected, key=lambda x: x["dev"]):
        stages = "+".join(sorted(e["stages"]))
        kws    = ", ".join(e["keywords"])
        date_s = e["date_label"] or "?"
        print(f"  [{e['tier']:>3}] {e['dev_km']:>6.1f} km  {e['name']:<50}  "
              f"[{stages:<12}]  kw={kws}  ~{date_s}")

# ── LaTeX macros ──────────────────────────────────────────────────────────────
print()
print(SEP)
print("  LATEX MACROS  (raw sweep — Test 2b primary)")
print(SEP)
print()

n_mound_stage    = sum(1 for e in selected if "mound" in e["stages"])
n_stupa_stage    = sum(1 for e in selected if "stupa" in e["stages"])
n_dome_stage     = sum(1 for e in selected if "dome"  in e["stages"])
n_mound_stage_ap = sum(1 for e in selected if "mound" in e["stages"] and e["ap"])
n_stupa_stage_ap = sum(1 for e in selected if "stupa" in e["stages"] and e["ap"])
n_dome_stage_ap  = sum(1 for e in selected if "dome"  in e["stages"] and e["ap"])
n_mound_stage_a  = sum(1 for e in selected if "mound" in e["stages"] and is_a_or_better(e["tier"]))
n_stupa_stage_a  = sum(1 for e in selected if "stupa" in e["stages"] and is_a_or_better(e["tier"]))
n_dome_stage_a   = sum(1 for e in selected if "dome"  in e["stages"] and is_a_or_better(e["tier"]))

p_mound_stage = binomtest(n_mound_stage_ap, n_mound_stage, P_NULL_AP, "greater").pvalue if n_mound_stage else 1.0
p_stupa_stage = binomtest(n_stupa_stage_ap, n_stupa_stage, P_NULL_AP, "greater").pvalue if n_stupa_stage else 1.0
p_dome_stage  = binomtest(n_dome_stage_ap,  n_dome_stage,  P_NULL_AP, "greater").pvalue if n_dome_stage  else 1.0

# Fisher exact (one-sided: greater) — stage vs full corpus background
def _fisher_ap(n_stage, n_stage_ap):
    table = [[n_stage_ap,            n_stage - n_stage_ap],
             [K_AP_CORPUS - n_stage_ap, N_CORPUS - n_stage - (K_AP_CORPUS - n_stage_ap)]]
    _, p = fisher_exact(table, alternative="greater")
    return p

def _fisher_a(n_stage, n_stage_a):
    table = [[n_stage_a,            n_stage - n_stage_a],
             [K_A_CORPUS - n_stage_a, N_CORPUS - n_stage - (K_A_CORPUS - n_stage_a)]]
    _, p = fisher_exact(table, alternative="greater")
    return p

fp_mound_ap = _fisher_ap(n_mound_stage, n_mound_stage_ap)
fp_stupa_ap = _fisher_ap(n_stupa_stage, n_stupa_stage_ap)
fp_dome_ap  = _fisher_ap(n_dome_stage,  n_dome_stage_ap)
fp_mound_a  = _fisher_a(n_mound_stage,  n_mound_stage_a)
fp_stupa_a  = _fisher_a(n_stupa_stage,  n_stupa_stage_a)
fp_dome_a   = _fisher_a(n_dome_stage,   n_dome_stage_a)
# Combined
fp_all_ap   = _fisher_ap(N, n_ap)
fp_all_a    = _fisher_a(N, n_a)
# A binomial combined (already computed as bt_a.pvalue)

print(f"  \\newcommand{{\\NevoTotal}}{{{N}}}            % total dome-evolution corpus")
print(f"  \\newcommand{{\\NevoApp}}{{{n_app}}}            % A++ sites in dome-evolution corpus")
print(f"  \\newcommand{{\\evoAppRate}}{{{100*n_app/N:.1f}}}          % A++ rate, dome-evolution corpus (%)")
print(f"  \\newcommand{{\\evoEnrichApp}}{{{enr_app:.2f}}}         % enrichment ratio, A++ in dome corpus")
print(f"  \\newcommand{{\\pEvoApp}}{{{bt_app.pvalue:.4f}}}         % p-value, A++ binomial (dome corpus)")
print(f"  \\newcommand{{\\NevoAp}}{{{n_ap}}}             % A+ sites in dome-evolution corpus")
print(f"  \\newcommand{{\\evoApRate}}{{{100*n_ap/N:.1f}}}           % A+ rate, dome-evolution corpus (%)")
print(f"  \\newcommand{{\\evoEnrichAp}}{{{enr_ap:.2f}}}          % enrichment ratio, A+ in dome corpus")
print(f"  \\newcommand{{\\pEvoAp}}{{{bt_ap.pvalue:.4f}}}          % p-value, A+ binomial (dome corpus)")
print(f"  \\newcommand{{\\pEvoApFisher}}{{{fp_all_ap:.4f}}}      % p-value, A+ Fisher exact (dome corpus)")
# Fisher ORs for combined
_evo_or_ap, _ = fisher_exact([[n_ap, N - n_ap], [K_AP_CORPUS - n_ap, N_CORPUS - N - (K_AP_CORPUS - n_ap)]], alternative="greater")
_evo_or_a,  _ = fisher_exact([[n_a,  N - n_a],  [K_A_CORPUS  - n_a,  N_CORPUS - N - (K_A_CORPUS  - n_a)]],  alternative="greater")
print(f"  \\newcommand{{\\evoApFisherOR}}{{{_evo_or_ap:.2f}}}    % OR, A+ Fisher exact (dome corpus)")
print(f"  \\newcommand{{\\pEvoAFisher}}{{{fp_all_a:.4f}}}        % p-value, A Fisher exact (dome corpus)")
print(f"  \\newcommand{{\\evoAFisherOR}}{{{_evo_or_a:.2f}}}      % OR, A Fisher exact (dome corpus)")
print(f"  \\newcommand{{\\NevoA}}{{{n_a}}}             % Tier-A sites in dome-evolution corpus")
print(f"  \\newcommand{{\\evoARate}}{{{100*n_a/N:.1f}}}           % A rate, dome-evolution corpus (%)")
print(f"  \\newcommand{{\\evoEnrichA}}{{{enr_a:.2f}}}           % enrichment ratio, A in dome corpus")
print(f"  \\newcommand{{\\pEvoA}}{{{bt_a.pvalue:.4f}}}          % p-value, A binomial (dome corpus)")
print(f"  \\newcommand{{\\evoExpAp}}{{{N*P_NULL_AP:.2f}}}          % expected A+ count at {P_NULL_AP:.0%} null (dome-evolution corpus)")
print(f"  \\newcommand{{\\evoExpA}}{{{N*P_NULL_A:.2f}}}           % expected A count at {P_NULL_A:.0%} null (dome-evolution corpus)")
print(f"  \\newcommand{{\\NevoMound}}{{{n_mound_stage}}}           % sites with mound stage")
print(f"  \\newcommand{{\\NevoMoundAp}}{{{n_mound_stage_ap}}}            % A+ sites with mound stage")
print(f"  \\newcommand{{\\evoMoundApRate}}{{{100*n_mound_stage_ap/n_mound_stage:.1f}}}           % A+ rate, mound stage (%)")
print(f"  \\newcommand{{\\NevoMoundA}}{{{n_mound_stage_a}}}            % A-tier sites with mound stage")
print(f"  \\newcommand{{\\evoMoundARate}}{{{100*n_mound_stage_a/n_mound_stage:.1f}}}           % A-tier rate, mound stage (%)")
print(f"  \\newcommand{{\\evoMoundEnrich}}{{{(n_mound_stage_ap/n_mound_stage)/P_NULL_AP:.2f}}}          % enrichment ratio, mound stage")
print(f"  \\newcommand{{\\pEvoMound}}{{{p_mound_stage:.4f}}}          % p-value, A+ binomial in mound stage")
print(f"  \\newcommand{{\\pEvoMoundFisher}}{{{fp_mound_ap:.4f}}}      % p-value, A+ Fisher exact in mound stage")
_mound_or_ap, _ = fisher_exact([[n_mound_stage_ap, n_mound_stage - n_mound_stage_ap], [K_AP_CORPUS - n_mound_stage_ap, N_CORPUS - n_mound_stage - (K_AP_CORPUS - n_mound_stage_ap)]], alternative="greater")
_mound_or_a,  _ = fisher_exact([[n_mound_stage_a,  n_mound_stage - n_mound_stage_a],  [K_A_CORPUS  - n_mound_stage_a,  N_CORPUS - n_mound_stage - (K_A_CORPUS  - n_mound_stage_a)]],  alternative="greater")
print(f"  \\newcommand{{\\evoMoundFisherOR}}{{{_mound_or_ap:.2f}}}    % OR, A+ Fisher exact in mound stage")
print(f"  \\newcommand{{\\pEvoMoundAFisher}}{{{fp_mound_a:.4f}}}      % p-value, A Fisher exact in mound stage")
print(f"  \\newcommand{{\\evoMoundAFisherOR}}{{{_mound_or_a:.2f}}}    % OR, A Fisher exact in mound stage")
print(f"  \\newcommand{{\\NevoStupa}}{{{n_stupa_stage}}}            % sites with stupa stage")
print(f"  \\newcommand{{\\NevoStupaAp}}{{{n_stupa_stage_ap}}}             % A+ sites with stupa stage")
print(f"  \\newcommand{{\\evoStupaApRate}}{{{100*n_stupa_stage_ap/n_stupa_stage:.1f}}}            % A+ rate, stupa stage (%)")
print(f"  \\newcommand{{\\NevoStupaA}}{{{n_stupa_stage_a}}}             % A-tier sites with stupa stage")
print(f"  \\newcommand{{\\evoStupaARate}}{{{100*n_stupa_stage_a/n_stupa_stage:.1f}}}            % A-tier rate, stupa stage (%)")
_stupa_or_ap, _ = fisher_exact([[n_stupa_stage_ap, n_stupa_stage - n_stupa_stage_ap], [K_AP_CORPUS - n_stupa_stage_ap, N_CORPUS - n_stupa_stage - (K_AP_CORPUS - n_stupa_stage_ap)]], alternative="greater")
_dome_or_ap,  _ = fisher_exact([[n_dome_stage_ap,  n_dome_stage  - n_dome_stage_ap],  [K_AP_CORPUS - n_dome_stage_ap,  N_CORPUS - n_dome_stage  - (K_AP_CORPUS - n_dome_stage_ap)]],  alternative="greater")
print(f"  \\newcommand{{\\evoStupaEnrich}}{{{(n_stupa_stage_ap/n_stupa_stage)/P_NULL_AP:.2f}}}           % enrichment ratio, stupa stage")
print(f"  \\newcommand{{\\pEvoStupa}}{{{p_stupa_stage:.4f}}}          % p-value, A+ binomial in stupa stage")
print(f"  \\newcommand{{\\pEvoStupaFisher}}{{{fp_stupa_ap:.4f}}}      % p-value, A+ Fisher exact in stupa stage")
print(f"  \\newcommand{{\\evoStupaFisherOR}}{{{_stupa_or_ap:.2f}}}    % OR, A+ Fisher exact in stupa stage")
print(f"  \\newcommand{{\\pEvoStupaAFisher}}{{{fp_stupa_a:.4f}}}      % p-value, A Fisher exact in stupa stage")
_stupa_or_a, _ = fisher_exact([[n_stupa_stage_a, n_stupa_stage - n_stupa_stage_a], [K_A_CORPUS - n_stupa_stage_a, N_CORPUS - n_stupa_stage - (K_A_CORPUS - n_stupa_stage_a)]], alternative="greater")
print(f"  \\newcommand{{\\evoStupaAFisherOR}}{{{_stupa_or_a:.2f}}}    % OR, A Fisher exact in stupa stage")
print(f"  \\newcommand{{\\NevoDome}}{{{n_dome_stage}}}            % sites with dome stage")
print(f"  \\newcommand{{\\NevoDomeAp}}{{{n_dome_stage_ap}}}             % A+ sites with dome stage")
print(f"  \\newcommand{{\\evoDomeApRate}}{{{100*n_dome_stage_ap/n_dome_stage:.1f}}}            % A+ rate, dome stage (%)")
print(f"  \\newcommand{{\\NevoDomeA}}{{{n_dome_stage_a}}}             % A-tier sites with dome stage")
print(f"  \\newcommand{{\\evoDomeARate}}{{{100*n_dome_stage_a/n_dome_stage:.1f}}}            % A-tier rate, dome stage (%)")
print(f"  \\newcommand{{\\evoDomeEnrich}}{{{(n_dome_stage_ap/n_dome_stage)/P_NULL_AP:.2f}}}           % enrichment ratio, dome stage")
print(f"  \\newcommand{{\\pEvoDome}}{{{p_dome_stage:.4f}}}          % p-value, A+ binomial in dome stage")
print(f"  \\newcommand{{\\pEvoDomeFisher}}{{{fp_dome_ap:.4f}}}       % p-value, A+ Fisher exact in dome stage")
print(f"  \\newcommand{{\\evoDomeFisherOR}}{{{_dome_or_ap:.2f}}}     % OR, A+ Fisher exact in dome stage")
print(f"  \\newcommand{{\\pEvoDomeAFisher}}{{{fp_dome_a:.4f}}}       % p-value, A Fisher exact in dome stage")
_dome_or_a, _ = fisher_exact([[n_dome_stage_a, n_dome_stage - n_dome_stage_a], [K_A_CORPUS - n_dome_stage_a, N_CORPUS - n_dome_stage - (K_A_CORPUS - n_dome_stage_a)]], alternative="greater")
print(f"  \\newcommand{{\\evoDomeAFisherOR}}{{{_dome_or_a:.2f}}}     % OR, A Fisher exact in dome stage")
print(f"  \\newcommand{{\\NevoMoundOnly}}{{{len(mound_only)}}}          % mound-only sites (no later dome/stupa)")
print(f"  \\newcommand{{\\NevoMoundOnlyAp}}{{{n_mo_ap}}}           % A+ in mound-only sites")
print(f"  \\newcommand{{\\NevoOverlap}}{{{len(overlap)}}}            % sites with both mound and dome/stupa stages")

# ── Write to results store ────────────────────────────────────────────────────
ResultsStore().write_many({
    "pEvoApp":          bt_app.pvalue,   # binomial p, A++ (dome-evolution corpus)
    "pEvoAp":           bt_ap.pvalue,    # binomial p, A+ (dome-evolution corpus) — Test 2b
    "pEvoA":            bt_a.pvalue,     # binomial p, A
    "pEvoApFisher":     fp_all_ap,       # Fisher p, A+ vs full corpus
    "pEvoAFisher":      fp_all_a,        # Fisher p, A vs full corpus
    "pEvoMound":        p_mound_stage,   # A+ binomial, mound stage
    "pEvoStupa":        p_stupa_stage,   # A+ binomial, stupa stage
    "pEvoDome":         p_dome_stage,    # A+ binomial, dome stage
    "pEvoMoundFisher":  fp_mound_ap,     # Fisher p, A+, mound stage
    "pEvoMoundAFisher": fp_mound_a,      # Fisher p, A, mound stage
    "pEvoStupaFisher":  fp_stupa_ap,     # Fisher p, A+, stupa stage
    "pEvoDomeFisher":   fp_dome_ap,      # Fisher p, A+, dome stage
    "NevoTotal":        N,               # total dome-evolution corpus
    "NevoAp":           n_ap,            # A+ sites
})
