import re
import sys
import numpy as np
from pathlib import Path
from scipy.stats import binomtest, chisquare

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus
from lib.results_store import ResultsStore
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_A_MAX, TIER_B_MAX,
    P_NULL_AP, P_NULL_A,
    deviation as _beru_dev, tier_label, is_aplus, is_a_or_better,
)
from lib.stats import significance_label as sig

# (spherical_monument_test.py, simulation_null_model.py, unit_sweep_fill.py)
from lib.dome_filter import (
    UNAMBIGUOUS_KEYWORDS, AMBIGUOUS_KEYWORDS, FORM_KEYWORDS, FORM_KEYWORD_RES,
    NEGATIVE_CONTEXT, NEGATIVE_CONTEXT_RES,
    POSITIVE_CONTEXT, POSITIVE_CONTEXT_RES,
    validate_ambiguous_sentence, validate_keyword_match,
)

# longitude or beru deviation.
#
MANUAL_EXCLUDE = {}  # no exclusions

# ── Helpers ──────────────────────────────────────────────────────────────────
def strip_html(text: str) -> str:
    return re.sub(r"<[^>]+>", " ", text or "")

def beru_deviation(lon: float):
    arc      = abs(lon - GERIZIM)
    beru_val = arc / BERU
    nearest  = round(beru_val * 10) / 10
    dev      = _beru_dev(lon)
    return arc, beru_val, nearest, dev

def form_sentences(desc: str, keywords: list) -> list:
    sentences = re.split(r"(?<=[.!?])\s+", desc)
    out = []
    for s in sentences:
        if any(FORM_KEYWORD_RES[k].search(s) for k in keywords):
            out.append(s.strip())
    return out

def load_sites(audit_mode: bool = False):
    corpus = load_corpus()
    total = len(corpus)

    sites    = []
    excluded = []
    context_rejected = []
    no_match = 0

    for site_obj in corpus:
        name = site_obj.site
        cat  = site_obj.category
        desc = site_obj.short_description

        # F1: Cultural or Mixed only
        if cat == "Natural":
            no_match += 1
            continue

        # F2: must have coordinates
        if not site_obj.has_coords:
            no_match += 1
            continue

        lat = site_obj.latitude
        lon = site_obj.longitude

        full_text = site_obj.full_text  # lowercased already

        all_text = desc
        if site_obj.extended_description:
            all_text += " " + site_obj.extended_description

        # Stage 1: raw keyword match
        raw_matched = [k for k in FORM_KEYWORDS if FORM_KEYWORD_RES[k].search(full_text)]
        if not raw_matched:
            no_match += 1
            continue

        # Stage 2: context validation
        validated_keywords = []
        validation_notes = {}
        all_valid_sentences = []
        all_rejected_sentences = []
        
        for kw in raw_matched:
            is_valid, sentences, note = validate_keyword_match(all_text, kw)
            validation_notes[kw] = note
            if is_valid:
                validated_keywords.append(kw)
                all_valid_sentences.extend(sentences)
        
        if not validated_keywords:
            context_rejected.append({
                "name":    name,
                "matched": raw_matched,
                "notes":   validation_notes,
                "lon":     lon,
            })
            continue

        # F3: mechanical exclusions
        if name in MANUAL_EXCLUDE:
            excluded.append({
                "name":    name,
                "matched": validated_keywords,
                "reason":  MANUAL_EXCLUDE[name],
                "lon":     lon,
            })
            continue

        arc, beru_val, nearest, dev = beru_deviation(lon)
        sentences = form_sentences(all_text, validated_keywords)
        
        xml_text = (desc + " " + name).lower()
        xml_matched = [k for k in FORM_KEYWORDS if FORM_KEYWORD_RES[k].search(xml_text)]
        xml_validated = []
        xml_text_full = desc + " " + name
        for kw in xml_matched:
            is_valid, _, _ = validate_keyword_match(xml_text_full, kw)
            if is_valid:
                xml_validated.append(kw)
        match_source = "xml" if xml_validated else "extended"

        sites.append({
            "name":      name,
            "lat":       lat,
            "lon":       lon,
            "cat":       cat,
            "keywords":  validated_keywords,
            "raw_keywords": raw_matched,
            "validation": validation_notes,
            "sentences": sentences,
            "match_source": match_source,
            "arc":       arc,
            "beru":      beru_val,
            "nearest":   nearest,
            "dev":       dev,
            "dev_km":    dev * BERU * 111.0,
            "tier":      tier_label(dev),
        })

    return sites, excluded, context_rejected, total

# ── Audit mode ───────────────────────────────────────────────────────────────
if "--audit" in sys.argv:
    sites, excluded, context_rejected, total = load_sites(audit_mode=True)
    print(f"\nUNESCO XML KEYWORD AUDIT — spherical_monument_test.py")
    print(f"Keyword set: {FORM_KEYWORDS}")
    print(f"Unambiguous: {UNAMBIGUOUS_KEYWORDS}")
    print(f"Ambiguous (context-validated): {AMBIGUOUS_KEYWORDS}")
    print(f"Total entries in XML : {total}")
    print(f"Sites INCLUDED       : {len(sites)}")
    print(f"Sites CONTEXT-REJECTED: {len(context_rejected)}")
    print(f"Sites EXCLUDED       : {len(excluded)}")
    print()
    print("INCLUDED SITES:")
    print("─" * 100)
    for s in sorted(sites, key=lambda x: x["name"]):
        print(f"  {s['name']}")
        print(f"    KEYS   : {s['keywords']}  (raw: {s['raw_keywords']})")
        print(f"    VALID  : {s['validation']}")
        for sent in s["sentences"][:2]:
            print(f"    FORM   : \"{sent[:200]}\"")
        print()
    print("CONTEXT-REJECTED SITES:")
    print("─" * 100)
    for e in context_rejected:
        print(f"  {e['name']}")
        print(f"    RAW KEYS : {e['matched']}")
        print(f"    NOTES    : {e['notes']}")
        print()
    print("EXCLUDED SITES (with reasons):")
    print("─" * 100)
    for e in excluded:
        print(f"  {e['name']}")
        print(f"    KEYS   : {e['matched']}")
        print(f"    REASON : {e['reason']}")
        print()
    sys.exit(0)

# ── Main analysis ────────────────────────────────────────────────────────────
sites, excluded, context_rejected, total_xml = load_sites()

print("=" * 100)
print("  UNESCO WHC — DOMED / SPHERICAL / STUPA / THOLOS MONUMENT TEST")
print("  MINIMAL MORPHOLOGICAL KEYWORD SET: {stupa, tholos, dome, domed, domes, spherical}")
print(f"  Anchor: Gerizim {GERIZIM}°E  |  BERU={BERU}°  |  "
      f"Tier-A+ ≤ {TIER_APLUS} beru (≤ {TIER_APLUS*BERU*111:.0f} km)  |  "
      f"Tier-A ≤ {TIER_A_MAX} beru (≤ {TIER_A_MAX*BERU*111:.0f} km)")
print(f"  Source: UNESCO WHC XML + Extended OUV Descriptions ({total_xml} total entries)")
print(f"  Keyword filter: {len(FORM_KEYWORDS)} morphological keywords")
print(f"  Mechanical exclusions: {len(excluded)} sites (documented in --audit output)")
print(f"  Context-rejected: {len(context_rejected)} sites (ambiguous keyword, no architectural context)")
n_xml = sum(1 for s in sites if s["match_source"] == "xml")
n_ext = sum(1 for s in sites if s["match_source"] == "extended")
print(f"  Matched via XML: {n_xml}  |  Matched via extended OUV text only: {n_ext}")
print("=" * 100)

print(f"""
─────────────────────────────────────────────────────────────────────────
           (site name + short_description + extended OUV web-page text).
           Keywords: {FORM_KEYWORDS}
           sentence-by-sentence against architectural context lists.
           {len(context_rejected)} sites rejected (keyword in non-architectural context).
  Step 2b. Mechanical exclusions applied: {len(excluded)} sites removed.
  
        the full Outstanding Universal Value text, which includes
        keywords like "stupa" and "dome" that often do not appear in
        the abbreviated XML short_description.
        e.g., Kathmandu Valley's XML has no "stupa" but the OUV text
        describes "the great stupas of Swayambhunath and Bauddhanath".

  {len(context_rejected)} context-rejected sites (keyword matched but no architectural context):""")
for cr in context_rejected:
    print(f"    ✗  {cr['name'][:70]}")
    print(f"       Raw keywords: {cr['matched']}  Notes: {cr['notes']}")

if excluded:
    print(f"\n  {len(excluded)} excluded sites:")
    for e in excluded:
        print(f"    ✗  {e['name'][:70]}")
        print(f"       Reason: {e['reason'][:90]}")

print(f"""
  {len(sites)} INCLUDED sites:""")

# ── Per-site table ────────────────────────────────────────────────────────────
print()
print("─" * 130)
print(f"  {'Site':<52s}  {'Lon':>8}  {'Arc°':>6}  {'Beru':>7}  "
      f"{'Near':>5}  {'Dev':>7}  {'km':>6}  T  Src  Keywords")
print("─" * 130)
for s in sorted(sites, key=lambda x: x["dev"]):
    mark = " ◀◀ A+" if is_aplus(s["tier"]) else (" ◀ A" if s["tier"] == "A" else "")
    kw_short = ", ".join(s["keywords"])[:35]
    src = "ext" if s["match_source"] == "extended" else "xml"
    print(
        f"  {s['name']:<52s}  {s['lon']:>8.4f}  {s['arc']:>6.3f}  "
        f"{s['beru']:>7.4f}  {s['nearest']:>5.1f}  {s['dev']:>7.5f}  "
        f"{s['dev_km']:>6.1f}  {s['tier']}{mark}  {src}  [{kw_short}]"
    )

# ── Statistical tests ─────────────────────────────────────────────────────────
N        = len(sites)
n_Ap     = sum(1 for s in sites if is_aplus(s["tier"]))
n_A      = sum(1 for s in sites if is_a_or_better(s["tier"]))
n_B      = sum(1 for s in sites if s["tier"] == "B")
n_C      = sum(1 for s in sites if s["tier"] == "C")
mean_dev = sum(s["dev"] for s in sites) / N if N else 0

bt_Ap = binomtest(n_Ap, N, P_NULL_AP, alternative="greater")
bt_A  = binomtest(n_A,  N, P_NULL_A,  alternative="greater")

obs = [0] * 5
for s in sites:
    obs[min(int(s["dev"] / 0.010), 4)] += 1
chi_stat, chi_p = chisquare(obs, f_exp=[N / 5.0] * 5)

print()
print("=" * 100)
print("  STATISTICAL TESTS")
print("  H₀: domed-monument longitudes are UNIFORM w.r.t. 0.1-beru harmonics")
print("  H₁: hit rate > geometric null  (one-tailed binomial)")
print("  Null rates are PURELY GEOMETRIC — no archaeological priors")
print("=" * 100)
print()
print(f"  {'Metric':<45s}  {'Obs':>5}  {'Exp(H₀)':>9}  {'Ratio':>7}  {'p':>8}  Sig")
print(f"  {'-'*90}")
print(f"  {'Tier-A+ (≤0.002 beru, ≤6.7 km)':<45s}  {n_Ap:>5d}  "
      f"{N*P_NULL_AP:>9.2f}  {n_Ap/max(N*P_NULL_AP,0.001):>6.2f}×  "
      f"{bt_Ap.pvalue:>8.4f}  {sig(bt_Ap.pvalue)}")
print(f"  {'Tier-A  (≤0.010 beru, ≤33 km)  ◀ PRIMARY':<45s}  {n_A:>5d}  "
      f"{N*P_NULL_A:>9.2f}  {n_A/max(N*P_NULL_A,0.001):>6.2f}×  "
      f"{bt_A.pvalue:>8.4f}  {sig(bt_A.pvalue)}")
print(f"  {'Tier-B  (≤0.050 beru)':<45s}  {n_B:>5d}")
print(f"  {'Tier-C  (>0.050 beru)':<45s}  {n_C:>5d}")
print(f"  {'Mean deviation':<45s}  {mean_dev:>5.4f} beru  ({mean_dev*BERU*111:.1f} km)")
print(f"  {'χ²-uniform (5 bins, 0–0.05 beru, df=4)':<45s}  "
      f"{'':>5}  {'':>9}  {'':>7}  {chi_p:>8.4f}  {sig(chi_p)}")

# Histogram
print()
bins = ["0–0.01", "0.01–0.02", "0.02–0.03", "0.03–0.04", "0.04–0.05"]
for lbl, o in zip(bins, obs):
    bar  = "█" * int(round(o * 40 / max(obs + [1])))
    note = "  ← Tier-A  (A+ = leftmost 20% of this bin)" if lbl == "0–0.01" else ""
    print(f"    {lbl}:  {o:3d}  ({100*o/N if N else 0:>4.0f}%)  {bar}{note}")

print(f"\n  Tier-A+ hits ({n_Ap}) — form sentence from UNESCO description:")
for s in sorted([s for s in sites if is_aplus(s["tier"])], key=lambda x: x["dev"]):
    print(f"    {s['name']:<52s}  dev={s['dev']:.5f}  {s['dev_km']:>5.1f} km  "
          f"lon={s['lon']:.4f}°  beru→{s['nearest']:.1f}")
    for sent in s["sentences"][:1]:
        print(f"      UNESCO: \"{sent[:180]}\"")

print(f"\n  Tier-A hits ({n_A - n_Ap}, excl. A+) — form sentence from UNESCO description:")
for s in sorted([s for s in sites if s["tier"] == "A"], key=lambda x: x["dev"]):
    print(f"    {s['name']:<52s}  dev={s['dev']:.5f}  {s['dev_km']:>5.1f} km  "
          f"lon={s['lon']:.4f}°  beru→{s['nearest']:.1f}")
    for sent in s["sentences"][:1]:
        print(f"      UNESCO: \"{sent[:180]}\"")

# ── Anchor sweep ──────────────────────────────────────────────────────────────
print()
print("=" * 100)
print("  ANCHOR SWEEP  34.0–37.0°E @ 0.001° resolution")
print("  Tests whether Gerizim is a special anchor vs any random longitude")
print("=" * 100)

sweep = np.arange(34.0, 37.001, 0.001)

def count_at(anchor, pop, threshold):
    n = 0
    for s in pop:
        arc  = abs(s["lon"] - anchor)
        bv   = arc / BERU
        near = round(bv * 10) / 10
        if abs(bv - near) <= threshold:
            n += 1
    return n

cnt_Ap = np.array([count_at(a, sites, TIER_APLUS) for a in sweep])
cnt_A  = np.array([count_at(a, sites, TIER_A_MAX)  for a in sweep])

g_Ap = count_at(GERIZIM, sites, TIER_APLUS)
g_A  = count_at(GERIZIM, sites, TIER_A_MAX)
j_Ap = count_at(35.235,  sites, TIER_APLUS)   # Jerusalem
j_A  = count_at(35.235,  sites, TIER_A_MAX)

pct_Ap   = float(np.mean(cnt_Ap <= g_Ap))
pct_A    = float(np.mean(cnt_A  <= g_A))
j_pct_Ap = float(np.mean(cnt_Ap <= j_Ap))
j_pct_A  = float(np.mean(cnt_A  <= j_A))

gi = int(round((GERIZIM - 34.0) / 0.001))
W  = 20
lmax_Ap = int(np.max(cnt_Ap[max(0, gi-W): gi+W+1]))
lmax_A  = int(np.max(cnt_A [max(0, gi-W): gi+W+1]))
gmax_Ap = int(np.max(cnt_Ap)); gmax_A = int(np.max(cnt_A))
gmax_loc_Ap = float(sweep[int(np.argmax(cnt_Ap))])
gmax_loc_A  = float(sweep[int(np.argmax(cnt_A))])

is_lmax_Ap = "✓ local max" if g_Ap >= lmax_Ap else f"✗ (local max={lmax_Ap})"
is_lmax_A  = "✓ local max" if g_A  >= lmax_A  else f"✗ (local max={lmax_A})"

print(f"\n  {'Anchor':<26}  {'A+(≤6.7km)':>11}  {'Pctile':>7}  {'A(≤33km)':>9}  {'Pctile':>7}")
print(f"  {'-'*72}")
print(f"  {'Gerizim (35.272°E)':<26}  {g_Ap:>11}  {pct_Ap*100:>6.0f}th  "
      f"{g_A:>9}  {pct_A*100:>6.0f}th")
print(f"  {'Jerusalem (35.235°E)':<26}  {j_Ap:>11}  {j_pct_Ap*100:>6.0f}th  "
      f"{j_A:>9}  {j_pct_A*100:>6.0f}th")
print(f"  Global max in sweep: A+={gmax_Ap} at {gmax_loc_Ap:.3f}°E  |  "
      f"A={gmax_A} at {gmax_loc_A:.3f}°E")
print(f"  Gerizim local-max status (±{W*0.001:.3f}°): "
      f"A+ {is_lmax_Ap}  |  A {is_lmax_A}")

# ── Keyword composition ──────────────────────────────────────────────────────
print()
print("=" * 100)
print("  KEYWORD COMPOSITION")
print("=" * 100)
from collections import Counter
kw_counts = Counter()
for s in sites:
    for k in s["keywords"]:
        kw_counts[k] += 1
print()
for kw, count in kw_counts.most_common():
    print(f"    {kw:<30s}  {count} site(s)")

# ── Verdict ───────────────────────────────────────────────────────────────────
enr_Ap = (n_Ap / N) / P_NULL_AP if N else 0
enr_A  = (n_A  / N) / P_NULL_A  if N else 0

print()
print("=" * 100)
print("  VERDICT")
print("=" * 100)

print(f"""
  SELECTION INDEPENDENCE
  ──────────────────────────────────────────────────────────────────────────────
  All {N} sites selected by UNESCO designation (a UN body with zero reference
  to Mount Gerizim) plus keyword match on MORPHOLOGICAL FORM.
  The keyword set {{stupa, tholos, dome, domed, domes, spherical}} is the
  minimal English-language morphological set for hemispherical/domed
  monumental architecture — pre-specifiable from any dictionary.

  PRIMARY RESULT — Tier-A+ binomial (H₁: rate > {P_NULL_AP:.0%})
  ──────────────────────────────────────────────────────────────────────────────
  N = {N}  |  Tier-A+: {n_Ap}/{N} = {100*n_Ap/N:.1f}%  |  {enr_Ap:.2f}× vs null 4%
  p = {bt_Ap.pvalue:.4f}  {sig(bt_Ap.pvalue)}

  SECONDARY — Tier-A (H₁: rate > {P_NULL_A:.0%})
  ──────────────────────────────────────────────────────────────────────────────
  N = {N}  |  Tier-A: {n_A}/{N} = {100*n_A/N:.1f}%  |  {enr_A:.2f}× vs null 20%
  p = {bt_A.pvalue:.4f}  {sig(bt_A.pvalue)}

  χ²-UNIFORM — (5 bins, df=4)
  ──────────────────────────────────────────────────────────────────────────────
  p = {chi_p:.4f}  {sig(chi_p)}

  ──────────────────────────────────────────────────────────────────────────────
  With N ≈ {N}, the binomial tests have moderate power.
  all metrics (binomial, χ², enrichment ratio, anchor sweep).

  ANCHOR SWEEP
  ──────────────────────────────────────────────────────────────────────────────
  Gerizim percentile: A+ = {pct_Ap*100:.0f}th  |  A = {pct_A*100:.0f}th
""")

# ── LaTeX macros ──────────────────────────────────────────────────────────────
print("=" * 100)
print("  LATEX MACROS  (context-validated — Exploratory robustness check)")
print("=" * 100)
print()
print(f"  \\newcommand{{\\pCircApValidated}}{{{bt_Ap.pvalue:.4f}}}  % p-value, A+ binomial (context-validated dome corpus)")
print(f"  \\newcommand{{\\pCircAValidated}}{{{bt_A.pvalue:.4f}}}   % p-value, A binomial (context-validated dome corpus)")
print(f"  \\newcommand{{\\circEnrichApValidated}}{{{enr_Ap:.2f}}}  % enrichment ratio, A+ (context-validated dome corpus)")
print()

# ── Write to results store ────────────────────────────────────────────────────
ResultsStore().write_many({
    "pCircAp_validated": bt_Ap.pvalue,   # binomial p, A+ (context-validated dome corpus) — Exploratory 2x
    "pCircA_validated":  bt_A.pvalue,    # binomial p, A  (context-validated dome corpus)
})
