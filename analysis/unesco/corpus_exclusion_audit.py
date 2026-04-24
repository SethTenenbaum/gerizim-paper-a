"""
corpus_exclusion_audit.py
=========================
Audits the two classes of sites excluded from the working corpus:
  (1) Natural-only sites (category = "Natural")
  (2) Cultural/Mixed sites with missing geocoordinates

Generates LaTeX macros for inclusion in the paper and prints a
reviewer-facing summary of why the exclusions do not bias the result.
"""

import sys
from pathlib import Path
from scipy.stats import binomtest, fisher_exact

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus
from lib.beru import deviation, tier_label, is_aplus, is_aplusplus, P_NULL_AP
from lib.stats import significance_label as sig


def fmt_p(p):
    return "< 0.001" if p < 0.001 else f"{p:.3f}"


corpus = load_corpus()

# ── Category / coord breakdown ────────────────────────────────────────────────
all_sites       = corpus
cultural_mixed  = [s for s in corpus if s.category in ("Cultural", "Mixed")]
natural_only    = [s for s in corpus if s.category == "Natural"]

cm_with_coords  = [s for s in cultural_mixed if s.has_coords]   # working corpus
cm_no_coords    = [s for s in cultural_mixed if not s.has_coords]
nat_with_coords = [s for s in natural_only   if s.has_coords]
nat_no_coords   = [s for s in natural_only   if not s.has_coords]

# ── Natural site tier distribution ───────────────────────────────────────────
nat_ap  = [s for s in nat_with_coords if is_aplus(tier_label(deviation(s.longitude)))]
nat_app = [s for s in nat_with_coords if is_aplusplus(tier_label(deviation(s.longitude)))]

n_nat    = len(nat_with_coords)
n_nat_ap = len(nat_ap)
rate_nat = 100 * n_nat_ap / n_nat

# ── Cultural A+ rate ──────────────────────────────────────────────────────────
n_cult    = len(cm_with_coords)
n_cult_ap = sum(1 for s in cm_with_coords if is_aplus(tier_label(deviation(s.longitude))))
rate_cult = 100 * n_cult_ap / n_cult

# ── Fisher: natural vs cultural A+ ───────────────────────────────────────────
ct = [[n_nat_ap, n_nat - n_nat_ap], [n_cult_ap, n_cult - n_cult_ap]]
or_nat_cult, p_nat_cult = fisher_exact(ct)  # two-sided

# ── Sensitivity: include natural in Test 1 ───────────────────────────────────
all_with = [s for s in corpus if s.has_coords]
n_all    = len(all_with)
n_all_ap = sum(1 for s in all_with if is_aplus(tier_label(deviation(s.longitude))))
p_incl   = binomtest(n_all_ap, n_all, P_NULL_AP, alternative="greater").pvalue
p_excl   = binomtest(n_cult_ap, n_cult, P_NULL_AP, alternative="greater").pvalue

# ── Print summary ─────────────────────────────────────────────────────────────
SEP = "─" * 72
print(SEP)
print("  CORPUS EXCLUSION AUDIT")
print(SEP)
print(f"  Total XML sites:              {len(corpus)}")
print(f"  Cultural/Mixed:               {len(cultural_mixed)}")
print(f"    with coords (working):      {n_cult}")
print(f"    missing coords (excluded):  {len(cm_no_coords)}")
print(f"  Natural only:                 {len(natural_only)}")
print(f"    with coords:                {n_nat}")
print(f"    missing coords:             {len(nat_no_coords)}")
print()
print("  Cultural/Mixed sites EXCLUDED for missing coordinates:")
for s in cm_no_coords:
    print(f"    [{s.id_number}] {s.site} ({s.states})")
print()
print(f"  Natural site A+ rate:   {n_nat_ap}/{n_nat} = {rate_nat:.1f}%")
print(f"  Cultural site A+ rate:  {n_cult_ap}/{n_cult} = {rate_cult:.1f}%")
print(f"  Fisher OR (nat/cult):   {or_nat_cult:.2f},  p = {fmt_p(p_nat_cult)} {sig(p_nat_cult)}")
print()
print("  Natural A+ sites (all geographic features, no domed monuments):")
for s in sorted(nat_ap, key=lambda x: deviation(x.longitude)):
    d = deviation(s.longitude)
    t = tier_label(d)
    print(f"    [{t}] {s.site[:60]}  lon={s.longitude:.2f}")
print()
print("  Sensitivity — Test 1 result if Natural sites included:")
print(f"    Standard (Cultural/Mixed, N={n_cult}):  p = {fmt_p(p_excl)} {sig(p_excl)}")
print(f"    Including Natural  (N={n_all}):        p = {fmt_p(p_incl)} {sig(p_incl)}")
print(SEP)

# ── LaTeX macros ──────────────────────────────────────────────────────────────
print()
print("% LaTeX macros — corpus_exclusion_audit.py")
print(f"\\newcommand{{\\auditTotalXML}}{{{len(corpus)}}}")
print(f"\\newcommand{{\\auditNaturalN}}{{{n_nat}}}")
print(f"\\newcommand{{\\auditNaturalApN}}{{{n_nat_ap}}}")
print(f"\\newcommand{{\\auditNaturalApRate}}{{{rate_nat:.1f}}}")
print(f"\\newcommand{{\\auditNaturalFisherOR}}{{{or_nat_cult:.2f}}}")
print(f"\\newcommand{{\\auditNaturalFisherP}}{{{fmt_p(p_nat_cult)}}}")
print(f"\\newcommand{{\\auditCMNoCoordN}}{{{len(cm_no_coords)}}}")
print(f"\\newcommand{{\\auditTest1InclNaturalP}}{{{fmt_p(p_incl)}}}")
