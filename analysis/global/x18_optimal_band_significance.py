"""
x18_optimal_band_significance.py
GROUP: 11d

SIGNIFICANCE TEST FOR THE OPTIMAL x.18-DEGREE BAND
====================================================
The global anchor sweep (anchor_uniqueness_audit.py, GROUP 11) already
identifies the optimal x.18°E phase as producing the maximum observed A+
count across all 36,000 trial anchors.  This script reads those already-
computed values from the results store and runs one-sided binomial tests
against the 4% geometric null to report the significance of the optimal band.

No re-scanning is performed.  All inputs come from the results store.

Outputs macros:
  \\xOptBandApCount      — observed A+ at the x.18° optimum (Sweep A, N=1010)
  \\xOptBandApCountB     — observed A+ at the x.18° optimum (Sweep B, N=1011)
  \\xOptBandN            — corpus N for the binomial test (Sweep A)
  \\xOptBandBinomP       — one-sided binomial p (Sweep A)
  \\xOptBandBinomPB      — one-sided binomial p (Sweep B)
  \\xOptBandEnrich       — enrichment ratio (observed/expected) at optimum (Sweep A)
  \\xOptBandEnrichB      — enrichment ratio (Sweep B)
====================================================
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from scipy.stats import binomtest
from lib.results_store import ResultsStore
from lib.beru import P_NULL_AP

rs = ResultsStore()

# ── Read pre-computed values from store ──────────────────────────────────────
obs_A = int(rs.read("anchorSweepGlobalMaxA"))   # 59  (Sweep A, Jerusalem removed)
obs_B = int(rs.read("anchorSweepGlobalMaxB"))   # 60  (Sweep B, Jerusalem kept)
N_A   = int(rs.read("anchorSweepNsweep"))       # 1010

# Sweep B corpus is N_A + 1 (Jerusalem is in the corpus, self-excluded → 1010)
N_B   = N_A   # both sweeps exclude the focal anchor, leaving 1010 sites

# ── Binomial tests ────────────────────────────────────────────────────────────
result_A = binomtest(obs_A, N_A, P_NULL_AP, alternative="greater")
result_B = binomtest(obs_B, N_B, P_NULL_AP, alternative="greater")

p_A = float(result_A.pvalue)
p_B = float(result_B.pvalue)

expected   = P_NULL_AP * N_A
enrich_A   = round(obs_A / expected, 2)
enrich_B   = round(obs_B / expected, 2)

def _sig(p):
    return ("***" if p < 0.001 else
            "**"  if p < 0.01  else
            "*"   if p < 0.05  else
            "~"   if p < 0.10  else "ns")

print("=" * 70)
print("  x.18° OPTIMAL BAND — BINOMIAL SIGNIFICANCE TEST")
print(f"  Null rate P0 = {P_NULL_AP:.0%}  (geometric 4% null)")
print("=" * 70)
print(f"  Sweep A (Jerusalem removed, N={N_A})")
print(f"    Observed A+ at optimum : {obs_A}")
print(f"    Expected under null    : {expected:.1f}")
print(f"    Enrichment             : {enrich_A}x")
print(f"    Binomial p (one-sided) : {p_A:.6f}  {_sig(p_A)}")
print()
print(f"  Sweep B (Jerusalem kept, N={N_B})")
print(f"    Observed A+ at optimum : {obs_B}")
print(f"    Expected under null    : {expected:.1f}")
print(f"    Enrichment             : {enrich_B}x")
print(f"    Binomial p (one-sided) : {p_B:.6f}  {_sig(p_B)}")
print()

print("  LATEX MACROS (GROUP 11d -- x.18 optimal band significance):")
print("=" * 70)
macro_pairs = [
    ("xOptBandApCount",  f"{obs_A}"),
    ("xOptBandApCountB", f"{obs_B}"),
    ("xOptBandN",        f"{N_A}"),
    ("xOptBandBinomP",   f"{p_A:.4f}"),
    ("xOptBandBinomPB",  f"{p_B:.4f}"),
    ("xOptBandEnrich",   f"{enrich_A}"),
    ("xOptBandEnrichB",  f"{enrich_B}"),
]
for name, val in macro_pairs:
    print(f"\\newcommand{{\\{name}}}{{{val}}}")

ResultsStore().write_many({
    "xOptBandApCount":  obs_A,
    "xOptBandApCountB": obs_B,
    "xOptBandN":        N_A,
    "xOptBandBinomP":   p_A,
    "xOptBandBinomPB":  p_B,
    "xOptBandEnrich":   enrich_A,
    "xOptBandEnrichB":  enrich_B,
})
print("Results written to data/store/results.json")
