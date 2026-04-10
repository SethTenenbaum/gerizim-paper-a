"""
bonferroni_correction.py — Bonferroni correction for manuscript p-values.

USAGE
-----
    cd /path/to/gerizim-paper-a
    python3 analysis/unesco/bonferroni_correction.py

OUTPUT
------
    Plain-text table of Bonferroni-adjusted p-values for all pre-specified
    (confirmatory) tests, plus a summary section for exploratory tests.
    Prints LaTeX \\newcommand lines for each adjusted p-value.

HOW P-VALUES ARE SOURCED
------------------------
    The test registry (which tests are confirmatory, which store key holds
    each raw p-value) lives in config.json → "tests".

    Each analysis script writes its p-value to the shared results store
    (data/store/results.json) alongside its \\newcommand output.

    This script reads from that store — NO hardcoded p-values here.
    If a value is missing, the script tells you which analysis script to run.

TEST NUMBERING (from config.json)
----------------------------------
    Test 1   — Global Enrichment (full corpus A+)
    Test 2   — Domed and Spherical Monuments  [RAW sweep]
    Test 3   — Cluster Asymmetry
    Test 4   — Temporal Gradient (Cochran-Armitage)
    Test 2b  — Hemispherical Mound Evolution  [sensitivity variant, not in family]
    Test E   — Founding-Site Enrichment       [exploratory]
"""

from __future__ import annotations
import json
from pathlib import Path
import sys

# ── Load shared config and results store ─────────────────────────────────────
_REPO = Path(__file__).parent.parent.parent
sys.path.insert(0, str(_REPO))

from lib.results_store import ResultsStore

_CONFIG_PATH = _REPO / "config.json"
with open(_CONFIG_PATH) as f:
    CONFIG = json.load(f)

TESTS_CFG    = CONFIG["tests"]
FAMILY_ALPHA = TESTS_CFG["family_alpha"]
store        = ResultsStore()

# ── Build test lists from config (never hardcoded here) ───────────────────────
def _load_p(entry: dict) -> float:
    """Fetch p-value from results store; abort with a clear message if missing."""
    key = entry["p_key"]
    val = store.read(key, default=None)
    if val is None:
        print(f"\n  ERROR: results store missing key '{key}'")
        print(f"  Run:   python3 {entry['script']}")
        print(f"  Then re-run this script.\n")
        sys.exit(1)
    return float(val)

confirmatory = [
    (e["id"], e["label"], _load_p(e), e["macro"], e["script"])
    for e in TESTS_CFG["confirmatory"]
]
sensitivity  = [
    (e["id"], e["label"], _load_p(e), e["macro"], e["script"], e.get("note",""))
    for e in TESTS_CFG.get("sensitivity", [])
]
exploratory  = [
    (e["id"], e["label"], _load_p(e), e.get("macro"), e["script"], e.get("note",""))
    for e in TESTS_CFG.get("exploratory", [])
    if store.has(e["p_key"])   # skip if the script hasn't been run yet
]

K = len(confirmatory)
PER_TEST_THRESHOLD = FAMILY_ALPHA / K


def holm_adjusted(raw_ps: list[float]) -> list[float]:
    """
    Holm (1979) step-down procedure.
    Returns adjusted p-values in the ORIGINAL order.
    Step-down rule: sort ascending; multiply p_(i) by (k - i + 1);
    enforce monotonicity by taking running maximum; cap at 1.
    """
    n     = len(raw_ps)
    order = sorted(range(n), key=lambda i: raw_ps[i])
    p_adj = [0.0] * n
    running = 0.0
    for rank, idx in enumerate(order):
        val     = raw_ps[idx] * (n - rank)
        running = max(running, val)
        p_adj[idx] = min(running, 1.0)
    return p_adj


# ── Compute adjusted p-values ─────────────────────────────────────────────────
raw_ps   = [raw_p for _, _, raw_p, _, _ in confirmatory]
bonf_adj = [min(p * K, 1.0) for p in raw_ps]
holm_adj = holm_adjusted(raw_ps)

results = []
for i, (tid, lbl, raw_p, macro, script) in enumerate(confirmatory):
    p_bonf   = bonf_adj[i]
    p_holm   = holm_adj[i]
    survives = p_holm < FAMILY_ALPHA      # Holm is the primary decision rule
    results.append((tid, lbl, raw_p, p_bonf, p_holm, survives, macro, script))

# ── Pretty-print ──────────────────────────────────────────────────────────────
SEP = "=" * 76

print()
print(SEP)
print("  MULTIPLE-COMPARISONS — Gerizim/UNESCO Analysis (Bonferroni & Holm)")
print(SEP)
print(f"  Confirmatory family size : k = {K}")
print(f"  Family-wise alpha        : {FAMILY_ALPHA}")
print(f"  Per-test threshold       : alpha/k = {FAMILY_ALPHA}/{K} = {PER_TEST_THRESHOLD:.4f}")
print(f"  Primary procedure        : Holm (1979) step-down (controls FWER exactly)")
print(f"  Secondary (for ref)      : Bonferroni (conservative upper bound)")
print(f"  P-values sourced from    : data/store/results.json")
print()
print(f"  {'Test':<6} {'Raw p':>9} {'Bonf. adj':>10} {'Holm adj':>9}  {'Result (Holm)':<15}  Label")
print("  " + "-" * 80)
for tid, lbl, raw_p, p_bonf, p_holm, survives, macro, _ in results:
    tag = ("*** SURVIVES" if p_holm < 0.001
           else "** SURVIVES"  if p_holm < 0.01
           else "*  SURVIVES"  if p_holm < 0.05
           else "   ns")
    print(f"  {tid:<6} {raw_p:>9.4f} {p_bonf:>10.4f} {p_holm:>9.4f}  {tag:<15}  {lbl}")

print()
print(SEP)
print("  SENSITIVITY VARIANTS (excluded from Bonferroni family)")
print(SEP)
for tid, lbl, raw_p, macro, script, note in sensitivity:
    print(f"  Test {tid}:  {lbl}")
    print(f"          Raw p = {raw_p}  |  [SENSITIVITY — not an independent test]")
    if note:
        print(f"          Note: {note}")
    print(f"          Source: {script}")
print()

print(SEP)
print("  EXPLORATORY (excluded from Bonferroni family)")
print(SEP)
for tid, lbl, raw_p, macro, script, note in exploratory:
    print(f"  Test {tid}:  {lbl}")
    print(f"          Raw p = {raw_p}  |  [NO BONFERRONI ADJUSTMENT — Exploratory]")
    if note:
        print(f"          Note: {note}")
    print(f"          Source: {script}")
print()

print(SEP)
print("  MANUSCRIPT MACRO VALUES  (copy into \\newcommand block)")
print(SEP)
print()
print(f"  %% Bonferroni family size")
print(f"  \\newcommand{{\\BonfK}}{{{K}}}  % k = confirmatory tests (Tests {', '.join(t[0] for t in confirmatory)})")
print()
for tid, lbl, raw_p, p_bonf, p_holm, survives, macro, _ in results:
    # Bonferroni macro
    p_bonf_rounded = round(p_bonf, 4)
    p_bonf_3 = f"{p_bonf_rounded:.3f}".rstrip('0').rstrip('.')
    if '.' not in p_bonf_3:
        p_bonf_3 += '.0'
    # Holm macro (3 sig figs for display)
    p_holm_rounded = round(p_holm, 4)
    p_holm_3 = f"{p_holm_rounded:.3f}".rstrip('0').rstrip('.')
    if '.' not in p_holm_3:
        p_holm_3 += '.0'
    # Survival decision uses Holm (primary), but T4 is labelled descriptive-only
    # because it is partially dependent on regional composition (see manuscript).
    # We keep its Holm-adj for transparency but flag it in the comment.
    descriptive_only = (tid == "4")
    status = ("SURVIVES (descriptive-only; see §temporal gradient)"
              if descriptive_only and survives
              else "SURVIVES" if survives else "ns")
    print(f"  % Test {tid} — {lbl}")
    print(f"  %   raw p = {raw_p:.6f}")
    print(f"  %   Bonferroni: {raw_p:.6f} × {K} = {p_bonf_rounded:.4f}")
    print(f"  %   Holm (step-down): adj p = {p_holm_rounded:.4f}  [{status}]")
    print(f"  \\newcommand{{\\pAdjTest{macro}}}{{{p_bonf_3}}}  % Bonferroni-adj p, Test {tid} ({lbl})")
    print(f"  \\newcommand{{\\pHolmTest{macro}}}{{{p_holm_3}}}  % Holm-adj p, Test {tid} ({lbl})")
    print()

print(SEP)
print()

# ── Supplementary Bonferroni-adjusted values (single tests, not in family) ───
# \pAdjClusterMW = cluster Mann-Whitney p × k (same k as main family)
_mw_raw = store.read("clusterMWp", default=None)
if _mw_raw is not None:
    _adj_mw = min(float(_mw_raw) * K, 1.0)
    print(f"  % Cluster asymmetry Mann-Whitney (Bonferroni-adj, k={K})")
    print(f"  \\newcommand{{\\pAdjClusterMW}}{{{_adj_mw:.3f}}}  % Bonf-adj p, cluster MW (= {_mw_raw:.4f} × {K})")
    print()
