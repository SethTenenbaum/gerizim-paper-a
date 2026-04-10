"""
fdr_multiple_comparisons.py — FDR and family-wise corrections across ALL tests.

HOW P-VALUES ARE SOURCED
-------------------------
The test registry lives in config.json → "tests" → "fdr_all" → "tests".
Each entry maps a human label to a store key and a category (C/E).

Each analysis script writes its p-value to data/store/results.json.
This script reads from that store — NO hardcoded p-values here.

If a key is missing (i.e. the producing script hasn't been run yet),
that test is skipped and a warning is printed.  Run all analysis scripts
first via:  bash manuscript/reproduce_all_macros.sh

CORRECTIONS APPLIED
--------------------
  • Bonferroni (k = all tests present)
  • Holm step-down (k = all tests present)
  • Benjamini-Hochberg FDR
  • Confirmatory-subset BH FDR (C-labelled tests only)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).parent.parent.parent
sys.path.insert(0, str(_REPO))

from lib.results_store import ResultsStore

# ── Load config ───────────────────────────────────────────────────────────────
_CONFIG_PATH = _REPO / "config.json"
with open(_CONFIG_PATH) as _f:
    _CONFIG = json.load(_f)

# ── Correction helpers ────────────────────────────────────────────────────────

def benjamini_hochberg(pvals, alpha=0.05):
    m = len(pvals)
    if m == 0:
        return np.array([])
    sorted_idx = np.argsort(pvals)
    sorted_pvals = np.array(pvals)[sorted_idx]
    q_values = np.zeros(m)
    q_values[m - 1] = sorted_pvals[m - 1]
    for i in range(m - 2, -1, -1):
        q_values[i] = min(q_values[i + 1], sorted_pvals[i] * m / (i + 1))
    q_out = np.zeros(m)
    for i, idx in enumerate(sorted_idx):
        q_out[idx] = q_values[i]
    return q_out

def holm_correction(pvals):
    m = len(pvals)
    if m == 0:
        return np.array([])
    sorted_idx = np.argsort(pvals)
    sorted_pvals = np.array(pvals)[sorted_idx]
    adj = np.zeros(m)
    adj[0] = sorted_pvals[0] * m
    for i in range(1, m):
        adj[i] = max(adj[i - 1], sorted_pvals[i] * (m - i))
    adj = np.minimum(adj, 1.0)
    out = np.zeros(m)
    for i, idx in enumerate(sorted_idx):
        out[idx] = adj[i]
    return out

# ── Build test list from config + results store ───────────────────────────────

def main():
    store = ResultsStore()
    registry = _CONFIG["tests"]["fdr_all"]["tests"]

    tests   = []   # (label, p_value, category)
    skipped = []   # labels of tests whose key was not found in the store

    for entry in registry:
        label = entry["label"]
        key   = entry["p_key"]
        cat   = entry.get("cat", "C")
        val   = store.read(key, default=None)
        if val is None:
            skipped.append((label, key))
            continue
        tests.append((label, float(val), cat))

    if skipped:
        print(f"\n  WARNING: {len(skipped)} test(s) missing from results store "
              f"(run their producing scripts first):")
        for lbl, k in skipped:
            print(f"    • [{k}]  {lbl}")
        print()

    if not tests:
        print("\n  ERROR: results store is empty.  Run all analysis scripts first:")
        print("    bash manuscript/reproduce_all_macros.sh")
        sys.exit(1)

    names = [t[0] for t in tests]
    pvals = np.array([t[1] for t in tests])
    cats  = [t[2] for t in tests]
    m     = len(tests)

    # ── Corrections ───────────────────────────────────────────────────────────
    bonf  = np.minimum(pvals * m, 1.0)
    holm  = holm_correction(pvals)
    bh_q  = benjamini_hochberg(pvals)

    conf_idx = [i for i, c in enumerate(cats) if c == "C"]
    expl_idx = [i for i, c in enumerate(cats) if c == "E"]

    conf_pvals = pvals[conf_idx]
    bh_conf = benjamini_hochberg(conf_pvals) if len(conf_pvals) else np.array([])

    # ── Table ─────────────────────────────────────────────────────────────────
    print("=" * 110)
    print("  MULTIPLE COMPARISONS: FDR AND FAMILY-WISE CORRECTIONS")
    print(f"  P-values sourced from   : data/store/results.json")
    print(f"  Total tests present     : {m}  |  Confirmatory: {len(conf_idx)}  |  Exploratory: {len(expl_idx)}")
    if skipped:
        print(f"  Tests skipped (missing) : {len(skipped)}")
    print("=" * 110)

    header = f"  {'Test':<45} {'Cat':>3} {'p_raw':>9} {'Bonf':>9} {'Holm':>9} {'BH q':>9} {'Sig':>5}"
    print(f"\n{header}")
    print(f"  {'-'*45} {'-'*3:>3} {'-'*9:>9} {'-'*9:>9} {'-'*9:>9} {'-'*9:>9} {'-'*5:>5}")

    for i in range(m):
        sig_lbl = ("***" if bh_q[i] < 0.001 else "**" if bh_q[i] < 0.01
                   else "*" if bh_q[i] < 0.05 else "~" if bh_q[i] < 0.10 else "ns")
        print(f"  {names[i]:<45} {cats[i]:>3} {pvals[i]:9.4f} {bonf[i]:9.4f} "
              f"{holm[i]:9.4f} {bh_q[i]:9.4f} {sig_lbl:>5}")

    # ── Summary counts ────────────────────────────────────────────────────────
    print(f"\n" + "─" * 110)
    print(f"  SUMMARY: How many tests survive each correction at α = 0.05?")
    print(f"  {'Correction':<35} {'Significant':>12} {'of':>4} {m}")
    print(f"  {'-'*35} {'-'*12:>12} {'-'*4:>4}")
    print(f"  {'Raw (uncorrected)':<35} {sum(pvals < 0.05):>12}")
    print(f"  {'Bonferroni (k={})'.format(m):<35} {sum(bonf < 0.05):>12}")
    print(f"  {'Holm step-down (k={})'.format(m):<35} {sum(holm < 0.05):>12}")
    print(f"  {'Benjamini-Hochberg FDR (k={})'.format(m):<35} {sum(bh_q < 0.05):>12}")

    # ── Confirmatory subset ───────────────────────────────────────────────────
    if len(conf_idx) > 0:
        print(f"\n" + "─" * 110)
        print(f"  CONFIRMATORY TESTS ONLY (pre-specified, N={len(conf_idx)})")
        print(f"  {'Correction':<35} {'Significant':>12} {'of':>4} {len(conf_idx)}")
        print(f"  {'-'*35} {'-'*12:>12} {'-'*4:>4}")
        bonf_conf = np.minimum(conf_pvals * len(conf_idx), 1.0)
        holm_conf = holm_correction(conf_pvals)
        print(f"  {'Raw (uncorrected)':<35} {sum(conf_pvals < 0.05):>12}")
        print(f"  {'Bonferroni (k={})'.format(len(conf_idx)):<35} {sum(bonf_conf < 0.05):>12}")
        print(f"  {'Holm step-down (k={})'.format(len(conf_idx)):<35} {sum(holm_conf < 0.05):>12}")
        print(f"  {'Benjamini-Hochberg FDR (k={})'.format(len(conf_idx)):<35} {sum(bh_conf < 0.05):>12}")

    # ── Paper's pre-specified Bonferroni (k=4, from config) ──────────────────
    bonf_family = _CONFIG["tests"]["confirmatory"]
    family_alpha = _CONFIG["tests"]["family_alpha"]
    k_paper = len(bonf_family)
    print(f"\n" + "─" * 110)
    print(f"  PAPER'S BONFERRONI (k={k_paper}, family_alpha={family_alpha}) — READ FROM CONFIG + STORE")
    for e in bonf_family:
        raw_p = store.read(e["p_key"], default=None)
        if raw_p is None:
            print(f"  Test {e['id']}: MISSING in store — run {e['script']}")
            continue
        p_adj = min(float(raw_p) * k_paper, 1.0)
        star = "*" if p_adj < family_alpha else "ns"
        print(f"  Test {e['id']:>3}  {e['label']:<50}  "
              f"p_raw={float(raw_p):.4f}  p_adj={p_adj:.4f}  {star}")

    # ── Interpretation ────────────────────────────────────────────────────────
    n_bh_sig   = int(sum(bh_q < 0.05))
    n_raw_sig  = int(sum(pvals < 0.05))
    print(f"\n" + "=" * 110)
    print(f"""
  {n_raw_sig} tests are significant at α = 0.05 (uncorrected).
  {n_bh_sig} survive Benjamini-Hochberg FDR correction across all {m} tests.

  Key results surviving full FDR correction (q < 0.05):""")
    for i in range(m):
        if bh_q[i] < 0.05:
            print(f"    • {names[i]:<45s}  p = {pvals[i]:.4f}  q = {bh_q[i]:.4f}")

    # ── LaTeX macros (GROUP 19) ───────────────────────────────────────────────
    n_conf        = len(conf_idx)
    n_expl        = m - n_conf
    n_survive_fdr  = int(np.sum(bh_q < 0.05))
    n_survive_bonf = int(np.sum(bonf < 0.05))

    print()
    print("  % LaTeX macros (GROUP 19):")
    print(f"  \\newcommand{{\\NtotalTests}}{{{m}}}           % total tests in FDR analysis")
    print(f"  \\newcommand{{\\NconfTests}}{{{n_conf}}}           % confirmatory tests")
    print(f"  \\newcommand{{\\NexplTests}}{{{n_expl}}}           % exploratory tests")
    print(f"  \\newcommand{{\\NsurviveFDR}}{{{n_survive_fdr}}}           % tests surviving BH FDR (q<0.05)")
    print(f"  \\newcommand{{\\NsurviveBonfAll}}{{{n_survive_bonf}}}           % tests surviving Bonferroni (all tests)")

    # ── Write summary counts to store ─────────────────────────────────────────
    ResultsStore().write_many({
        "NtotalTests":     m,
        "NconfTests":      n_conf,
        "NsurviveFDR":     n_survive_fdr,
        "NsurviveBonfAll": n_survive_bonf,
    })
    print()
    return 0

if __name__ == "__main__":
    main()
