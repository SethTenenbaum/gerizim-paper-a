#!/usr/bin/env python3
r"""Emit LaTeX macros for the full-corpus A+ Clopper-Pearson confidence interval.

This script reads the analysis results store (or falls back to the generated_macros.tex
file), computes the 95% Clopper-Pearson interval for \NclusterAp/\NclusterTotal, and
emits corresponding \newcommand macros.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

try:
    from scipy.stats import beta
except ImportError as exc:
    raise SystemExit("scipy is required to run this script. Install it with pip.") from exc

from lib.results_store import ResultsStore


MACROS_PATH = ROOT / "manuscript" / "generated_macros.tex"
STORE_PATH = ROOT / "data" / "store" / "results.json"


def read_counts_from_store() -> tuple[int, int] | None:
    store = ResultsStore()
    if not STORE_PATH.exists():
        return None
    try:
        n_ap = store.read("NclusterAp")
        n_tot = store.read("NclusterTotal")
    except KeyError:
        return None
    return int(n_ap), int(n_tot)


def read_counts_from_macros() -> tuple[int, int]:
    text = MACROS_PATH.read_text()
    m_ap = re.search(r"\\newcommand\{\\NclusterAp\}\{(\d+)\}", text)
    m_tot = re.search(r"\\newcommand\{\\NclusterTotal\}\{(\d+)\}", text)
    if not (m_ap and m_tot):
        raise SystemExit("Could not find \\NclusterAp or \\NclusterTotal in generated_macros.tex")
    return int(m_ap.group(1)), int(m_tot.group(1))


def clopper_pearson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if k == 0:
        lo = 0.0
    else:
        lo = beta.ppf(alpha / 2, k, n - k + 1)
    if k == n:
        hi = 1.0
    else:
        hi = beta.ppf(1 - alpha / 2, k + 1, n - k)
    return lo, hi


def format_pct(x: float) -> str:
    return f"{x * 100:.1f}"


def main() -> None:
    counts = read_counts_from_store() or read_counts_from_macros()
    n_ap, n_tot = counts
    lo, hi = clopper_pearson_ci(n_ap, n_tot)
    lo_pct = format_pct(lo)
    hi_pct = format_pct(hi)

    # If the store is available, write the macros there as well.
    if STORE_PATH.exists():
        store = ResultsStore()
        store.write_many({
            "clusterApCIlo": lo_pct,
            "clusterApCIhi": hi_pct,
        })

    print(f"  \\newcommand{{\\clusterApCIlo}}{{{lo_pct}}}  % Clopper-Pearson 95\\% CI lower, A+ rate (%)")
    print(f"  \\newcommand{{\\clusterApCIhi}}{{{hi_pct}}}  % Clopper-Pearson 95\\% CI upper, A+ rate (%)")


if __name__ == "__main__":
    main()
