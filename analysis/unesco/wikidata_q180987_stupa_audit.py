"""
wikidata_q180987_stupa_audit.py
================================
Beru-harmonic audit for the Wikidata Q180987 stupa corpus — a pipeline
portability demonstration on a corpus entirely outside the UNESCO WHC XML.

Q180987 is the Wikidata class identifier for stupa (a hemispherical
Buddhist/Jain monument).

Corpus  : data/store/unesco/wikidata_stupas_q180987.csv
Fetch   : data/scripts/fetch_wikidata_q180987.py

TESTS
-----
1. A-tier enrichment : one-sided exact binomial (H1: rate > 20% null)
2. A+ enrichment     : one-sided exact binomial (H1: rate > 4% null)
3. Circular-shift permutation (100,000 trials, seed 42):
   shift all site longitudes by the same uniform-random offset in
   [0, 3°) — one full harmonic period — and count A-tier hits.
   Z and one-sided p vs observed.
4. Cluster asymmetry : mean site count at nodes with ≥1 A-tier site
   versus nodes with no A-tier site.
5. Java-region spotlight: sites within 107–115°E (2.5-beru node area).
6. Named-site deviations: Borobudur, Candi Ngawen, Lumbini.

USAGE
-----
    cd /path/to/gerizim-paper-a
    python3 analysis/unesco/wikidata_q180987_stupa_audit.py

MANUSCRIPT MACROS PRODUCED
---------------------------
    \\wikiStupaTotal          — N = 229
    \\wikiStupaATierCount     — Tier-A (A++ + A+ + A) count
    \\wikiStupaATierRate      — Tier-A rate (%)
    \\wikiStupaATierEnrich    — enrichment factor vs 20% null
    \\wikiStupaATierBinomP    — binomial p-value for A-tier enrichment
    \\wikiStupaATierBinomStar — significance label for A-tier test
    \\wikiStupaPermNperms     — number of permutation trials
    \\wikiStupaPermZ          — circular-shift permutation Z-score
    \\wikiStupaPermP          — circular-shift permutation p-value
    \\wikiStupaPermStar       — significance label for permutation test
    \\wikiStupaApCount        — combined A+/A++ count
    \\wikiStupaApRate         — combined A+/A++ rate (%)
    \\wikiStupaApEnrich       — A+ enrichment vs 4% null
    \\wikiStupaApBinomP       — binomial p-value for A+ enrichment
    \\wikiStupaApBinomStar    — significance label for A+ test
    \\wikiStupaJavaTotal      — sites in Java-region window (107–115°E)
    \\wikiStupaJavaATierCount — A-tier sites in Java region
    \\wikiStupaJavaATierRate  — A-tier rate in Java region (%)
    \\wikiStupaClusterRatio   — mean(A-node count) / mean(non-A-node count)
    \\wikiStupaClusterAMean   — mean site count at A-tier-bearing nodes
    \\wikiStupaClusterNonAMean — mean site count at non-A-tier nodes
    \\wikiStupaBorobudurDevKm — Borobudur deviation in km
    \\wikiStupaNgawenDevKm    — Candi Ngawen deviation in km (in metres: ×1000)
    \\wikiStupaNgawenDevM     — Candi Ngawen deviation in metres
    \\wikiStupaLumbiniDevKm   — Lumbini deviation in km
    \\wikiStupaLumbiniDevM    — Lumbini deviation in metres
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, fisher_exact

np.random.seed(42)

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from lib.beru import (
    GERIZIM, BERU,
    TIER_A_MAX, TIER_APLUS,
    P_NULL_AP, P_NULL_A,
    deviation, tier_label, full_calculation,
)
from lib.results_store import ResultsStore

CORPUS_CSV = REPO_ROOT / "data" / "store" / "unesco" / "wikidata_stupas_q180987.csv"

# Java-region longitude window (2.5-beru harmonic area)
JAVA_LON_MIN = 107.0
JAVA_LON_MAX = 115.0

# Permutation parameters
N_PERMS = 100_000
HARMONIC_PERIOD_DEG = BERU * 0.1   # 3.0°  — one full grid period

# Named sites to spotlight
NAMED_SITES = {
    "borobudur": ("Borobudur", "wikiStupaBorobudurDevKm"),
    "ngawen":    ("Ngawen",    "wikiStupaNgawenDevKm"),
    "lumbini":   ("Lumbini",   "wikiStupaLumbiniDevKm"),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sig_label(p: float) -> str:
    """Return a LaTeX-safe significance label."""
    if p < 0.001:  return "***"
    if p < 0.01:   return "**"
    if p < 0.05:   return "*"
    if p < 0.10:   return r"\ensuremath{\sim}"
    return "ns"


def fmt_p(p: float) -> str:
    """Format p-value for LaTeX macro (no math-mode symbols)."""
    if p < 0.001:
        return "< 0.001"
    return f"{p:.3f}"


# ---------------------------------------------------------------------------
# Load corpus
# ---------------------------------------------------------------------------

def load_corpus() -> list[dict]:
    sites = []
    with open(CORPUS_CSV, newline="") as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    with open(CORPUS_CSV, newline="") as _f:
        pass  # just to satisfy the with-block; we'll use lines directly

    reader = csv.DictReader(lines)
    for row in reader:
        try:
            lat = float(row["lat"])
            lon = float(row["lon"])
        except (ValueError, KeyError):
            continue
        r = full_calculation(lon)
        sites.append({
            "qid":    row.get("qid", ""),
            "name":   row.get("name", ""),
            "lat":    lat,
            "lon":    lon,
            "dev":    r["dev"],
            "dev_km": r["dev_km"],
            "tier":   r["tier"],
        })
    return sites


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    if not CORPUS_CSV.exists():
        print(f"ERROR: corpus not found at {CORPUS_CSV}", file=sys.stderr)
        print("Run:  python3 data/scripts/fetch_wikidata_q180987.py", file=sys.stderr)
        sys.exit(1)

    sites = load_corpus()
    N = len(sites)
    lons = [s["lon"] for s in sites]

    # ── Tier counts ──────────────────────────────────────────────────────────
    is_a_tier   = lambda s: s["tier"] in ("A++", "A+", "A")
    is_ap_tier  = lambda s: s["tier"] in ("A++", "A+")

    n_a_tier  = sum(1 for s in sites if is_a_tier(s))
    n_ap_tier = sum(1 for s in sites if is_ap_tier(s))
    rate_a    = 100.0 * n_a_tier  / N
    rate_ap   = 100.0 * n_ap_tier / N
    enrich_a  = (n_a_tier  / N) / P_NULL_A   # vs 20% null
    enrich_ap = (n_ap_tier / N) / P_NULL_AP  # vs 4% null

    binom_a  = binomtest(n_a_tier,  N, P_NULL_A,  alternative="greater")
    binom_ap = binomtest(n_ap_tier, N, P_NULL_AP, alternative="greater")

    # ── Circular-shift permutation ───────────────────────────────────────────
    null_counts = np.empty(N_PERMS, dtype=int)
    shifts = np.random.uniform(0.0, HARMONIC_PERIOD_DEG, N_PERMS)
    for i, shift in enumerate(shifts):
        null_counts[i] = sum(
            1 for lon in lons
            if deviation(lon + shift) <= TIER_A_MAX
        )
    null_mean = null_counts.mean()
    null_std  = null_counts.std(ddof=1)
    z_perm    = (n_a_tier - null_mean) / null_std if null_std > 0 else 0.0
    p_perm    = float((null_counts >= n_a_tier).mean())

    # ── Java-region spotlight ────────────────────────────────────────────────
    java_sites  = [s for s in sites if JAVA_LON_MIN <= s["lon"] <= JAVA_LON_MAX]
    n_java      = len(java_sites)
    n_java_a    = sum(1 for s in java_sites if is_a_tier(s))
    rate_java_a = 100.0 * n_java_a / n_java if n_java > 0 else 0.0

    # ── Geographic breakdown (for paper macros) ──────────────────────────────
    is_c_band   = lambda s: s["tier"] in ("C", "C-", "C--")

    K_AP = sum(1 for s in sites if is_ap_tier(s))
    K_A  = sum(1 for s in sites if is_a_tier(s))
    K_C  = sum(1 for s in sites if is_c_band(s))

    def _fisher(sub, k_sub, K_bg, alternative="greater"):
        n = len(sub)
        if n == 0:
            return float("nan"), float("nan")
        table = [[k_sub, n - k_sub], [K_bg - k_sub, N - n - (K_bg - k_sub)]]
        or_, p = fisher_exact(table, alternative=alternative)
        return round(or_, 2), round(p, 4)

    # Java/Sumatra node 105–115°E
    java115 = [s for s in sites if 105.0 <= s["lon"] <= 115.0]
    n_java115    = len(java115)
    nap_java115  = sum(1 for s in java115 if is_ap_tier(s))
    na_java115   = sum(1 for s in java115 if is_a_tier(s))
    or_a_java115,  p_a_java115  = _fisher(java115, na_java115,  K_A)
    or_ap_java115, p_ap_java115 = _fisher(java115, nap_java115, K_AP)

    # Java tight 107–112°E
    java112 = [s for s in sites if 107.0 <= s["lon"] <= 112.0]
    n_java112    = len(java112)
    nap_java112  = sum(1 for s in java112 if is_ap_tier(s))
    na_java112   = sum(1 for s in java112 if is_a_tier(s))
    or_a_java112,  p_a_java112  = _fisher(java112, na_java112,  K_A)
    or_ap_java112, p_ap_java112 = _fisher(java112, nap_java112, K_AP)

    # Myanmar/Cambodia 90–105°E
    myanmar = [s for s in sites if 90.0 <= s["lon"] <= 105.0]
    n_myanmar  = len(myanmar)
    nc_myanmar = sum(1 for s in myanmar if is_c_band(s))
    or_c_myanmar, p_c_myanmar = _fisher(myanmar, nc_myanmar, K_C)

    # Heartland 70–105°E combined
    heartland70 = [s for s in sites if 70.0 <= s["lon"] <= 105.0]
    n_hl70  = len(heartland70)
    nc_hl70 = sum(1 for s in heartland70 if is_c_band(s))
    or_c_hl70, p_c_hl70 = _fisher(heartland70, nc_hl70, K_C)

    # ── Cluster asymmetry ────────────────────────────────────────────────────
    # Group sites by nearest harmonic node (0.1-beru step)
    node_sites: dict[float, list] = defaultdict(list)
    for s in sites:
        arc      = abs(s["lon"] - GERIZIM)
        beru_val = arc / BERU
        nearest  = round(beru_val / 0.1) * 0.1
        node_sites[round(nearest, 2)].append(s)

    a_node_counts     = []
    non_a_node_counts = []
    for _node, ss in node_sites.items():
        count = len(ss)
        if any(is_a_tier(s) for s in ss):
            a_node_counts.append(count)
        else:
            non_a_node_counts.append(count)

    mean_a_nodes     = float(np.mean(a_node_counts))     if a_node_counts     else 0.0
    mean_non_a_nodes = float(np.mean(non_a_node_counts)) if non_a_node_counts else 0.0
    cluster_ratio    = (mean_a_nodes / mean_non_a_nodes
                        if mean_non_a_nodes > 0 else float("inf"))

    # ── Named-site deviations ────────────────────────────────────────────────
    named_devs: dict[str, float] = {}
    for key, (fragment, _macro) in NAMED_SITES.items():
        matches = [s for s in sites if fragment.lower() in s["name"].lower()]
        if matches:
            # Pick the one with smallest deviation
            best = min(matches, key=lambda s: s["dev"])
            named_devs[key] = best["dev_km"]
        else:
            named_devs[key] = float("nan")

    # ── Print report ─────────────────────────────────────────────────────────
    sep  = "=" * 78
    dash = "─" * 78
    print(sep)
    print("  WIKIDATA Q180987 STUPA CORPUS — BERU HARMONIC AUDIT")
    print(f"  Q180987: Wikidata class identifier for stupa structures")
    print(sep)
    print(f"  Corpus  : {CORPUS_CSV.relative_to(REPO_ROOT)}")
    print(f"  Anchor  : Mount Gerizim  {GERIZIM}°E")
    print(f"  Beru    : {BERU}°  |  harmonic step = 0.1 beru = 3.0°")
    print(f"  N       : {N}")
    print()

    print(dash)
    print("  TIER DISTRIBUTION")
    print(dash)
    for tier in ("A++", "A+", "A", "B", "C"):
        cnt = sum(1 for s in sites if s["tier"] == tier)
        print(f"    {tier:<4}  {cnt:>3}  ({100*cnt/N:.1f}%)")
    print()

    print(dash)
    print("  A-TIER ENRICHMENT  (A++ + A+ + A  vs 20% geometric null)")
    print(dash)
    print(f"  Observed  : {n_a_tier}/{N}  = {rate_a:.1f}%")
    print(f"  Null      : {100*P_NULL_A:.0f}%")
    print(f"  Enrichment: {enrich_a:.2f}×")
    print(f"  Binomial p: {binom_a.pvalue:.4f}  {sig_label(binom_a.pvalue)}")
    print()

    print(dash)
    print("  CIRCULAR-SHIFT PERMUTATION TEST  "
          f"(N={N_PERMS:,}, seed=42, shift∈[0,{HARMONIC_PERIOD_DEG}°))")
    print(dash)
    print(f"  Observed A-tier   : {n_a_tier}")
    print(f"  Null mean ± SD    : {null_mean:.2f} ± {null_std:.2f}")
    print(f"  Z                 : {z_perm:.3f}")
    print(f"  p (one-sided, ≥obs): {p_perm:.4f}  {sig_label(p_perm)}")
    print()

    print(dash)
    print("  A+ ENRICHMENT  (A++ + A+  vs 4% geometric null)")
    print(dash)
    print(f"  Observed  : {n_ap_tier}/{N}  = {rate_ap:.1f}%")
    print(f"  Null      : {100*P_NULL_AP:.0f}%")
    print(f"  Enrichment: {enrich_ap:.2f}×")
    print(f"  Binomial p: {binom_ap.pvalue:.4f}  {sig_label(binom_ap.pvalue)}")
    print()

    print(dash)
    print(f"  JAVA-REGION SPOTLIGHT  ({JAVA_LON_MIN}–{JAVA_LON_MAX}°E, "
          f"near 2.5-beru node = {GERIZIM + 2.5*BERU:.3f}°E)")
    print(dash)
    print(f"  N in region       : {n_java}")
    print(f"  A-tier in region  : {n_java_a}  ({rate_java_a:.1f}%)")
    print()

    print(dash)
    print("  CLUSTER ASYMMETRY  (sites per harmonic node)")
    print(dash)
    print(f"  Nodes with ≥1 A-tier site : {len(a_node_counts)}  "
          f"(mean {mean_a_nodes:.2f} sites/node)")
    print(f"  Nodes with no A-tier site : {len(non_a_node_counts)}  "
          f"(mean {mean_non_a_nodes:.2f} sites/node)")
    print(f"  Ratio                     : {cluster_ratio:.2f}×")
    print()

    print(dash)
    print("  NAMED-SITE DEVIATIONS")
    print(dash)
    for key, (fragment, _) in NAMED_SITES.items():
        km = named_devs.get(key, float("nan"))
        m  = km * 1000
        print(f"  {fragment:<30}  dev = {km:.3f} km  = {m:.0f} m  "
              f"tier={([s['tier'] for s in sites if fragment.lower() in s['name'].lower()] or ['?'])[0]}")
    print()

    # ── LaTeX macros ─────────────────────────────────────────────────────────
    print(sep)
    print("  LATEX MACROS")
    print(sep)

    macros = [
        ("wikiStupaTotal",           str(N)),
        ("wikiStupaATierCount",      str(n_a_tier)),
        ("wikiStupaATierRate",       f"{rate_a:.1f}"),
        ("wikiStupaATierEnrich",     f"{enrich_a:.2f}"),
        ("wikiStupaATierBinomP",     fmt_p(binom_a.pvalue)),
        ("wikiStupaATierBinomStar",  sig_label(binom_a.pvalue)),
        ("wikiStupaPermNperms",      f"{N_PERMS:,}"),
        ("wikiStupaPermZ",           f"{z_perm:.2f}"),
        ("wikiStupaPermP",           f"{p_perm:.3f}"),
        ("wikiStupaPermStar",        sig_label(p_perm)),
        ("wikiStupaApCount",         str(n_ap_tier)),
        ("wikiStupaApRate",          f"{rate_ap:.1f}"),
        ("wikiStupaApEnrich",        f"{enrich_ap:.2f}"),
        ("wikiStupaApBinomP",        fmt_p(binom_ap.pvalue)),
        ("wikiStupaApBinomStar",     sig_label(binom_ap.pvalue)),
        ("wikiStupaJavaTotal",       str(n_java)),
        ("wikiStupaJavaATierCount",  str(n_java_a)),
        ("wikiStupaJavaATierRate",   f"{rate_java_a:.1f}"),
        ("wikiStupaClusterRatio",    f"{cluster_ratio:.2f}"),
        ("wikiStupaClusterAMean",    f"{mean_a_nodes:.1f}"),
        ("wikiStupaClusterNonAMean", f"{mean_non_a_nodes:.1f}"),
        ("wikiStupaBorobudurDevKm",  f"{named_devs.get('borobudur', 0):.2f}"),
        ("wikiStupaNgawenDevKm",     f"{named_devs.get('ngawen', 0):.3f}"),
        ("wikiStupaNgawenDevM",      f"{named_devs.get('ngawen', 0)*1000:.0f}"),
        ("wikiStupaLumbiniDevKm",    f"{named_devs.get('lumbini', 0):.3f}"),
        ("wikiStupaLumbiniDevM",     f"{named_devs.get('lumbini', 0)*1000:.0f}"),
        # Geographic breakdown macros
        ("wikiJavaNodeN",            str(n_java115)),
        ("wikiJavaNodeATierN",       str(na_java115)),
        ("wikiJavaNodeATierRate",    f"{100*na_java115/n_java115:.1f}" if n_java115 else "0"),
        ("wikiJavaNodeATierOR",      f"{or_a_java115:.2f}"),
        ("wikiJavaNodeATierP",       f"{p_a_java115:.3f}"),
        ("wikiJavaTightN",           str(n_java112)),
        ("wikiJavaTightATierN",      str(na_java112)),
        ("wikiJavaTightATierRate",   f"{100*na_java112/n_java112:.1f}" if n_java112 else "0"),
        ("wikiJavaTightATierOR",     f"{or_a_java112:.2f}"),
        ("wikiJavaTightATierP",      f"{p_a_java112:.3f}"),
        ("wikiJavaTightApN",         str(nap_java112)),
        ("wikiJavaTightApRate",      f"{100*nap_java112/n_java112:.1f}" if n_java112 else "0"),
        ("wikiJavaTightApOR",        f"{or_ap_java112:.2f}"),
        ("wikiJavaTightApP",         f"{p_ap_java112:.3f}"),
        ("wikiMyanmarN",             str(n_myanmar)),
        ("wikiMyanmarCbandN",        str(nc_myanmar)),
        ("wikiMyanmarCbandRate",     f"{100*nc_myanmar/n_myanmar:.1f}" if n_myanmar else "0"),
        ("wikiMyanmarCbandOR",       f"{or_c_myanmar:.2f}"),
        ("wikiMyanmarCbandP",        f"{p_c_myanmar:.3f}"),
        ("wikiHeartlandN",           str(n_hl70)),
        ("wikiHeartlandCbandN",      str(nc_hl70)),
        ("wikiHeartlandCbandRate",   f"{100*nc_hl70/n_hl70:.1f}" if n_hl70 else "0"),
        ("wikiHeartlandCbandOR",     f"{or_c_hl70:.2f}"),
        ("wikiHeartlandCbandP",      f"{p_c_hl70:.3f}"),
    ]

    for name, val in macros:
        print(f"\\newcommand{{\\{name}}}{{{val}}}")

    # ── Results store ─────────────────────────────────────────────────────────
    store_data = {k: v for k, v in [
        ("wikiStupaTotal",           float(N)),
        ("wikiStupaATierCount",      float(n_a_tier)),
        ("wikiStupaATierRate",       round(rate_a, 1)),
        ("wikiStupaATierEnrich",     round(enrich_a, 2)),
        ("wikiStupaATierBinomP",     round(binom_a.pvalue, 6)),
        ("wikiStupaPermNperms",      float(N_PERMS)),
        ("wikiStupaPermZ",           round(z_perm, 2)),
        ("wikiStupaPermP",           round(p_perm, 3)),
        ("wikiStupaApCount",         float(n_ap_tier)),
        ("wikiStupaApRate",          round(rate_ap, 1)),
        ("wikiStupaApEnrich",        round(enrich_ap, 2)),
        ("wikiStupaApBinomP",        round(binom_ap.pvalue, 4)),
        ("wikiStupaJavaTotal",       float(n_java)),
        ("wikiStupaJavaATierCount",  float(n_java_a)),
        ("wikiStupaJavaATierRate",   round(rate_java_a, 1)),
        ("wikiStupaClusterRatio",    round(cluster_ratio, 2)),
        ("wikiStupaClusterAMean",    round(mean_a_nodes, 1)),
        ("wikiStupaClusterNonAMean", round(mean_non_a_nodes, 1)),
        ("wikiStupaBorobudurDevKm",  round(named_devs.get("borobudur", 0), 3)),
        ("wikiStupaNgawenDevKm",     round(named_devs.get("ngawen", 0), 4)),
        ("wikiStupaLumbiniDevKm",    round(named_devs.get("lumbini", 0), 4)),
        # Geographic breakdown
        ("wikiJavaNodeN",            float(n_java115)),
        ("wikiJavaNodeATierN",       float(na_java115)),
        ("wikiJavaNodeATierRate",    round(100*na_java115/n_java115, 1) if n_java115 else 0.0),
        ("wikiJavaNodeATierOR",      float(or_a_java115)),
        ("wikiJavaNodeATierP",       float(p_a_java115)),
        ("wikiJavaTightN",           float(n_java112)),
        ("wikiJavaTightATierN",      float(na_java112)),
        ("wikiJavaTightATierRate",   round(100*na_java112/n_java112, 1) if n_java112 else 0.0),
        ("wikiJavaTightATierOR",     float(or_a_java112)),
        ("wikiJavaTightATierP",      float(p_a_java112)),
        ("wikiJavaTightApN",         float(nap_java112)),
        ("wikiJavaTightApRate",      round(100*nap_java112/n_java112, 1) if n_java112 else 0.0),
        ("wikiJavaTightApOR",        float(or_ap_java112)),
        ("wikiJavaTightApP",         float(p_ap_java112)),
        ("wikiMyanmarN",             float(n_myanmar)),
        ("wikiMyanmarCbandN",        float(nc_myanmar)),
        ("wikiMyanmarCbandRate",     round(100*nc_myanmar/n_myanmar, 1) if n_myanmar else 0.0),
        ("wikiMyanmarCbandOR",       float(or_c_myanmar)),
        ("wikiMyanmarCbandP",        float(p_c_myanmar)),
        ("wikiHeartlandN",           float(n_hl70)),
        ("wikiHeartlandCbandN",      float(nc_hl70)),
        ("wikiHeartlandCbandRate",   round(100*nc_hl70/n_hl70, 1) if n_hl70 else 0.0),
        ("wikiHeartlandCbandOR",     float(or_c_hl70)),
        ("wikiHeartlandCbandP",      float(p_c_hl70)),
    ]}
    ResultsStore().write_many(store_data)
    print()
    print("Results written to data/store/results.json")


if __name__ == "__main__":
    main()
