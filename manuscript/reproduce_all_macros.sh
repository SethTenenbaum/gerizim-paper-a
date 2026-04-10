#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# reproduce_all_macros.sh
# ══════════════════════════════════════════════════════════════════════════════
#
# Runs every analysis script that produces numeric macros for
# paper_a_primary_unesco.tex, in dependency order.
#
# DESIGN
# ──────
# 1. The results store (data/store/results.json) is CLEARED at the start.
# 2. Every analysis script is run in dependency order.
#    Each script writes its p-values and statistics to the store via
#    lib/results_store.py — NO manual copying of values between scripts.
# 3. Summary scripts (bonferroni_correction.py, fdr_multiple_comparisons.py)
#    run LAST and read all values from the store.
#    Adding a test or changing a keyword only requires re-running this script.
#
# USAGE:
#   cd /path/to/gerizim-paper-a
#   bash manuscript/reproduce_all_macros.sh
#   bash manuscript/reproduce_all_macros.sh --macros-only   # LaTeX output only
#
# PREREQUISITES:
#   - Python 3.10+ with numpy, scipy, pandas
#   - pip install -r requirements.txt
#   - UNESCO XML:  data/store/unesco/unesco.xml
#   - P1435 CSV:   data/store/wikidata/p1435_global_control.csv
#
# ══════════════════════════════════════════════════════════════════════════════

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

echo "═══════════════════════════════════════════════════════════════════"
echo "  REPRODUCE ALL MACROS — Paper A (Primary UNESCO Analysis)"
echo "  Repository: $REPO_ROOT"
echo "  Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "═══════════════════════════════════════════════════════════════════"
echo ""

# ── Step 0: Clear the results store ──────────────────────────────────────────
# This ensures every value in the store was computed THIS run,
# not left over from a previous run with different code or data.
echo "Clearing results store..."
python3 - <<'PYEOF'
import sys
from pathlib import Path
sys.path.insert(0, str(Path(".").resolve()))
from lib.results_store import ResultsStore
ResultsStore().clear()
print("  ✓ data/store/results.json cleared")
PYEOF
echo ""

# ── --macros-only: run all scripts, collect only \newcommand lines ────────────
if [ "${1:-}" = "--macros-only" ]; then
    OUT_FILE="manuscript/generated_macros.tex"
    echo "% Generated macros — $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$OUT_FILE"
    echo "% Source: bash manuscript/reproduce_all_macros.sh --macros-only" >> "$OUT_FILE"
    echo "% IMPORTANT: values are computed fresh each run — do NOT hardcode" >> "$OUT_FILE"
    echo "" >> "$OUT_FILE"

    # DEPENDENCY ORDER:
    #   Primary analysis scripts first (they write to the store)
    #   Summary scripts last (they read from the store)
    SCRIPTS=(
        # ── Constants (config.json + corpus size) — no corpus analysis needed ──
        analysis/global/emit_constants.py                    # → GerizimLon, NwhcTotal, etc.
        analysis/global/geodesic_sensitivity.py              # → GeodesicApCurrent, GeodesicDropOut, GeodesicGainIn, etc.
        # ── Primary analysis (write to store) ─────────────────────────────
        analysis/unesco/spherical_monument_raw_sweep.py      # → pCircAp, pCircA, pCircChi
        analysis/unesco/spherical_monument_test.py            # → pCircAp_validated (context-validated, Exploratory 2x)
        analysis/unesco/cluster_asymmetry_test.py            # → clusterApBinom, clusterPermP
        analysis/unesco/harmonic_density_attractor_test.py
        analysis/unesco/unesco_buddhist_heritage_test.py     # → pBudAp, pBudA
        analysis/unesco/origin_sites_test.py                 # → pCanon, pPreTwoK, pModern
        analysis/unesco/founding_sites_analysis.py           # → pFoundKwFisher
        analysis/unesco/sacred_origin_test.py
        analysis/unesco/meta_keyword_test.py
        analysis/unesco/temporal_gradient_test.py            # → pCochranThree, pCochranFive
        analysis/unesco/deep_temporal_analysis.py            # → peakSigP, pFoundDateSpearman
        analysis/global/anchor_uniqueness_audit.py
        analysis/unesco/verify_x18_periodicity.py
        analysis/global/peak_geography_audit.py
        analysis/global/landmark_anchor_ranking.py
        analysis/global/global_corridor_comparison.py
        analysis/global/corridor_precision_test.py           # → pCorridorBinom, pCorridorFisher
        analysis/global/anchor_site_comparison.py
        analysis/unesco/spatial_independence_test.py         # → pNeffQuarter, pNeffHalf
        analysis/unesco/regional_temporal_gradient.py        # → pCMH
        analysis/unesco/tumulus_dome_evolution_raw_sweep.py  # → pEvoAp (Test 2b)
        analysis/unesco/tumulus_dome_evolution_test.py        # → pEvoAp_validated (context-validated, Exploratory 2bx)
        analysis/unesco/leave_one_out_sensitivity.py         # → LOOdomeN/Ap/Rate/Enrich/P, LOOstupaN/Ap/Rate/P
        analysis/global/wikidata_p1435_control_analysis.py
        analysis/global/dome_periodicity_audit.py
        analysis/americas/americas_harmonic_depletion_audit.py
        # ── Simulation (slow — writes simDomePermP etc.) ──────────────────
        analysis/unesco/simulation_null_model.py
        # ── Summary scripts (read from store) ─────────────────────────────
        analysis/unesco/fdr_multiple_comparisons.py          # reads all keys from config+store
        analysis/unesco/bonferroni_correction.py             # reads confirmatory keys from store
    )

    for s in "${SCRIPTS[@]}"; do
        # Skip comment lines
        [[ "$s" == \#* ]] && continue
        echo "Running $s ..."
        echo "% --- $s" >> "$OUT_FILE"
        if [ -f "$s" ]; then
            python3 "$s" 2>&1 | grep '\\newcommand' >> "$OUT_FILE" || true
        else
            echo "% SKIPPED — not found" >> "$OUT_FILE"
        fi
        echo "" >> "$OUT_FILE"
    done

    n_macros=$(grep -c '\\newcommand' "$OUT_FILE" || true)
    echo ""
    echo "═══════════════════════════════════════════════════════════════════"
    echo "  Wrote $n_macros \\newcommand lines to $OUT_FILE"
    echo "═══════════════════════════════════════════════════════════════════"
    exit 0
fi

# ── Data prerequisites ────────────────────────────────────────────────────────
echo "Checking data prerequisites..."
if [ ! -f "data/store/unesco/unesco.xml" ]; then
    echo "  ERROR: data/store/unesco/unesco.xml not found."
    echo "  Run: python3 data/scripts/fetch_extended.py"
    exit 1
fi
echo "  ✓ UNESCO XML found"

P1435_CSV="data/store/wikidata/p1435_global_control.csv"
HAS_P1435=true
if [ ! -f "$P1435_CSV" ]; then
    echo "  WARNING: $P1435_CSV not found."
    echo "  Wikidata P1435 and Americas scripts will be skipped."
    echo "  Run: python3 data/scripts/fetch_p1435_global.py"
    HAS_P1435=false
fi
echo ""

# ── Helper: run a script, print output ───────────────────────────────────────
run_script() {
    local label="$1"
    local script="$2"
    echo "─── $label ───"
    if [ -f "$script" ]; then
        python3 "$script" 2>&1 | tail -60
    else
        echo "  SKIPPED — $script not found"
    fi
    echo ""
}

# ════════════════════════════════════════════════════════════════════════════
# PRIMARY ANALYSIS SCRIPTS
# (Each writes its p-values and statistics to data/store/results.json)
# ════════════════════════════════════════════════════════════════════════════

run_script "GROUP 0: Constants (config.json + corpus size)"        \
           analysis/global/emit_constants.py

run_script "GROUP 0b: Geodesic latitude-sensitivity check (§3.1)"  \
           analysis/global/geodesic_sensitivity.py

run_script "GROUP 1: Dome/Spherical Monument Raw Sweep (Test 2)"  \
           analysis/unesco/spherical_monument_raw_sweep.py

run_script "GROUP 1b: Dome/Spherical Monument Context-Validated (Exploratory 2x)" \
           analysis/unesco/spherical_monument_test.py

run_script "GROUP 2: Cluster Asymmetry (Tests 1 & 3)"             \
           analysis/unesco/cluster_asymmetry_test.py

run_script "GROUP 3: Harmonic Density Attractor"                   \
           analysis/unesco/harmonic_density_attractor_test.py

run_script "GROUP 4: Buddhist Heritage"                            \
           analysis/unesco/unesco_buddhist_heritage_test.py

run_script "GROUP 5: Origin Sites (canon/pre-2000/modern)"        \
           analysis/unesco/origin_sites_test.py

run_script "GROUP 6: Founding Sites Analysis (keyword enrichment)" \
           analysis/unesco/founding_sites_analysis.py

run_script "GROUP 7: Sacred Origin Test"                           \
           analysis/unesco/sacred_origin_test.py

run_script "GROUP 8: Meta-Keyword Test"                            \
           analysis/unesco/meta_keyword_test.py

run_script "GROUP 9: Temporal Gradient (Test 4)"                   \
           analysis/unesco/temporal_gradient_test.py

run_script "GROUP 10: Deep Temporal Analysis"                      \
           analysis/unesco/deep_temporal_analysis.py

run_script "GROUP 11: Global Anchor Sweep"                         \
           analysis/global/anchor_uniqueness_audit.py

run_script "GROUP 12: x.18° Periodicity Verification"             \
           analysis/unesco/verify_x18_periodicity.py

run_script "GROUP 12b: Peak Geography Audit"                       \
           analysis/global/peak_geography_audit.py

run_script "GROUP 13: Landmark Anchor Ranking"                     \
           analysis/global/landmark_anchor_ranking.py

run_script "GROUP 14: Global Corridor Comparison"                  \
           analysis/global/global_corridor_comparison.py

run_script "GROUP 15: Corridor Precision Test"                     \
           analysis/global/corridor_precision_test.py

run_script "GROUP 16 + 22: Anchor Site Comparison"                 \
           analysis/global/anchor_site_comparison.py

run_script "GROUP 17: Spatial Independence / Autocorrelation"      \
           analysis/unesco/spatial_independence_test.py

run_script "GROUP 18: Regional Temporal Gradient (CMH)"            \
           analysis/unesco/regional_temporal_gradient.py

run_script "GROUP 20: Hemispherical Evolution Raw Sweep (Test 2b)" \
           analysis/unesco/tumulus_dome_evolution_raw_sweep.py

run_script "GROUP 20b: Hemispherical Evolution Context-Validated (Exploratory 2bx)" \
           analysis/unesco/tumulus_dome_evolution_test.py

run_script "GROUP 20c: Leave-One-Out Sensitivity (Tests 2 & 2b)" \
           analysis/unesco/leave_one_out_sensitivity.py

echo "─── GROUP 23: Wikidata P1435 Control Analysis ───"
if [ "$HAS_P1435" = true ]; then
    python3 analysis/global/wikidata_p1435_control_analysis.py 2>&1 | tail -40
else
    echo "  SKIPPED — P1435 CSV not found"
fi
echo ""

run_script "GROUP 24: Dome Periodicity Audit"                      \
           analysis/global/dome_periodicity_audit.py

echo "─── GROUP 25: Simulation Null Model (may take several minutes) ───"
python3 analysis/unesco/simulation_null_model.py 2>&1 | tail -40
echo ""

echo "─── GROUP 26: Americas Harmonic Depletion ───"
if [ "$HAS_P1435" = true ]; then
    python3 analysis/americas/americas_harmonic_depletion_audit.py 2>&1 | tail -30
else
    echo "  SKIPPED — P1435 CSV not found"
fi
echo ""

run_script "GROUP 27: Unit Sweep (spacing sensitivity)"            \
           analysis/unesco/unit_sweep_fill.py

# ════════════════════════════════════════════════════════════════════════════
# SUMMARY SCRIPTS — run AFTER all data-producing scripts have written the store
# ════════════════════════════════════════════════════════════════════════════

echo "═══════════════════════════════════════════════════════════════════"
echo "  SUMMARY SCRIPTS (read from data/store/results.json)"
echo "═══════════════════════════════════════════════════════════════════"
echo ""

run_script "GROUP 19: FDR Multiple Comparisons (reads store)"      \
           analysis/unesco/fdr_multiple_comparisons.py

run_script "GROUP 21: Bonferroni Correction (reads store)"         \
           analysis/unesco/bonferroni_correction.py

echo "═══════════════════════════════════════════════════════════════════"
echo "  ALL SCRIPTS COMPLETE"
echo "  Results store: data/store/results.json"
echo "  To regenerate LaTeX macros file:"
echo "    bash manuscript/reproduce_all_macros.sh --macros-only"
echo "═══════════════════════════════════════════════════════════════════"
