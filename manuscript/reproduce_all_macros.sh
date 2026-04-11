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
        analysis/unesco/dome_leave_k_out.py                  # → LKOsiteOne/Two/ThreeName/DevKm, LKOtwo*/three*/worstTwo/ThreeP
        analysis/global/wikidata_p1435_control_analysis.py
        analysis/global/dome_periodicity_audit.py
        analysis/global/x18_periodicity_formal_test.py       # → rayleighR/Z/PermP, fullRayleighR/Z/PermP, maxApCount/Z/PermP (GROUP 11b)
        analysis/americas/americas_harmonic_depletion_audit.py
        # ── Reviewer robustness checks ────────────────────────────────────
        analysis/unesco/dome_footprint_window_sensitivity.py   # → geoNullDomeTwoDeg*, geoNullDomeTenDeg*
        analysis/unesco/stupa_coordinate_perturbation.py       # → stupaCoordPerturb*
        analysis/unesco/americas_directional_test.py           # → AmericasN, AmericasApCount, AmericasOneSidedP, AmericasDirectional
        # ── Simulation (slow — writes simDomePermP etc.) ──────────────────
        analysis/unesco/simulation_null_model.py
        # ── Cluster asymmetry + conditional sensitivity subtest ───────────
        # (cluster_asymmetry_test.py now includes the conditional subtest at the end;
        #  it emits clusterCondN, clusterCondClusterOk, clusterCondPeakFrac, etc.)

        analysis/unesco/unit_sweep_fill.py
        analysis/unesco/sensitivity_slope_permutation_test.py  # → permSlopeNperms, permSlopePcanon, permSlopePsharp
        analysis/unesco/sensitivity_slope_specificity_test.py  # → permSlopeCanonRankObs, permSlopeCanonBestJoint, permSlopeCanonRankP, permSlopeCanonCondFraction
        analysis/unesco/dome_geographic_concentration_test.py  # → geoNullDomeBootP/Z/Mean, geoNullDomeRestrictedP/Z/Mean/N, domeEurasianFraction
        analysis/unesco/stupa_geographic_concentration_test.py # → stupaGeoBootP/Z/Mean, stupaRegionP/Z/N, stupaGeoRestrictedP/Z/N (GROUP 29b)
        analysis/unesco/dome_founding_stratification.py        # → domeStratNfs, domeStratNnfs, domeStratF/NfsAp/Rate/Enrich/P, domeStratFisherOR/P
        # ── Geographic concentration robustness ───────────────────────────
        analysis/unesco/region_conditioned_permutation.py    # → regionCondPermP/Z/Obs/Mean/Std/Nregions/Nperms (GROUP 35)
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

    n_macros_raw=$(grep -c '\\newcommand' "$OUT_FILE" || true)

    # ── Deduplicate: keep LAST definition of each \newcommand ─────────────────
    # Multiple scripts may emit the same macro name (e.g. summary scripts
    # re-emit values already written by primary scripts).  LaTeX requires each
    # \newcommand to be unique.  We post-process the file: scan from the bottom,
    # keep only the first (= last in file order) occurrence of each macro name,
    # preserve all non-\newcommand lines (comments, blanks, section headers).
    python3 - "$OUT_FILE" <<'PYEOF'
import sys, re
from pathlib import Path

path = Path(sys.argv[1])
lines = path.read_text().splitlines(keepends=True)

# Regex to extract the macro name from a \newcommand line
_RE = re.compile(r'\\newcommand\{(\\[A-Za-z@]+)\}')

seen = set()
kept = []
# Scan in reverse so we always keep the LAST definition
for line in reversed(lines):
    m = _RE.search(line)
    if m:
        name = m.group(1)
        if name in seen:
            continue          # drop earlier (duplicate) definition
        seen.add(name)
    kept.append(line)

# Restore original order
kept.reverse()
path.write_text("".join(kept))
print(f"  ✓ Deduplicated: kept {len(seen)} unique macros "
      f"(removed {len(lines)-len(kept)} duplicate lines)")
PYEOF

    n_macros=$(grep -c '\\newcommand' "$OUT_FILE" || true)
    echo ""
    echo "═══════════════════════════════════════════════════════════════════"
    echo "  Raw \\newcommand lines : $n_macros_raw"
    echo "  After deduplication  : $n_macros  (in $OUT_FILE)"
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

run_script "GROUP 2: Cluster Asymmetry (Tests 1 & 3) + conditional sensitivity subtest" \
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

run_script "GROUP 20d: Leave-2-out / Leave-3-out Stability (dome corpus)" \
           analysis/unesco/dome_leave_k_out.py

echo "─── GROUP 23: Wikidata P1435 Control Analysis ───"
if [ "$HAS_P1435" = true ]; then
    python3 analysis/global/wikidata_p1435_control_analysis.py 2>&1 | tail -40
else
    echo "  SKIPPED — P1435 CSV not found"
fi
echo ""

run_script "GROUP 24: Dome Periodicity Audit"                      \
           analysis/global/dome_periodicity_audit.py

echo "─── GROUP 24b: x.18° Periodicity Formal Test (may take ~30 seconds) ───"
python3 analysis/global/x18_periodicity_formal_test.py 2>&1 | tail -40
echo ""

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

run_script "GROUP 26b: Dome footprint window sensitivity (±2°/±10°)"  \
           analysis/unesco/dome_footprint_window_sensitivity.py

run_script "GROUP 26c: Stupa coordinate-perturbation sensitivity"      \
           analysis/unesco/stupa_coordinate_perturbation.py

run_script "GROUP 26d: Americas UNESCO directional test"               \
           analysis/unesco/americas_directional_test.py

run_script "GROUP 27: Unit Sweep (spacing sensitivity)"            \
           analysis/unesco/unit_sweep_fill.py

echo "─── GROUP 28: Sensitivity Slope Permutation Test (may take several minutes) ───"
python3 analysis/unesco/sensitivity_slope_permutation_test.py 2>&1 | tail -40
echo ""

echo "─── GROUP 28b: Sensitivity Slope Specificity Test — canonical-unit rank (may take several minutes) ───"
python3 analysis/unesco/sensitivity_slope_specificity_test.py 2>&1 | tail -50
echo ""

echo "─── GROUP 29: Dome Geographic-Concentration Null (may take ~2 minutes) ───"
python3 analysis/unesco/dome_geographic_concentration_test.py 2>&1 | tail -40
echo ""

echo "─── GROUP 29b: Stupa Geographic-Concentration Null (may take ~2 minutes) ───"
python3 analysis/unesco/stupa_geographic_concentration_test.py 2>&1 | tail -40
echo ""

run_script "GROUP 29c: Dome × Founding-Category Stratification"   \
           analysis/unesco/dome_founding_stratification.py

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
