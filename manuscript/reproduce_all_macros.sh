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
#   bash manuscript/reproduce_all_macros.sh          # fast: LaTeX macros only (default)
#   bash manuscript/reproduce_all_macros.sh --full   # slow: runs all analysis scripts
#
# PREREQUISITES:
#   - Python 3.10+ with numpy, scipy, pandas
#   - pip install -r requirements.txt
#   - UNESCO XML:  data/store/unesco/unesco.xml
#   - UNESCO XML: data/store/unesco/unesco.xml
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

# ── default (no --full): run all scripts, collect only \newcommand lines ─────
if [ "${1:-}" != "--full" ]; then
    OUT_FILE="manuscript/generated_macros.tex"
    echo "% Generated macros -- $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$OUT_FILE"
    echo "% Source: bash manuscript/reproduce_all_macros.sh" >> "$OUT_FILE"
    echo "% IMPORTANT: values are computed fresh each run -- do NOT hardcode" >> "$OUT_FILE"
    echo "" >> "$OUT_FILE"

    # DEPENDENCY ORDER:
    #   Primary analysis scripts first (they write to the store)
    #   Summary scripts last (they read from the store)
    SCRIPTS=(
        # ── Constants (config.json + corpus size) — no corpus analysis needed ──
        analysis/global/emit_constants.py                    # → GerizimLon, NwhcTotal, etc.
        analysis/global/vonmises_threshold_validation.py     # → vmKappaOwtrad, vmKappaStupa, vmSigmaBeruStupa, etc.
        analysis/global/geodesic_sensitivity.py              # → GeodesicApCurrent, GeodesicDropOut, GeodesicGainIn, etc.
        analysis/global/owtrad_route_alignment.py            # → pOwtradMidApp/Ap, NOwtradEdges/Vertices, etc.
        analysis/global/clustered_control_sweep.py           # → clusterCtrlPermP, clusterCtrlNsim, clusterCtrlTobs, etc.
        # ── Primary analysis (write to store) ─────────────────────────────
        analysis/unesco/spherical_monument_test.py            # → pCircApValidated, NcircValidated (context-validated, Exploratory 2x) — runs FIRST so raw_sweep can read it
        analysis/unesco/spherical_monument_raw_sweep.py      # → pCircAp, pCircA, pCircChi (reads validated counts from store for comparison table)
        analysis/unesco/cluster_asymmetry_test.py            # → clusterApBinom, clusterPermP, clusterHarmonicMWpStr
        analysis/global/emit_cluster_ap_ci.py                # → clusterApCIlo, clusterApCIhi
        analysis/unesco/harmonic_density_attractor_test.py
        analysis/unesco/origin_sites_test.py                 # → pCanon, pPreTwoK, pModern, pReligUnion, pChrist, pBudRelig, pJudATierFisher
        analysis/unesco/founding_sites_analysis.py           # → pFoundKwFisher
        analysis/unesco/sacred_origin_test.py
        analysis/unesco/meta_keyword_test.py
        analysis/unesco/keyword_sensitivity_stupa.py         # → kwSensStupa/Mound/NoStupa N/App/Ap counts+p
        analysis/unesco/temporal_gradient_test.py            # → pCochranThree, pCochranFive
        analysis/unesco/deep_temporal_analysis.py            # → peakSigP, pFoundDateSpearman
        analysis/global/anchor_uniqueness_audit.py
        analysis/unesco/verify_phase_peak_periodicity.py
        analysis/global/peak_geography_audit.py
        analysis/global/landmark_anchor_ranking.py
        analysis/global/anchor_site_comparison.py            # → NjerSelfExclCorpus, NgerizimSweepCorpus
        analysis/unesco/spatial_independence_test.py         # → pNeffQuarter, pNeffHalf
        analysis/unesco/regional_temporal_gradient.py        # → pCMH
        analysis/unesco/tumulus_dome_evolution_test.py        # → pEvoAp_validated, NevoValidTotal (context-validated, Exploratory 2bx) — runs FIRST so raw_sweep can read it
        analysis/unesco/tumulus_dome_evolution_raw_sweep.py  # → pEvoAp (Test 2b) (reads validated counts from store for comparison table)
        analysis/unesco/mound_keyword_context_audit.py        # → NmoundRaw/Accepted/Rejected, moundFPRate, moundKwFPRate, pmoundAccAp, moundEnrichAcc
        analysis/unesco/leave_one_out_sensitivity.py         # → LOOdomeN/Ap/Rate/Enrich/P, LOOstupaN/Ap/Rate/P
        analysis/unesco/dome_leave_k_out.py                  # → LKOsiteOne/Two/ThreeName/DevKm, LKOtwo*/three*/worstTwo/ThreeP
        analysis/unesco/evo_leave_k_out.py                   # → EvoLKOnCombosTwo, EvoLKOnCombosThree
        analysis/global/dome_periodicity_audit.py
        analysis/global/phase_peak_periodicity_formal_test.py       # → rayleighR/Z/PermP, fullRayleighR/Z/PermP, targetedPdome/Wdome/Kdome (GROUP 11b)
        analysis/unesco/periodogram_test.py                          # → periodogramPermP, periodogramPeakTDome, periodogramRank3Dome, etc.
        analysis/unesco/joint_periodogram_test.py                    # → npcMonK, npcMonPermP, npcMonGridDeltaE, npcOwR, npcOwRayleighP, etc.
        analysis/unesco/circular_two_sample_tests.py                 # → circConcStupaR/P, circConcUnescoR/P, circWatsonPUnescoStupa, circKuiperPUnescoStupa, etc.
        analysis/global/phase_peak_max_permutation_test.py          # → anchorMaxPermObsMax/NullMu/NullSD/Z/P/Nperms/BootMu/BootSD/BootZ/BootP (GROUP 11c)
        analysis/global/phase_envelope_analysis.py                  # → phaseFull*/phaseDome*/phaseGerizimPhase (phase-marginalised robustness)
        # ── Reviewer robustness checks ────────────────────────────────────
        analysis/unesco/dome_footprint_window_sensitivity.py   # → geoNullDomeTwoDeg*, geoNullDomeTenDeg*
        analysis/unesco/stupa_coordinate_perturbation.py       # → stupaCoordPerturb*
        analysis/unesco/americas_directional_test.py           # → AmericasN, AmericasApCount, AmericasOneSidedP, AmericasDirectional
        analysis/unesco/wikidata_q180987_stupa_audit.py        # → wikiStupaTotal, wikiStupaATierCount/Rate/Enrich/BinomP, wikiStupaPermZ/P, wikiStupaAp*, wikiStupaJava*, wikiStupaCluster*, named-site dev_km
        analysis/unesco/multiscale_combined_p.py               # → multiscaleDomeCombinedP, multiscaleJointMaxP, multiscaleNaturalKm, etc.
        # ── Simulation (slow — writes simDomePermP etc.) ──────────────────
        analysis/unesco/simulation_null_model.py
        # ── Cluster asymmetry + conditional sensitivity subtest ───────────
        # (cluster_asymmetry_test.py now includes the conditional subtest at the end;
        #  it emits clusterCondN, clusterCondClusterOk, clusterCondPeakFrac, etc.)

        analysis/unesco/unit_sweep_fill.py
        analysis/unesco/unit_sweep_montecarlo.py              # → mcDomeNsigObs, mcDomePObsGe, mcFullNsigObs, mcNperms
        analysis/unesco/sensitivity_slope_permutation_test.py  # → permSlopeNperms, permSlopePcanon, permSlopePsharp
        analysis/unesco/sensitivity_slope_specificity_test.py  # → permSlopeCanonBestJoint{APP,AP,A}, permSlopeCanonRankObs{APP,AP,A}, permSlopeCanonCondFraction{APP,AP,A}
        analysis/unesco/dome_geographic_concentration_test.py  # → geoNullDomeBootP/Z/Mean, geoNullDomeRestrictedP/Z/Mean/N, nullBBonfK/Alpha, domeEurasianFraction
        analysis/unesco/stupa_geographic_concentration_test.py # → stupaGeoBootP/Z/Mean, stupaRegionP/Z/N, stupaGeoRestrictedP/Z/N (GROUP 29b)
        analysis/unesco/osm_stupa_audit.py                    # → osmStupaTotal, osmStupaRayleighR/P, osmStupaATierCount, etc. (GROUP 29d)
        analysis/unesco/cross_corpus_geometric_test.py        # → geoDome/Wiki/OsmApp/Ap/A N/Hits/Exp/Enr/P/Sig, geoFisherApp/Ap/A Chi/P/Sig (PRIMARY)
        analysis/unesco/dome_founding_stratification.py        # → domeStratNfs, domeStratNnfs, domeStratF/NfsAp/Rate/Enrich/P, domeStratFisherOR/P
        # ── Geographic concentration robustness ───────────────────────────
        analysis/unesco/region_conditioned_permutation.py    # → regionCondPermP/Z/Obs/Mean/Std/Nregions/Nperms (GROUP 35)
        # ── Global corridor / anchor comparison ───────────────────────────
        analysis/global/anchor_corridor_audit.py             # → corridorSharedAp, corridorUniqueGer, corridorUniqueJer
        # ── Tier log-sweep sensitivity ────────────────────────────────────
        analysis/unesco/tier_logsweep_sensitivity.py         # → logSweepNtaus, logSweepAregionNsig/Ntaus, logSweepOpt*
        # ── Summary scripts (read from store) ─────────────────────────────
        analysis/unesco/fdr_multiple_comparisons.py          # reads all keys from config+store
        analysis/unesco/bonferroni_correction.py             # reads confirmatory keys from store
        # ── Sig-label companions (MUST be last — reads all p-values from store) ─
        analysis/global/emit_sig_macros.py                   # → \macroNameSig auto-generated from store
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
_RE = re.compile(r'\\newcommand\{(\\[A-Za-z@]+)\}')
seen = set()
kept = []
for line in reversed(lines):
    m = _RE.search(line)
    if m:
        name = m.group(1)
        if name in seen:
            continue
        seen.add(name)
    kept.append(line)
kept.reverse()
path.write_text("".join(kept))
print(f"  Deduplicated: kept {len(seen)} unique macros (removed {len(lines)-len(kept)} duplicate lines)")
PYEOF

    # ── Strip non-ASCII from comment portions of generated_macros.tex ─────────
    # LaTeX reports "Missing character" warnings for Unicode in % comments even
    # when inputenc is loaded.  Sanitise by replacing non-ASCII chars in comment
    # tails (after the first %) with their closest ASCII equivalents.
    python3 - "$OUT_FILE" <<'PYEOF2'
import sys, re
from pathlib import Path

REPLACEMENTS = {
    '\u2014': '--',   # em dash
    '\u2013': '-',    # en dash
    '\u2212': '-',    # minus sign
    '\u00B0': 'deg',  # degree sign
    '\u00D7': 'x',    # multiplication sign
    '\u2265': '>=',   # greater-or-equal
    '\u2264': '<=',   # less-or-equal
    '\u00E9': 'e',    # e acute
    '\u0113': 'e',    # e macron
    '\u00A7': 'S',    # section sign
}

path = Path(sys.argv[1])
lines = path.read_text(encoding='utf-8').splitlines(keepends=True)
out = []
for line in lines:
    # Find the first % that starts a comment (not inside {})
    pct = line.find('%')
    if pct != -1:
        before = line[:pct]
        comment = line[pct:]
        for ch, asc in REPLACEMENTS.items():
            comment = comment.replace(ch, asc)
        # Drop any remaining non-ASCII in comments
        comment = comment.encode('ascii', 'replace').decode('ascii')
        line = before + comment
    out.append(line)
path.write_text(''.join(out), encoding='utf-8')
print(f"  Sanitised non-ASCII from comment lines in {path.name}")
PYEOF2

    n_macros=$(grep -c '\\newcommand' "$OUT_FILE" || true)
    echo ""
    echo "═══════════════════════════════════════════════════════════════════"
    echo "  Raw \\newcommand lines : $n_macros_raw"
    echo "  After deduplication  : $n_macros  (in $OUT_FILE)"
    echo "═══════════════════════════════════════════════════════════════════"
    echo ""
    echo "  Regenerating manuscript figures to match updated macros..."
    python3 manuscript/generate_figures.py 2>&1 | grep -E "✓|Error|WARNING" || true
    echo ""
    echo "  ✓ Figures regenerated. Caption annotations and image annotations"
    echo "    now derive from the same live pipeline run."
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

run_script "GROUP 0c: Von Mises threshold validation (§3.4)"       \
           analysis/global/vonmises_threshold_validation.py

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

run_script "GROUP 4: Origin Sites (canon/pre-2000/modern/religion)" \
           analysis/unesco/origin_sites_test.py

run_script "GROUP 5: Founding Sites Analysis (keyword enrichment)" \
           analysis/unesco/founding_sites_analysis.py

run_script "GROUP 7: Sacred Origin Test"                           \
           analysis/unesco/sacred_origin_test.py

run_script "GROUP 8: Meta-Keyword Test"                            \
           analysis/unesco/meta_keyword_test.py

run_script "GROUP 8b: Keyword Sensitivity / Leave-Stupa-Out"       \
           analysis/unesco/keyword_sensitivity_stupa.py

run_script "GROUP 9: Temporal Gradient (Test 4)"                   \
           analysis/unesco/temporal_gradient_test.py

run_script "GROUP 10: Deep Temporal Analysis"                      \
           analysis/unesco/deep_temporal_analysis.py

run_script "GROUP 11: Global Anchor Sweep"                         \
           analysis/global/anchor_uniqueness_audit.py

run_script "GROUP 12: phase-peak° Periodicity Verification"             \
           analysis/unesco/verify_phase_peak_periodicity.py

run_script "GROUP 12b: Peak Geography Audit"                       \
           analysis/global/peak_geography_audit.py

run_script "GROUP 13: Landmark Anchor Ranking"                     \
           analysis/global/landmark_anchor_ranking.py

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

echo "─── GROUP 20c: Mound Keyword Context Audit (FP rate, Test 2b-x) ───"
python3 analysis/unesco/mound_keyword_context_audit.py --latex 2>&1 | tail -60
echo ""

run_script "GROUP 20d: Leave-One-Out Sensitivity (Tests 2 & 2b)" \
           analysis/unesco/leave_one_out_sensitivity.py

run_script "GROUP 20d: Leave-2-out / Leave-3-out Stability (dome corpus)" \
           analysis/unesco/dome_leave_k_out.py

run_script "GROUP 24: Dome Periodicity Audit"                      \
           analysis/global/dome_periodicity_audit.py

echo "─── GROUP 24b: phase-peak° Periodicity Formal Test (may take ~30 seconds) ───"
python3 analysis/global/phase_peak_periodicity_formal_test.py 2>&1 | tail -40
echo ""

echo "─── GROUP 24c: phase-peak° Discovery-Corrected Max-Permutation Test (may take several minutes) ───"
python3 analysis/global/phase_peak_max_permutation_test.py 2>&1 | tail -40
echo ""

echo "─── GROUP 24d: Phase-Envelope Analysis (anchor-free robustness) ───"
python3 analysis/global/phase_envelope_analysis.py 2>&1 | tail -30
echo ""
echo "─── GROUP 25: Simulation Null Model (may take several minutes) ───"
python3 analysis/unesco/simulation_null_model.py 2>&1 | tail -40
echo ""

run_script "GROUP 26b: Dome footprint window sensitivity (±2°/±10°)"  \
           analysis/unesco/dome_footprint_window_sensitivity.py

run_script "GROUP 26c: Stupa coordinate-perturbation sensitivity"      \
           analysis/unesco/stupa_coordinate_perturbation.py

run_script "GROUP 26d: Americas UNESCO directional test"               \
           analysis/unesco/americas_directional_test.py

run_script "GROUP 26e: Wikidata Q180987 stupa corpus portability audit" \
           analysis/unesco/wikidata_q180987_stupa_audit.py

run_script "GROUP 26f: Circular two-sample tests (Dome/Stupa + UNESCO/Stupa)" \
           analysis/unesco/circular_two_sample_tests.py

run_script "GROUP 26g: Multi-scale combined enrichment (4-corpus joint)" \
           analysis/unesco/multiscale_combined_p.py

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

echo "─── GROUP 29a: Null B — Uniform Random Draw (may take ~1 minute) ───"
python3 analysis/unesco/null_b_uniform_draw.py 2>&1 | tail -20
echo ""

echo "─── GROUP 29b: Stupa Geographic-Concentration Null (may take ~2 minutes) ───"
python3 analysis/unesco/stupa_geographic_concentration_test.py 2>&1 | tail -40
echo ""

run_script "GROUP 29c: Dome × Founding-Category Stratification"   \
           analysis/unesco/dome_founding_stratification.py

run_script "GROUP 29d: OSM Stupa Audit" \
           analysis/unesco/osm_stupa_audit.py

run_script "GROUP 29e: Cross-Corpus Geometric Null Test (PRIMARY)" \
           analysis/unesco/cross_corpus_geometric_test.py

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
echo "  REGENERATING MANUSCRIPT FIGURES"
echo "  (must run after macros so figure annotations match pipeline values)"
echo "═══════════════════════════════════════════════════════════════════"
python3 manuscript/generate_figures.py 2>&1 | tail -20
echo ""

echo "═══════════════════════════════════════════════════════════════════"
echo "  ALL SCRIPTS COMPLETE"
echo "  Results store: data/store/results.json"
echo "  To regenerate LaTeX macros file (fast, default):"
echo "    bash manuscript/reproduce_all_macros.sh"
echo ""
echo "  RECOMMENDED FINAL SUBMISSION SEQUENCE:"
echo "    1. bash manuscript/reproduce_all_macros.sh --full   # full pipeline (slow)"
echo "    2. bash manuscript/reproduce_all_macros.sh          # regen macros (fast)"
echo "    3. cd manuscript/archaeometry && pdflatex paper_a_archaeometry.tex     # compile"
echo "    4. pdflatex paper_a_archaeometry.tex                      # resolve refs"
echo "    5. cd ../primary && pdflatex paper_a_primary_unesco.tex   # compile primary"
echo "═══════════════════════════════════════════════════════════════════"
