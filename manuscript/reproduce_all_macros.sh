#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# reproduce_all_macros.sh
# ══════════════════════════════════════════════════════════════════════════════
#
# Runs every analysis script that produces numeric macros for
# paper_a_primary_unesco.tex, in dependency order.
#
# USAGE:
#   cd /path/to/gerizim-paper-a
#   bash manuscript/reproduce_all_macros.sh
#
# PREREQUISITES:
#   - Python 3.10+ with numpy, scipy, pandas
#   - pip install -r requirements.txt
#   - UNESCO XML:  data/store/unesco/unesco.xml  (fetched by data/scripts/fetch_extended.py)
#   - P1435 CSV:   data/store/wikidata/p1435_global_control.csv
#                  (fetched by data/scripts/fetch_p1435_global.py)
#
# Each script prints its macro values to stdout.  Redirect to a file to
# compare against the values in paper_a_primary_unesco.tex.
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

# ── Data prerequisites ────────────────────────────────────────────────────────
echo "Checking data prerequisites..."
if [ ! -f "data/store/unesco/unesco.xml" ]; then
    echo "  ERROR: data/store/unesco/unesco.xml not found."
    echo "  Run: python3 data/scripts/fetch_extended.py"
    exit 1
fi
echo "  ✓ UNESCO XML found"

if [ ! -f "data/store/wikidata/p1435_global_control.csv" ]; then
    echo "  WARNING: data/store/wikidata/p1435_global_control.csv not found."
    echo "  Groups 23 and 26 will fail."
    echo "  Run: python3 data/scripts/fetch_p1435_global.py"
fi
echo ""

# ── Group 1: Dome / Spherical Monument Raw Sweep (Test 2) ────────────────────
echo "─── GROUP 1: Dome/Spherical Monument Raw Sweep (Test 2) ───"
python3 analysis/unesco/spherical_monument_raw_sweep.py 2>&1 | tail -60
echo ""

# ── Group 2: Cluster Asymmetry (Test 3) ──────────────────────────────────────
echo "─── GROUP 2: Cluster Asymmetry (Test 3) ───"
python3 analysis/unesco/cluster_asymmetry_test.py 2>&1 | tail -60
echo ""

# ── Group 3: Harmonic Density Attractor ──────────────────────────────────────
echo "─── GROUP 3: Harmonic Density Attractor ───"
python3 analysis/unesco/harmonic_density_attractor_test.py 2>&1 | tail -40
echo ""

# ── Group 4: Buddhist Heritage ───────────────────────────────────────────────
echo "─── GROUP 4: Buddhist Heritage ───"
python3 analysis/unesco/unesco_buddhist_heritage_test.py 2>&1 | tail -60
echo ""

# ── Group 5: Founding / Origin (Test E) ──────────────────────────────────────
echo "─── GROUP 5: Origin Sites (Test E) ───"
python3 analysis/unesco/origin_sites_test.py 2>&1 | tail -60
echo ""

# ── Group 6: Founding Sites Analysis (keyword enrichment) ────────────────────
echo "─── GROUP 6: Founding Sites Analysis ───"
python3 analysis/unesco/founding_sites_analysis.py 2>&1 | tail -40
echo ""

# ── Group 7: Sacred Origin Test ──────────────────────────────────────────────
echo "─── GROUP 7: Sacred Origin Test ───"
python3 analysis/unesco/sacred_origin_test.py 2>&1 | tail -20
echo ""

# ── Group 8: Meta-Keyword Test ───────────────────────────────────────────────
echo "─── GROUP 8: Meta-Keyword Test ───"
python3 analysis/unesco/meta_keyword_test.py 2>&1 | tail -20
echo ""

# ── Group 9: Temporal Gradient (Test 4) ──────────────────────────────────────
echo "─── GROUP 9: Temporal Gradient (Test 4) ───"
python3 analysis/unesco/temporal_gradient_test.py 2>&1 | tail -40
echo ""

# ── Group 10: Deep Temporal Analysis ─────────────────────────────────────────
echo "─── GROUP 10: Deep Temporal Analysis ───"
python3 analysis/unesco/deep_temporal_analysis.py 2>&1 | tail -80
echo ""

# ── Group 11: Global Anchor Sweep ────────────────────────────────────────────
echo "─── GROUP 11: Global Anchor Sweep ───"
python3 analysis/global/anchor_uniqueness_audit.py 2>&1 | tail -40
echo ""

# ── Group 12: x.18°E Periodicity Verification ───────────────────────────────
echo "─── GROUP 12: x.18° Periodicity ───"
python3 analysis/unesco/verify_x18_periodicity.py 2>&1 | tail -20
echo ""

# ── Group 12 (cont): Peak Geography Audit ────────────────────────────────────
echo "─── GROUP 12 (cont): Peak Geography Audit ───"
python3 analysis/global/peak_geography_audit.py 2>&1 | tail -30
echo ""

# ── Group 13: Landmark Anchor Ranking ────────────────────────────────────────
echo "─── GROUP 13: Landmark Anchor Ranking ───"
python3 analysis/global/landmark_anchor_ranking.py 2>&1 | tail -20
echo ""

# ── Group 14: Global Corridor Comparison ─────────────────────────────────────
echo "─── GROUP 14: Global Corridor Comparison ───"
python3 analysis/global/global_corridor_comparison.py 2>&1 | tail -30
echo ""

# ── Group 15: Corridor Precision Test ────────────────────────────────────────
echo "─── GROUP 15: Corridor Precision Test ───"
python3 analysis/global/corridor_precision_test.py 2>&1 | tail -60
echo ""

# ── Group 16 + 22: Anchor Site Comparison (Lumbini, Takht, Khoja, WPP) ──────
echo "─── GROUP 16 + 22: Anchor Site Comparison ───"
python3 analysis/global/anchor_site_comparison.py 2>&1
echo ""

# ── Group 17: Spatial Independence / Autocorrelation ─────────────────────────
echo "─── GROUP 17: Spatial Independence ───"
python3 analysis/unesco/spatial_independence_test.py 2>&1 | tail -30
echo ""

# ── Group 18: Regional Temporal Gradient (Mantel-Haenszel) ───────────────────
echo "─── GROUP 18: Regional Temporal Gradient ───"
python3 analysis/unesco/regional_temporal_gradient.py 2>&1 | tail -20
echo ""

# ── Group 19: FDR Multiple Comparisons ───────────────────────────────────────
echo "─── GROUP 19: FDR Multiple Comparisons ───"
python3 analysis/unesco/fdr_multiple_comparisons.py 2>&1 | tail -20
echo ""

# ── Group 20: Hemispherical Mound Evolution Raw Sweep (Test 2b) ──────────────
echo "─── GROUP 20: Hemispherical Evolution Raw Sweep (Test 2b) ───"
python3 analysis/unesco/tumulus_dome_evolution_raw_sweep.py 2>&1 | tail -40
echo ""

# ── Group 21: Bonferroni-Adjusted p-values ───────────────────────────────────
echo "─── GROUP 21: Bonferroni Correction ───"
python3 analysis/unesco/bonferroni_correction.py 2>&1 | tail -20
echo ""

# ── Group 23: Wikidata P1435 Control Analysis ────────────────────────────────
echo "─── GROUP 23: Wikidata P1435 Control Analysis ───"
if [ -f "data/store/wikidata/p1435_global_control.csv" ]; then
    python3 analysis/global/wikidata_p1435_control_analysis.py 2>&1
else
    echo "  SKIPPED — P1435 CSV not found"
fi
echo ""

# ── Group 24: Dome Periodicity Audit ─────────────────────────────────────────
echo "─── GROUP 24: Dome Periodicity Audit ───"
python3 analysis/global/dome_periodicity_audit.py 2>&1 | tail -40
echo ""

# ── Group 25: Simulation-Based Null Model ────────────────────────────────────
echo "─── GROUP 25: Simulation Null Model (this may take a few minutes) ───"
python3 analysis/unesco/simulation_null_model.py 2>&1 | tail -40
echo ""

# ── Group 26: Americas P1435 Control ─────────────────────────────────────────
echo "─── GROUP 26: Americas Harmonic Depletion ───"
if [ -f "data/store/wikidata/p1435_global_control.csv" ]; then
    python3 analysis/americas/americas_harmonic_depletion_audit.py 2>&1 | tail -30
else
    echo "  SKIPPED — P1435 CSV not found"
fi
echo ""

# ── Group 27: Unit Sweep (Tables 6-7) ────────────────────────────────────────
echo "─── GROUP 27: Unit Sweep (spacing sensitivity) ───"
python3 analysis/unesco/unit_sweep_fill.py 2>&1 | tail -40
echo ""

echo "═══════════════════════════════════════════════════════════════════"
echo "  ALL SCRIPTS COMPLETE"
echo "  Compare the output above with the macro values in"
echo "  manuscript/paper_a_primary_unesco.tex"
echo "═══════════════════════════════════════════════════════════════════"
