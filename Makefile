# Makefile — Reproduce all Paper A results
# Usage: make all

PYTHON ?= python3

.PHONY: all test figures analyses primary robustness global clean

all: test analyses figures
	@echo ""
	@echo "✓ All analyses and figures complete."

# ── Unit tests ────────────────────────────────────────────────────────────────
test:
	$(PYTHON) -m pytest tests/ -v

# ── Primary confirmatory tests ────────────────────────────────────────────────
primary:
	@echo "\n══ Test 2: Dome enrichment (confirmatory raw sweep) ══"
	$(PYTHON) analysis/unesco/spherical_monument_raw_sweep.py
	@echo "\n══ Test 2x: Dome enrichment (context-validated) ══"
	$(PYTHON) analysis/unesco/spherical_monument_test.py
	@echo "\n══ Test 3: Cluster asymmetry ══"
	$(PYTHON) analysis/unesco/cluster_asymmetry_test.py
	@echo "\n══ Test 4: Temporal gradient ══"
	$(PYTHON) analysis/unesco/temporal_gradient_test.py
	@echo "\n══ Test E: Founding sites ══"
	$(PYTHON) analysis/unesco/origin_sites_test.py
	$(PYTHON) analysis/unesco/founding_sites_analysis.py
	@echo "\n══ Sacred origin test ══"
	$(PYTHON) analysis/unesco/sacred_origin_test.py
	@echo "\n══ Meta-keyword test ══"
	$(PYTHON) analysis/unesco/meta_keyword_test.py
	@echo "\n══ Buddhist heritage test ══"
	$(PYTHON) analysis/unesco/unesco_buddhist_heritage_test.py
	@echo "\n══ Harmonic density attractor ══"
	$(PYTHON) analysis/unesco/harmonic_density_attractor_test.py
	@echo "\n══ Hemispherical mound evolution ══"
	$(PYTHON) analysis/unesco/tumulus_dome_evolution_raw_sweep.py
	@echo "\n══ Deep temporal analysis ══"
	$(PYTHON) analysis/unesco/deep_temporal_analysis.py

# ── Robustness and corrections ────────────────────────────────────────────────
robustness:
	@echo "\n══ Unit sensitivity sweep ══"
	$(PYTHON) analysis/unesco/unit_sweep_fill.py
	@echo "\n══ Bonferroni correction ══"
	$(PYTHON) analysis/unesco/bonferroni_correction.py
	@echo "\n══ FDR multiple comparisons ══"
	$(PYTHON) analysis/unesco/fdr_multiple_comparisons.py
	@echo "\n══ Spatial independence ══"
	$(PYTHON) analysis/unesco/spatial_independence_test.py
	@echo "\n══ Regional temporal gradient ══"
	$(PYTHON) analysis/unesco/regional_temporal_gradient.py
	@echo "\n══ Simulation null model ══"
	$(PYTHON) analysis/unesco/simulation_null_model.py
	@echo "\n══ x.18° periodicity verification ══"
	$(PYTHON) analysis/unesco/verify_x18_periodicity.py

# ── Global checks ─────────────────────────────────────────────────────────────
global:
	@echo "\n══ Global anchor uniqueness audit ══"
	$(PYTHON) analysis/global/anchor_uniqueness_audit.py
	@echo "\n══ Peak geography audit ══"
	$(PYTHON) analysis/global/peak_geography_audit.py
	@echo "\n══ Dome periodicity audit ══"
	$(PYTHON) analysis/global/dome_periodicity_audit.py
	@echo "\n══ Global corridor comparison ══"
	$(PYTHON) analysis/global/global_corridor_comparison.py
	@echo "\n══ Corridor precision test ══"
	$(PYTHON) analysis/global/corridor_precision_test.py
	@echo "\n══ Landmark anchor ranking ══"
	$(PYTHON) analysis/global/landmark_anchor_ranking.py
	@echo "\n══ Americas depletion control ══"
	$(PYTHON) analysis/americas/americas_harmonic_depletion_audit.py

# ── All analyses ──────────────────────────────────────────────────────────────
analyses: primary robustness global

# ── Figures ───────────────────────────────────────────────────────────────────
figures:
	@echo "\n══ Generating manuscript figures ══"
	$(PYTHON) manuscript/generate_figures.py

# ── LaTeX compilation (optional — requires pdflatex) ──────────────────────────
pdf:
	cd manuscript && pdflatex paper_a_primary_unesco.tex

# ── Clean generated files ─────────────────────────────────────────────────────
clean:
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	rm -f manuscript/*.aux manuscript/*.log manuscript/*.out manuscript/*.toc
