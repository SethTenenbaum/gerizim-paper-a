# Supplementary S1 — Analysis Script Index

This file documents every analysis, audit, and figure-generation script in
the repository, grouped by role.  All scripts read from and write to
`data/store/results.json` (the shared results store) unless noted otherwise.
The **"Macros it produces"** column lists representative LaTeX macro names
emitted to `manuscript/generated_macros.tex`; these can be cross-referenced
against the macro file to confirm which section of the paper uses each value.
Macro regeneration and figure generation are fully automated; no results are
entered by hand.

---

## Reproduction Quick-Start

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Run the full macro pipeline
bash manuscript/reproduce_all_macros.sh

# 3. Regenerate all audit documents
bash tools/update_all_audits.sh

# 4. Regenerate figures
python manuscript/generate_figures.py
```

Individual scripts can be run directly with `python <path>` from the
repository root; each prints a summary and writes results to the store.

---

## Core Statistical Tests (`analysis/unesco/`)

| Script | Macros it produces | Paper section |
|--------|-------------------|---------------|
| `spherical_monument_test.py` | `pCircApValidated`, `pCircAppValidatedFisher`, `NcircValidatedAp`, … | §4.2 Domed and Spherical Monuments (Test 2, Primary) |
| `spherical_monument_raw_sweep.py` | `NcircTotal`, `pCircAp`, `pCircApp`, `circEnrichAp`, `NcircTierApp`, … | §4.2 Domed and Spherical Monuments (Test 2, Primary) |
| `cluster_asymmetry_test.py` | `clusterApBinom`, `clusterPermP`, `clusterHarmonicRatio`, `clusterHarmonicMWp`, … | §4.4 Cluster Asymmetry (Test 3, Primary) |
| `origin_sites_test.py` | `pOriginAll`, `clusterApRate`, full-corpus cohort p-values, all religion-group p-values | §4.1 Global Enrichment (Test 1, Primary) |
| `temporal_gradient_test.py` | `ZcochranThree`, `pCochranThree`, `pCanonApp`, cohort-rate macros | §4.5 Temporal Gradient (Test 4, Descriptive) |
| `regional_temporal_gradient.py` | `pLogRegCohort`, `MHcommonOR`, `pCMH`, regional A+ rates | §4.5 Temporal Gradient (Test 4, Descriptive) |
| `deep_temporal_analysis.py` | `NdatedSites`, `pPreCEbinom`, `peakSigP`, `pFirstHalf`, block p-values | §4.5 Temporal Gradient (Test 4, Descriptive) |
| `tumulus_dome_evolution_test.py` | *(no macros — raw output only; see raw_sweep below)* | §4.3 Hemispherical Mound Evolution (Test 2b) |
| `tumulus_dome_evolution_raw_sweep.py` | `NevoTotal`, `pEvoAp`, `pEvoApp`, `evoStupaAppFisherOR`, tier-by-stage p-values | §4.3 Hemispherical Mound Evolution (Test 2b) |
| `dome_geographic_concentration_test.py` | `geoNullDomeRestrictedP`, `geoNullDomeBootP` (Null A / Null C) | §4.2 Domed and Spherical Monuments (Test 2) / §6.2 |
| `stupa_geographic_concentration_test.py` | `stupaGeoBootP`, `stupaGeoRestrictedP` (stupa Null C) | §4.3 Hemispherical Mound Evolution (Test 2b) |
| `harmonic_density_attractor_test.py` | `attractorMWp`, `attractorHarmonicRatio`, `attractorHarmonicMWp` | §5.4 Spatial Autocorrelation |
| `spatial_independence_test.py` | `DEFFhalf`, `NeffHalf`, `pNeffHalf`, `blockBootCIlo`, `blockBootCIhi` | §5.4 Spatial Autocorrelation |
| `simulation_null_model.py` | `simDomePermP`, `simDomeBootP`, `simKDEZ`, `simKDEP`, `simNperms` | §5.5 Null Model Validation |
| `unit_sweep_fill.py` | All `sweepZeroX…` macros (0.05–0.20 beru hits and p-values) | §5.1 Unit Sensitivity Sweep |
| `unit_sweep_montecarlo.py` | `mcDomePObsGe`, `mcDomeNsigCollapsed`, `mcDomeCollapsedPEqOne` | §6.2 Geographic and Cultural Selection Bias |
| `sensitivity_slope_permutation_test.py` | `permSlopePcanon`, `permSlopePsharp` | §5.2 Sensitivity Slope at the Canonical Unit |
| `sensitivity_slope_specificity_test.py` | `permSlopeCanonBestJointDomeAPP`, `pEvoApValidated`, `NevoValidTotal`, … | §5.2 Sensitivity Slope at the Canonical Unit |
| `tier_logsweep_sensitivity.py` | `logSweepNtaus`, `logSweepNsig`, `logSweepOptApBeru`, `logSweepOptRatio` | §5.2 Sensitivity Slope at the Canonical Unit |
| `leave_one_out_sensitivity.py` | `LOOdomeN`, `LOOdomeP`, `LOOstupaP` | §4.2 Domed and Spherical Monuments (Test 2, Primary) |
| `dome_leave_k_out.py` | `LKOtwoP`, `LKOthreeP`, `LKOworstTwoP`, `LKOworstThreeP` | §4.2 Domed and Spherical Monuments (Test 2, Primary) |
| `evo_leave_k_out.py` | `EvoLKOworstTwoP`, `EvoLKOworstThreeP` | §4.3 Hemispherical Mound Evolution (Test 2b) |
| `stupa_coordinate_perturbation.py` | `stupaCoordPerturbMean`, `stupaCoordPerturbSigPct` | §5.6 Pipeline Portability: Wikidata Q180987 Stupa Corpus |
| `wikidata_q180987_stupa_audit.py` | `wikiStupaTotal`, `wikiStupaATierBinomP`, `wikiJavaNodeATierP`, … | §5.6 Pipeline Portability: Wikidata Q180987 Stupa Corpus |
| `americas_directional_test.py` | `AmericasN`, `AmericasApRate`, `AmericasOneSidedP` | §6.3 Americas Negative Control |
| `region_conditioned_permutation.py` | `regionCondPermP`, `regionCondPermMean` | §5.4 Spatial Autocorrelation |
| `founding_sites_analysis.py` | `NthematicFounding`, `foundKwFisherOR`, `pFoundKwFisher` | Appendix A.2 (Test 5, post-hoc) |
| `dome_founding_stratification.py` | `domeStratFsP`, `domeStratFisherOR` | Appendix A.2 (Test 5, post-hoc) |
| `sacred_origin_test.py` | `sacredHitOne`, `sacredHitTwo`, `sacredHitThree` | §4.2 Domed and Spherical Monuments (Test 2, Primary) |
| `meta_keyword_audit.py` | *(audit output only, no macros)* | Appendix A.1 Dome/Spherical Filter |
| `meta_keyword_test.py` | `NmetaKeywords`, `NmetaBeatingChi`, `metaBeatingChiPct` | Appendix A.1 Dome/Spherical Filter |
| `mound_keyword_context_audit.py` | `NmoundRaw`, `NmoundAccepted`, `moundFPRate`, `pmoundRawAp` | Appendix A.1 / §4.3 |
| `corpus_exclusion_audit.py` | *(audit output only, no macros)* | §3.1 Corpus Construction |
| `fine_sweep_audit.py` | *(no standalone macros — superseded by unit_sweep_fill)* | §5.2 Sensitivity Slope at the Canonical Unit |
| `dome_footprint_window_sensitivity.py` | `geoNullDomeTwoDegP`, `geoNullDomeTenDegP` | §4.2 Domed and Spherical Monuments (Test 2, Primary) |
| `bonferroni_correction.py` | `pAdjTestOne`, `pAdjTestTwo`, `pAdjTestThree`, `BonfK` | §3.6 Multiple-Comparisons Protocol + §4.6 FDR |
| `fdr_multiple_comparisons.py` | `NtotalTests`, `NsurviveFDR`, `NsurviveBonfAll`, all `FDRqSurvive…` | §4.6 FDR Audit and Bonferroni Correction |

---

## Global / Anchor Analyses (`analysis/global/`)

| Script | Macros it produces | Paper section |
|--------|-------------------|---------------|
| `emit_constants.py` | `GerizimLon`, `NwhcTotal`, `NclusterTotal`, `ApThreshBeru`, `NullRateAp`, tier thresholds, corpus-region %s | §3.1–3.4 (constants used throughout) |
| `geodesic_sensitivity.py` | `GeodesicApCurrent`, `GeodesicApPhysical`, `GeodesicDropOut`, `GeodesicGainIn` | §6.4 Other Limitations |
| `owtrad_route_alignment.py` | `NOwtradEdges`, `NOwtradVertices`, `pOwtradVertexApp`, `pOwtradDegWApp`, `owtradClustDegRatio`, … | §7.2 Architectural Stratification and Route-Segment Alignment |
| `emit_cluster_ap_ci.py` | `clusterApCIlo`, `clusterApCIhi` | §4.4 Cluster Asymmetry (Test 3, Primary) |
| `anchor_uniqueness_audit.py` | `anchorSweepNanchors`, `GerizimSweepAp`, `anchorSweepPctile`, `JerusalemSweepAp`, `JerusalemSweepPctile` — all cited | §5.3 Global Anchor Sweep |
| `verify_phase_peak_periodicity.py` | `xbandLabel`, `optimalPhase`, `GerizimPhase`, `phaseGap`, `phaseGapKm` — all cited | §5.3 Global Anchor Sweep |
| `peak_geography_audit.py` | `JerichoAnchorAp`, `DeadSeaAnchorAp`, `MtNeboAnchorAp`, `DamascusAnchorAp`, `localBestLon`, `localBestAp` — **stdout-only; not written to results store; none currently cited in manuscript prose** (audit/diagnostic) | §5.3 Global Anchor Sweep (audit only) |
| `landmark_anchor_ranking.py` | `GerizimTempleRank` — cited (§5.3); `GerizimTempleAp` — in generated_macros.tex but not cited in prose — **stdout-only output** | §5.3 Global Anchor Sweep |
| `anchor_site_comparison.py` | **Stdout-only macros (not cited):** `LumbiniDevKm`, `TakhtGerizimDevKm`, `KhojaGerizimDevBeru`, `topHit*`, `appSite*` — printed to stdout, not in results store, not cited in manuscript. **Results-store macros (cited):** `pMeccaAnchor`, `pMeruAnchor`, `pKailashAnchor` (§3.3 anchor discussion); `NgerizimSweepCorpus`, `JerusalemAp`, `pJerusalemAp` (§5.3) | §3.3 Anchor / §5.3 Global Anchor Sweep |
| `anchor_corridor_audit.py` | `corridorSharedAp`, `corridorUniqueGer`, `corridorUniqueJer` — all cited | §5.3 Global Anchor Sweep |
| `dome_periodicity_audit.py` | `domePhaseChiSq`, `domePhaseChiP`, `domeGerizimPctile` — all cited | §5.3 Global Anchor Sweep |
| `phase_peak_periodicity_formal_test.py` | `fullRayleighR`, `fullRayleighZ`, `fullRayleighPermP`, `anchorShiftPermP`, `anchorShiftZ`, `anchorShiftNullMu`, `anchorShiftApCount`, `targetedPfull`, `targetedPAp`, `targetedDfull`, `targetedWfull`, `targetedKfull` (Tests A, C, D) — all cited | §4.4 Cluster Asymmetry (Tests A/C/D) / §5.3 |
| `phase_peak_max_permutation_test.py` | `anchorMaxPermP` — cited in **§4.4 (Test E)**, **§5.3**, and **§6.1**; `anchorMaxPermZ` — in generated_macros.tex but not cited in prose | §4.4 Cluster Asymmetry (Test E) / §5.3 / §6.1 |
| `emit_sig_macros.py` | All `…Sig` star macros (aggregated from store) | Used throughout Results |

---

## Americas Analysis (`analysis/americas/`)

| Script | Role | Paper section |
|--------|------|---------------|
| `americas_harmonic_depletion_audit.py` | Audit of harmonic-node depletion in Americas sub-corpus (control region) | §6.3 Americas Negative Control |

---

## Audit Output Files (`supplementary/audit/`)

All files in `supplementary/audit/` are **auto-generated** — do not edit by
hand.  Regenerate them with `bash tools/update_all_audits.sh`.

| File | Contents | Paper section |
|------|----------|---------------|
| `README_audit.txt` | Overview of the audit archive: lists every file, the script that produced it, and a one-line description | — |
| `anchor_sweep_audit.txt` | Site-by-site A/A+ scores for every candidate anchor in the global sweep; Gerizim percentile rank | §5.3 Global Anchor Sweep |
| `aplus_sites_audit.txt` | Full listing of A+ UNESCO sites: name, coordinates, harmonic distance, tier | §4.1 Global Enrichment (Test 1, Primary) |
| `dome_keyword_audit.txt` | Every UNESCO site matched (or rejected) by the domed-monument keyword classifier; acceptance/rejection reason | §3.1 Corpus Construction / Appendix A.1 Dome/Spherical Filter |
| `dome_mound_keyword_audit.txt` | Site-by-site typology classification into dome vs. mound vs. flat-topped categories | §4.3 Hemispherical Mound Evolution (Test 2b) |
| `fdr.txt` | All \NtotalTests{} tests with raw p-values, BH-adjusted q-values, and Bonferroni thresholds | §4.6 FDR Audit and Bonferroni Correction |
| `fine_sweep_audit.txt` | p-values and enrichment ratios at each bandwidth step in the fine unit sweep | §5.2 Sensitivity Slope at the Canonical Unit |
| `founding_keyword_audit.txt` | Sites classified by founding-era keyword rules (post-hoc Test 5); includes false-positive review | Appendix A.2 (Test 5, post-hoc) |
| `interharmonic_audit.txt` | Harmonic-node spacing checks: distances between successive nodes at the canonical unit | §3.4 Tier Classification |
| `owtrad_audit.txt` | OWTRad trade-route vertices and edges that fall on A+ harmonic nodes; alignment statistics | §7.2 Architectural Stratification and Route-Segment Alignment |
| `religion_keyword_audit.txt` | Site-by-site religion/tradition keyword assignments; counts per tradition | §4.1 Global Enrichment (Test 1, Primary) |
| `site_as_anchor_audit.txt` | Every UNESCO site scored as a hypothetical anchor; Gerizim rank vs. all alternatives | §5.3 Global Anchor Sweep |
| `stupa_geo_audit.txt` | Geographic distribution summary for the UNESCO stupa sub-corpus | §5.6 Pipeline Portability: Wikidata Q180987 Stupa Corpus |
| `stupa_q180987_geo_audit.txt` | Geographic distribution summary for the Wikidata Q180987 stupa corpus | §5.6 Pipeline Portability: Wikidata Q180987 Stupa Corpus |
| `tier_sensitivity_audit.txt` | Tier-threshold sensitivity audit: 8 subsets × 6 tier levels × null-rate sweep 1%–10% | Appendix: Tier-Threshold Sensitivity Audit |
| `tier_sensitivity_part1.csv` | Machine-readable Part 1: primary null rates, all subsets × all tiers | Appendix: Tier-Threshold Sensitivity Audit |
| `tier_sensitivity_part2.csv` | Machine-readable Part 2: A+ null-rate sweep, all subsets | Appendix: Tier-Threshold Sensitivity Audit |

---

## Audit-Document Generators (`tools/`)

These scripts produce the HTML/Markdown audit documents in `supplementary/audit/`.
Run `bash tools/update_all_audits.sh` to regenerate all of them.

| Script | Audit document produced | Paper section |
|--------|-------------------------|---------------|
| `generate_audit_anchor_sweep.py` | Anchor-sweep sensitivity audit | §5.3 Global Anchor Sweep |
| `generate_audit_aplus_sites.py` | A+ site listing audit | §4.1 Global Enrichment (Test 1, Primary) |
| `generate_audit_corridor.py` | Corridor-geometry audit | §5.3 Global Anchor Sweep |
| `generate_audit_dome.py` | Domed-monument corpus audit | §3.1 Corpus Construction / Appendix A.1 |
| `generate_audit_dome_mound.py` | Dome / mound typology audit | §4.3 Hemispherical Mound Evolution (Test 2b) |
| `generate_audit_fdr.py` | FDR correction audit | §4.6 FDR Audit and Bonferroni Correction |
| `generate_audit_fine_sweep.py` | Fine-sweep bandwidth audit | §5.2 Sensitivity Slope at the Canonical Unit |
| `generate_audit_founding.py` | Founding-era classification audit | §3.5 Test Populations / Appendix A.2 |
| `generate_audit_interharmonic.py` | Inter-harmonic spacing audit | §3.4 Tier Classification |
| `generate_audit_owtrad.py` | OWTRad trade-route audit | §7.2 Architectural Stratification and Route-Segment Alignment |
| `generate_audit_religion.py` | Religion / tradition keyword audit | §4.1 Global Enrichment (Test 1, Primary) |
| `generate_audit_site_as_anchor.py` | Site-as-anchor scoring audit | §5.3 Global Anchor Sweep |
| `generate_audit_stupa_geo.py` | Stupa geographic distribution audit | §5.6 Pipeline Portability: Wikidata Q180987 Stupa Corpus |
| `generate_audit_stupa_q180987_geo.py` | Q180987 stupa geographic audit | §5.6 Pipeline Portability: Wikidata Q180987 Stupa Corpus |
| `generate_audit_tier_sensitivity.py` | Tier-threshold sensitivity audit (all subsets × tiers × null-rate sweep) | Appendix: Tier-Threshold Sensitivity Audit |
| `generate_interactive_map.py` | Interactive HTML map (`supplementary/interactive_map/`) | §7.2 / Figure 5 |
| `owtrad_tier_analysis.py` | OWTRad tier-level analysis | §7.2 Architectural Stratification and Route-Segment Alignment |
| `expand_and_convert.py` | Utility: expand Markdown and convert to PDF | — |

---

## Manuscript & Pipeline Scripts (`manuscript/`)

| Script / File | Role | Paper section |
|---------------|------|---------------|
| `generate_figures.py` | Generates all manuscript figures to `manuscript/figures/` | — |
| `reproduce_all_macros.sh` | Master shell script: runs the full analysis pipeline and regenerates all LaTeX macros | §3.7 Simulation Null Models (pipeline) |
| `generated_macros.tex` | **Auto-generated** — do not edit by hand; produced by `emit_sig_macros.py` and related scripts | — |
| `shared_content.tex` | Shared manuscript body included by both journal targets | — |

---

## Data & Library Modules

| Module | Role |
|--------|------|
| `data/unesco_corpus.py` | Loads and filters the UNESCO World Heritage corpus |
| `data/__init__.py` | Data package initialisation |
| `lib/beru.py` | Beru unit conversion and harmonic-node geometry |
| `lib/dome_filter.py` | Keyword-based domed-monument classifier |
| `lib/founding_filter.py` | Founding-era site classifier |
| `lib/landmarks.py` | Landmark / anchor-site definitions |
| `lib/reporting.py` | Shared reporting utilities (print tables, format p-values) |
| `lib/results_store.py` | Read/write interface for `data/store/results.json` |
| `lib/stats.py` | Core statistical functions (Fisher, Cochran–Armitage, permutation) |
| `lib/sweep.py` | Bandwidth-sweep infrastructure |
| `lib/units.py` | Unit definitions and constants |

---

## Test Suite (`tests/`)

```bash
pytest          # runs full test suite
```

| File | Coverage |
|------|----------|
| `tests/test_beru.py` | Beru geometry and unit conversion |
| `tests/test_config.py` | Configuration loading |
| `tests/test_dome_filter.py` | Dome-filter keyword classifier |
| `tests/test_stats.py` | Statistical functions |
| `tests/test_sweep.py` | Sweep infrastructure |

---

*This index is maintained alongside the code. If a script is added or renamed,
update this file in the same commit.*
