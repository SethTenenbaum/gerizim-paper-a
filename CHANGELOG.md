# Changelog

### v1.3.0 — 2026-04-25
**DOI:** (this release — Zenodo will assign DOI on publication of this tag)

#### Analysis changes
- **Religion keyword audit expanded:** `tools/generate_audit_religion.py` updated
  to include A++ and C-- counts, enrichment ratios, and Fisher p-values in all
  per-religion summary tables and site listings. Output regenerated at
  `supplementary/audit/religion_keyword_audit.txt`.
- **A+ site audit expanded:** `tools/generate_audit_aplus_sites.py` updated to
  include A++ p-values and significance stars in the summary table, and a
  chain-rule joint probability (P(A) × P(A+|A) × P(A++|A+)) for the nested
  A/A+/A++ tiers using conditional binomial tests for both Gerizim and
  Jerusalem anchors. Output regenerated at `supplementary/audit/aplus_sites_audit.txt`.

#### Supplementary audit files
- **Nine new audit files added** (all generated deterministically from scripts in
  `tools/`; regenerate with `bash tools/update_all_audits.sh`):
  - `aplus_sites_audit.txt` — A/A+/A++ tier listing with joint probabilities
  - `anchor_sweep_audit.txt` — global anchor sweep ranking (36,000 anchors)
  - `corridor_audit.txt` — Levantine corridor site breakdown
  - `fine_sweep_audit.txt` — fine unit sweep (±1% of 0.10 bēru)
  - `interharmonic_audit.txt` — C-band / inter-harmonic site listing
  - `religion_keyword_audit.txt` — religion keyword audit (A++/A+/A/C-- tiers)
  - `site_as_anchor_audit.txt` — site-as-own-anchor ranking
  - `stupa_geo_audit.txt` — stupa geographic concentration audit
  - `stupa_q180987_geo_audit.txt` — Wikidata Q180987 stupa geographic audit
- **`tools/update_all_audits.sh` added** — regenerates all audit files in one pass.
- **`tools/generate_interactive_map.py` added** — generates an interactive HTML
  map of tier-classified sites.

#### Manuscript / figures
- **Two new figures added:**
  - `fig_null_c.{pdf,png}` — Null C bootstrap distributions (Figure 4)
  - `fig_geo_trail.{pdf,png}` — Eurasian corridor tier distribution (Figure 5)
  - `fig_supp_silkroad_ac.{pdf,png}` — supplementary Silk Road A/C bimodal
- **`shared_content.tex` updated** with expanded religion sub-group results,
  Wikidata Q180987 stupa corpus section (Test 6), and updated robustness /
  sensitivity sections referencing the new audit outputs.

#### Housekeeping
- All version references (`CITATION.cff`, `README.md`, both `.tex` files)
  bumped to `v1.3.0`.
- README directory structure updated to reflect current `supplementary/audit/`,
  `tools/`, and `manuscript/figures/` contents.

---

### v1.2.0 — 2026-04-15
**DOI:** (this release — Zenodo will assign DOI automatically on publication)

#### Analysis changes
- **P1435 (Wikidata heritage property) corpus removed:** `data/store/wikidata/p1435_global_control.csv`,
  `data/store/wikidata/p1435_mesoamerica.csv`, `data/scripts/fetch_p1435_global.py`,
  `analysis/global/wikidata_p1435_control_analysis.py`, and
  `analysis/americas/americas_harmonic_depletion_audit.py` deleted from the
  repository. All manuscript references to the P1435 control analysis removed.
  `data/store/wikidata/` is now empty.
- **Ambiguous "mound" keyword removed:** `mound` removed from `keywords.json`
  keyword lists (was causing ambiguous matches in the tumulus→dome morphological
  evolution test). Manuscript prose for Test 2b updated to list only unambiguous
  keywords (`tumulus`, `tumuli`, `barrow`, `barrows`, `kofun`). All affected
  macros regenerated.
- **Missing macro emissions added** to `analysis/unesco/spherical_monument_raw_sweep.py`,
  `cluster_asymmetry_test.py`, `regional_temporal_gradient.py`,
  `dome_geographic_concentration_test.py`, `origin_sites_test.py`,
  `unit_sweep_fill.py`, `tumulus_dome_evolution_raw_sweep.py`,
  `analysis/global/anchor_site_comparison.py`, `anchor_uniqueness_audit.py`,
  `emit_constants.py`, and `x18_periodicity_formal_test.py`.
- **Wikidata Q180987 stupa fetch/validate script added:**
  `data/scripts/fetch_wikidata_q180987.py` — single authoritative source for
  `data/store/unesco/wikidata_stupas_q180987.csv` (229-site global stupa corpus).
  Supports `--validate` (default) and `--fetch` modes with checksum verification.
- **`analysis/unesco/wikidata_q180987_stupa_audit.py` added:** corpus audit
  script for the Q180987 stupa dataset used in Test 6.

#### Manuscript
- **`manuscript/archaeometry/paper_a_archaeometry.tex` updated:** P1435 control
  analysis section removed; Test 2b mound keyword prose updated; all macros
  now reflect the mound-free keyword set.
- **`manuscript/generated_macros.tex` regenerated** with updated values
  (mound-free tumulus/dome evolution macros; P1435 macros removed).

#### Housekeeping
- All version references (`CITATION.cff`, `README.md`, both `.tex` files)
  bumped to `v1.2.0`.
- `data/store/wikidata/` directory retained (empty) to preserve structure.

---

### v1.1.0 — 2026-04-14
**DOI:** [10.5281/zenodo.19574076](https://doi.org/10.5281/zenodo.19574076)

#### Analysis changes
- **Bonferroni primary family reduced from 4 → 3 tests:** Test 4 (temporal
  gradient / Cochran-Armitage) demoted from the primary Bonferroni family to
  exploratory in `config.json`, with an explicit note that it is a descriptive
  trend only (motivated by the *Archaeometry* submission framing).
- **All three manuscript figures regenerated** following the family change;
  figure annotations now read `ZcochranThree` / `pCochranThree` directly from
  `generated_macros.tex` rather than recomputing inline — ensuring the figure
  and manuscript text are always in sync.
- **`generate_figures.py` refactored:** removed the inline Cochran-Armitage
  computation from the figure script; figures now pull all annotated values
  from the canonical macro file.
- **`data/store/results.json` updated** to reflect the revised family
  membership and new macro values.

#### Manuscript
- **Archaeometry submission added:** `manuscript/archaeometry/paper_a_archaeometry.tex`
  — a parallel manuscript formatted for Wiley *Archaeometry* (USG template).
- **Manuscript directory restructured:** `manuscript/primary/` and
  `manuscript/archaeometry/` subdirectories created; each contains only its
  own `.tex`, `.pdf`, and template-specific files.
- **Wiley template files co-located:** `USG.cls`, all `.sty` files,
  `Wiley_logo.*`, `allergy.*`, and `images/` moved into `manuscript/archaeometry/`.
- **Shared resources** (`generated_macros.tex`, `figures/`, `generate_figures.py`,
  `reproduce_all_macros.sh`) remain at `manuscript/` root; both `.tex` files
  reference them via `../`.

#### Housekeeping
- All version references (`CITATION.cff`, `README.md`, both `.tex` files)
  bumped to `v1.1.0`.
- Merged `archaeometry` branch into `master` and tagged `v1.1.0`.

---

### v1.0.9 — 2026-04-13
- **Keyword filter restored:** `mound_positive_context` in `keywords.json`
  reverted to original broader context list (`earthwork`, `earthen`, `platform`,
  `ritual`, `ceremonial`, `sacred`, `archaeological`, `prehistoric`, `ancient`,
  `constructed`, etc.) — capturing spherical earthworks, platform mounds, and
  burial mounds that may have influenced hemispherical dome architecture.
- **Audit files regenerated:** `dome_mound_keyword_audit.txt` (117 included,
  14 A+ sites) and `dome_keyword_audit.txt` (83 included, 11 A+ sites).
- **README updated:** directory structure corrected to reflect current layout
  (`supplementary/audit/`, `supplementary/UNESCO/`, `tools/`, `results/`).

---

### v1.0.8 — 2026-04-13
- **Supplementary materials reorganized:** `supplementary/audit/` holds all
  keyword-classification audit files; `supplementary/UNESCO/` holds UNESCO
  source files (rendered HTML, PDF, PNG); `fetch_unesco_playwright.py` deleted.
- **Three keyword-audit files added** (all generated from deterministic scripts,
  reproducible with `python3 tools/generate_audit_<name>.py`):
  - `dome_keyword_audit.txt` — dome/spherical monument sweep (Test 2)
  - `dome_mound_keyword_audit.txt` — dome + mound evolution sweep (Test 2b)
  - `founding_keyword_audit.txt` — founding/sacred-origin classifier (Test 3)
- **Audit generation scripts added:** `tools/generate_audit_dome.py`,
  `tools/generate_audit_dome_mound.py`, `tools/generate_audit_founding.py`.
- **`README_audit.txt` updated** to document all three audit files and
  describe reproducibility instructions.
- All keyword lists remain in `keywords.json` (single source of truth);
  no changes to analysis logic.

---

### v1.0.7 — 2026-04-13
- **Abstract:** added inline Tier-A+ definition ($\lesssim$6.7 km of a 0.1-beru
  harmonic node); named Mount Gerizim explicitly; corrected sweep-discovery
  disclosure to reflect actual analysis sequence (Gerizim-first, sweep second).
- **§1.1 Circularity protocol:** integrated Gelman & Loken (2014) forking-paths
  citation; reordered and tightened the timestamp/fabrication argument for
  cleaner flow; removed redundant em-dash punctuation.
- **§2.1 Beru unit:** moved observational-capacity paragraphs (meridian transit,
  eclipse timing, gnomon) here from §2.2; removed duplicate closing scope
  statement; single logical arc from unit definition through territorial use
  to measurement capacity.
- **§2.2 Metrology:** replaced Kangle (1965) citation with first-hand source
  Shamasastry (1915) throughout.
- **§1.4 Interpretive framework:** geographic-concentration item now cites the
  actual corridor statistics (\domeEurasianFraction% vs \fullCorpusEurasianFraction%,
  20°–120°E) rather than asserting "Eurasian" without definition.
- **Table 1 caption:** "Confirmatory family" → "Primary test family" to match
  §1.1 framing that all results are exploratory.
- **Founding-origin classifier description:** compressed 9-line paragraph to 6
  lines; updated keyword-source reference from `config.json` to `keywords.json`.
- **`config.json`:** removed duplicate `"keywords"` stanza (~110 lines);
  replaced with `_keywords_moved` redirect. `keywords.json` is now the sole
  source of truth for all keyword content.
- **Bibliography:** added Gelman & Loken (2014) and Shamasastry (1915).

---

### v1.0.6 — 2026-04-11
- **Keywords revised:** replaced generic archaeoastronomy / history-of-astronomy
  terms with journal-targeted vocabulary for *Journal of Archaeological Science*
  (`longitude periodicity`, `monument distribution`, `binomial enrichment test`,
  `heritage site geography`, `archaeometry`, `Levantine metrology`,
  `quantitative archaeology`). Updated in `CITATION.cff` and manuscript.
- **Anchor citation completed:** archived rendered HTML, screenshot, and PDF of
  UNESCO Tentative List ref. 5706 saved to `supplementary/`. Full footnote with
  archived URL added at every `ref. 5706` occurrence in the manuscript
  (`sec:anchor`, `sec:methods`, landmark-ranking footnote).
- **Submission note updated:** TODO for web archiving marked complete in
  manuscript header comment.

---

### v1.0.5 — 2026-04-10
- **Manuscript revision:** abstract restructured to open with the primary
  discovery (global 3° periodic structure); introduction tightened and
  reordered to state the Gerizim hypothesis before the finding.
- **Prose trimming:** anchor section, interpretive framework, x.18°E
  robustness discussion, alternative-explanations section, limitations, and
  conclusion all condensed to remove repeated restatements of the same
  evidence. Overall manuscript reduced by ~230 lines.
- **Clarity:** added explicit framing that no natural geographic mechanism
  predicts monuments clustering at 3° beru intervals from the Levant;
  unit-sensitivity argument linked directly to the circular-phase objection
  in §2.3.
- **Version bump:** `\date` in manuscript, bibentry, and `CITATION.cff`
  updated to v1.0.5.

---

### v1.0.3 — 2026-04-10
- **Pipeline-driven manuscript:** all 439 manuscript macros are now emitted
  by analysis scripts via `manuscript/reproduce_all_macros.sh`; zero
  hand-curated values remain.
- Added `analysis/global/emit_constants.py` to emit pure constants (Gerizim/
  Jerusalem coordinates, DMS, separation, `\NwhcTotal`).
- Extended macro emission in `temporal_gradient_test.py`,
  `deep_temporal_analysis.py`, `tumulus_dome_evolution_raw_sweep.py`,
  `anchor_site_comparison.py`, `anchor_uniqueness_audit.py`, and
  `bonferroni_correction.py`.
- Added DMS and geometry fields to `config.json`.
- `manuscript/generated_macros.tex` is now the single source of truth for all
  numerical values; 493 macros emitted, 439 consumed by the manuscript.

---

### v1.0.0 — 2026-03-xx
- Initial preprint package: manuscript, all analysis scripts, figures,
  data, tests, CITATION.cff, and license files.
