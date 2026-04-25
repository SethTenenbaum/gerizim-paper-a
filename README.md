# Domed Monuments Cluster at Babylonian Bēru Harmonics: A Longitude Enrichment Test on the UNESCO World Heritage List

**Paper A — Primary UNESCO Analysis** · `v1.3.0`

Seth Tenenbaum · Independent Scholar
ORCID: [0009-0008-5797-2498](https://orcid.org/0009-0008-5797-2498)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19574076.svg)](https://doi.org/10.5281/zenodo.19574076)
> **v1.3.0 release in progress.** A new Zenodo DOI will be assigned on publication of this tag. The badge above resolves to the latest version via the concept DOI.

---

## Abstract

This repository contains the manuscript, analysis code, and data for
"Domed Monuments Cluster at Babylonian Bēru Harmonics: A Longitude
Enrichment Test on the UNESCO World Heritage List" — a statistical study testing whether domed and
spherical monumental heritage sites in the UNESCO World Heritage List are
non-randomly concentrated near integer-multiple longitudes of the Babylonian
*beru* (30° of arc) measured from Mount Gerizim (35.269°E).

The primary exploratory test finds that **domed/spherical UNESCO monuments** cluster
on beru harmonics at rates significantly exceeding the geometric null
(binomial *p* < 0.001), with supporting evidence from cluster asymmetry
(harmonics containing a site within ≤6.7 km host more total sites than
harmonics without one), a temporal gradient
(pre-2000 vs. post-2000 inscription cohorts), morphological evolution
(hemispherical mound → dome), unit-sensitivity sweep, and global robustness
checks.

## Repository Structure

```
gerizim-paper-a/
├── manuscript/
│   ├── primary/                              # Primary UNESCO manuscript
│   │   ├── paper_a_primary_unesco.tex        # LaTeX source
│   │   └── paper_a_primary_unesco.pdf        # Compiled PDF
│   │
│   ├── archaeometry/                         # Archaeometry submission manuscript
│   │   ├── paper_a_archaeometry.tex          # LaTeX source (Wiley/USG template)
│   │   ├── paper_a_archaeometry.pdf          # Compiled PDF
│   │   ├── USG.cls                           # Wiley Archaeometry document class
│   │   ├── NJDapacite.sty                    # Wiley citation style
│   │   ├── NJDnatbib.sty                     # Wiley natbib style
│   │   ├── algorithm.sty                     # Algorithm package
│   │   ├── algorithmicx.sty                  # Algorithmicx package
│   │   ├── amssymb.sty                       # AMS symbols
│   │   ├── appendix.sty                      # Appendix package
│   │   ├── lettersp.sty                      # Letter spacing package
│   │   ├── Wiley_logo.eps                    # Wiley logo (template asset)
│   │   ├── allergy.eps                       # Template asset
│   │   └── images/                           # Template icon assets (ORCID, logos, etc.)
│   │
│   ├── shared_content.tex            # Body content (single source of truth for both manuscripts)
│   ├── generated_macros.tex          # All pipeline-emitted LaTeX macros (single source of truth)
│   ├── reproduce_all_macros.sh       # Shell script to regenerate all macros + figures
│   ├── generate_figures.py           # Generates all 3 manuscript figures
│   └── figures/
│       ├── fig_devhist.{pdf,png}          # Figure 1: Beru deviation histogram
│       ├── fig_temporal.{pdf,png}         # Figure 2: Temporal gradient
│       ├── fig_unitsweep.{pdf,png}        # Figure 3: Unit sensitivity sweep
│       ├── fig_null_c.{pdf,png}           # Figure 4: Null C bootstrap distributions
│       ├── fig_geo_trail.{pdf,png}        # Figure 5: Eurasian corridor tier distribution
│       └── fig_supp_silkroad_ac.{pdf,png} # Supplementary: Silk Road A/C bimodal
│
├── analysis/
│   ├── unesco/                       # Primary UNESCO corpus tests
│   │   ├── spherical_monument_raw_sweep.py        # Test 2: Dome enrichment (primary)
│   │   ├── spherical_monument_test.py             # Test 2x: Dome enrichment (context-validated)
│   │   ├── cluster_asymmetry_test.py              # Test 3: Phase-split cluster asymmetry
│   │   ├── temporal_gradient_test.py              # Test 4: Temporal gradient
│   │   ├── deep_temporal_analysis.py              # Sequential inscription analysis
│   │   ├── origin_sites_test.py                   # Test E: Founding/origin sites
│   │   ├── founding_sites_analysis.py             # Founding sites detailed analysis
│   │   ├── dome_founding_stratification.py        # Dome × founding stratification
│   │   ├── sacred_origin_test.py                  # Sacred origin keyword test
│   │   ├── meta_keyword_test.py                   # Meta-keyword cross-check
│   │   ├── harmonic_density_attractor_test.py     # Harmonic density attractor
│   │   ├── unesco_buddhist_heritage_test.py       # Buddhist heritage keyword test
│   │   ├── tumulus_dome_evolution_raw_sweep.py    # Mound→dome morphological test
│   │   ├── tumulus_dome_evolution_test.py         # Mound→dome evolution (validated)
│   │   ├── unit_sweep_fill.py                     # Unit sensitivity sweep
│   │   ├── bonferroni_correction.py               # Bonferroni family correction
│   │   ├── simulation_null_model.py               # Permutation/bootstrap null models
│   │   ├── spatial_independence_test.py           # Spatial independence correction
│   │   ├── regional_temporal_gradient.py          # Regional temporal gradient
│   │   ├── region_conditioned_permutation.py      # Region-conditioned permutation test
│   │   ├── fdr_multiple_comparisons.py            # FDR multiple comparisons
│   │   ├── leave_one_out_sensitivity.py           # Leave-one-out sensitivity
│   │   ├── dome_leave_k_out.py                    # Dome leave-k-out robustness
│   │   ├── dome_footprint_window_sensitivity.py   # Footprint window sensitivity
│   │   ├── dome_geographic_concentration_test.py  # Geographic concentration test
│   │   ├── stupa_coordinate_perturbation.py       # Stupa coordinate perturbation
│   │   ├── stupa_geographic_concentration_test.py # Stupa geographic concentration
│   │   ├── sensitivity_slope_permutation_test.py  # Sensitivity slope permutation
│   │   ├── sensitivity_slope_specificity_test.py  # Sensitivity slope specificity
│   │   ├── americas_directional_test.py           # Americas directional control
│   │   ├── mound_keyword_context_audit.py         # Mound keyword context audit
│   │   ├── wikidata_q180987_stupa_audit.py        # Wikidata Q180987 stupa corpus audit (Test 6)
│   │   └── verify_x18_periodicity.py              # x.18°E artifact check
│   │
│   ├── global/                       # Global robustness checks
│   │   ├── anchor_uniqueness_audit.py             # Global anchor sweep
│   │   ├── anchor_site_comparison.py              # Anchor site comparison
│   │   ├── peak_geography_audit.py                # Peak geography audit
│   │   ├── dome_periodicity_audit.py              # Dome periodicity / x.18° artifact test
│   │   ├── global_corridor_comparison.py          # Corridor comparison
│   │   ├── corridor_precision_test.py             # Corridor precision test
│   │   ├── geodesic_sensitivity.py                # Geodesic vs. planar sensitivity
│   │   ├── landmark_anchor_ranking.py             # Gerizim vs. other anchors
│   │   ├── x18_periodicity_formal_test.py         # x.18°E formal periodicity test
│   │   ├── x18_max_permutation_test.py            # x.18°E max permutation test
│   │   ├── x18_optimal_band_significance.py       # x.18°E optimal band significance
│   │   └── emit_constants.py                      # Emit pure constants as LaTeX macros
│   │
│   └── americas/                     # Control comparison (directional test scripts only)
│
├── lib/                              # Shared analysis library
│   ├── beru.py               # Beru-unit calculations, tier classification
│   ├── dome_filter.py        # Context-aware dome/stupa keyword matching
│   ├── founding_filter.py    # Founding/origin site classifier
│   ├── stats.py              # Statistical test wrappers
│   ├── units.py              # Unit conversion utilities
│   ├── landmarks.py          # Landmark coordinate definitions
│   ├── reporting.py          # Output formatting, LaTeX macro generation
│   ├── results_store.py      # Persistent results store (JSON)
│   └── sweep.py              # Anchor sweep analysis
│
├── data/
│   ├── unesco_corpus.py              # Shim for backward-compatible imports
│   ├── scripts/
│   │   ├── unesco_corpus.py                # Canonical UNESCO data loader
│   │   ├── fetch_extended.py               # Fetch extended descriptions from UNESCO
│   │   └── fetch_wikidata_q180987.py       # Fetch / validate Wikidata Q180987 stupa corpus
│   └── store/
│       ├── results.json                    # Aggregated pipeline results
│       ├── unesco/
│       │   ├── unesco.xml                  # UNESCO WHC XML export (1,248 sites)
│       │   ├── extended_cache.json         # Cached extended descriptions
│       │   ├── meta_keyword_results.json   # Pre-computed keyword results
│       │   └── wikidata_stupas_q180987.csv # Wikidata Q180987 stupa corpus (229 sites)
│       └── wikidata/                       # (empty — P1435 corpus removed in v1.2.0)
│
├── supplementary/
│   ├── audit/                        # Keyword-classification and tier audit files (reproducible)
│   │   ├── dome_keyword_audit.txt         # Dome/spherical monument sweep (Test 2)
│   │   ├── dome_mound_keyword_audit.txt   # Dome + mound evolution sweep (Test 2b)
│   │   ├── founding_keyword_audit.txt     # Founding/sacred-origin classifier (Test 5)
│   │   ├── aplus_sites_audit.txt          # A/A+/A++ tier site listing (Gerizim & Jerusalem)
│   │   ├── anchor_sweep_audit.txt         # Global anchor sweep ranking
│   │   ├── corridor_audit.txt             # Levantine corridor site breakdown
│   │   ├── fine_sweep_audit.txt           # Fine unit sweep (±1% of 0.10 bēru)
│   │   ├── interharmonic_audit.txt        # C-band / inter-harmonic site listing
│   │   ├── religion_keyword_audit.txt     # Religion keyword audit (A++/A+/A/C tiers)
│   │   ├── site_as_anchor_audit.txt       # Site-as-own-anchor ranking
│   │   ├── stupa_geo_audit.txt            # Stupa geographic concentration audit
│   │   ├── stupa_q180987_geo_audit.txt    # Wikidata Q180987 stupa geographic audit
│   │   ├── fdr.txt                        # FDR multiple-comparisons output
│   │   └── README_audit.txt              # Audit provenance and reproducibility notes
│   └── UNESCO/                       # Archived source materials for anchor citation (ref. 5706)
│       ├── unesco_5706_rendered.html # Archived rendered UNESCO Tentative List page
│       ├── unesco_5706.pdf           # PDF snapshot
│       └── unesco_5706.png           # Screenshot
│
├── tools/                            # Audit generation scripts
│   ├── generate_audit_dome.py               # Regenerate dome_keyword_audit.txt
│   ├── generate_audit_dome_mound.py         # Regenerate dome_mound_keyword_audit.txt
│   ├── generate_audit_founding.py           # Regenerate founding_keyword_audit.txt
│   ├── generate_audit_aplus_sites.py        # Regenerate aplus_sites_audit.txt
│   ├── generate_audit_anchor_sweep.py       # Regenerate anchor_sweep_audit.txt
│   ├── generate_audit_corridor.py           # Regenerate corridor_audit.txt
│   ├── generate_audit_fine_sweep.py         # Regenerate fine_sweep_audit.txt
│   ├── generate_audit_interharmonic.py      # Regenerate interharmonic_audit.txt
│   ├── generate_audit_religion.py           # Regenerate religion_keyword_audit.txt
│   ├── generate_audit_site_as_anchor.py     # Regenerate site_as_anchor_audit.txt
│   ├── generate_audit_stupa_geo.py          # Regenerate stupa_geo_audit.txt
│   ├── generate_audit_stupa_q180987_geo.py  # Regenerate stupa_q180987_geo_audit.txt
│   ├── generate_audit_fdr.py                # Regenerate fdr.txt
│   ├── update_all_audits.sh                 # Regenerate all audit files in one pass
│   ├── generate_interactive_map.py          # Generate interactive HTML map
│   └── md2pdf.sh                            # Markdown → PDF utility
│
├── results/                          # Cached permutation null distributions
│   ├── x18_maxperm_null_A.npy
│   ├── x18_maxperm_null_B.npy
│   └── x18_maxperm_perm_max.npy
│
├── guide/                            # Reference guides
│   ├── statistical_methods_guide.md  # Statistical methods reference
│   └── statistical_methods_guide.pdf
│
├── tests/                            # Unit tests for shared library
│   ├── test_beru.py
│   ├── test_config.py
│   ├── test_dome_filter.py
│   ├── test_stats.py
│   └── test_sweep.py
│
├── keywords.json             # Single source of truth for all keyword lists
├── config.json               # Pipeline parameters, anchors, thresholds
├── conftest.py               # Pytest path configuration
├── pytest.ini                # Pytest settings
├── requirements.txt          # Python dependencies
├── Makefile                  # Reproduce all results with `make all`
├── LICENSE                   # MIT License (code)
├── LICENSE_MANUSCRIPT.txt    # CC BY 4.0 (manuscript text and figures)
├── LICENSE_DATA.txt          # Data licensing notes (UNESCO, Wikidata)
├── LICENSE_NOTES.md          # Detailed licensing breakdown
└── README.md                 # This file
```

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for the full version history.

## Reproducing the Results

### Prerequisites

- Python ≥ 3.10
- LaTeX distribution (for PDF compilation; optional)

### Setup

```bash
git clone https://github.com/SethTenenbaum/gerizim-paper-a.git
cd gerizim-paper-a
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Run all analyses

```bash
make all          # runs tests, all analysis scripts, and regenerates figures
```

Or individually:

```bash
# Unit tests
pytest

# Primary test (Test 2: Dome enrichment)
python3 analysis/unesco/spherical_monument_raw_sweep.py

# Temporal gradient (Test 4)
python3 analysis/unesco/temporal_gradient_test.py

# Bonferroni correction (all tests)
python3 analysis/unesco/bonferroni_correction.py

# Global anchor uniqueness audit
python3 analysis/global/anchor_uniqueness_audit.py

# Dome periodicity / x.18° artifact test
python3 analysis/global/dome_periodicity_audit.py

# Regenerate manuscript figures
python3 manuscript/generate_figures.py

# Compile LaTeX → PDF — primary manuscript
cd manuscript/primary && pdflatex paper_a_primary_unesco.tex

# Compile LaTeX → PDF — archaeometry submission
cd manuscript/archaeometry && pdflatex paper_a_archaeometry.tex
```

### Key outputs

Each analysis script prints results to stdout, including LaTeX `\newcommand`
macros. All macros are collected into `manuscript/generated_macros.tex` by
running the pipeline build script:

```bash
bash manuscript/reproduce_all_macros.sh          # fast: regenerates macros only (default)
bash manuscript/reproduce_all_macros.sh --full   # slow: reruns all analysis scripts first
```

Both manuscripts (`primary/` and `archaeometry/`) import `../generated_macros.tex`
so that **every number in the PDF is pipeline-driven** — no values are
hand-curated. As of v1.1.0 the manuscripts contain 439 macros, all emitted
by analysis scripts, with zero value mismatches between the manuscripts and
the pipeline.

> **Pipeline status:** macros emitted by analysis scripts · all used values in manuscript are pipeline-driven · 0 hardcoded values

## Data Sources

| Dataset | Source | N |
|---------|--------|---|
| UNESCO World Heritage List | [whc.unesco.org/en/list/xml](https://whc.unesco.org/en/list/xml) | 1,248 sites |
| Extended descriptions | Scraped from individual UNESCO site pages | 1,248 entries |
| Wikidata Q180987 stupa corpus | [query.wikidata.org](https://query.wikidata.org) (P31/P279* wd:Q180987, with P625) | 229 sites |

The UNESCO XML is included in `data/store/unesco/unesco.xml`. To refresh
the extended descriptions cache, run `python3 data/scripts/fetch_extended.py`
(requires internet; may take ~20 min due to rate limiting).

## Manuscripts

Both manuscripts share the same body content via `manuscript/shared_content.tex`.
Each `.tex` file is a thin wrapper providing only the template-specific
preamble, title/author block, and abstract. Edits to the paper body go in
`shared_content.tex` and appear in both PDFs.

### Primary UNESCO analysis
- **Wrapper:** `manuscript/primary/paper_a_primary_unesco.tex` (standard `article` class)
- **Compiled PDF:** `manuscript/primary/paper_a_primary_unesco.pdf`

### Archaeometry submission
- **Wrapper:** `manuscript/archaeometry/paper_a_archaeometry.tex` (Wiley USG.cls)
- **Compiled PDF:** `manuscript/archaeometry/paper_a_archaeometry.pdf`
- Template files (`USG.cls`, `*.sty`, `images/`) live alongside the `.tex` in `manuscript/archaeometry/`.

### Shared resources
- **Body content:** `manuscript/shared_content.tex` — single source of truth for all paper content
- **Macros:** `manuscript/generated_macros.tex` — pipeline-emitted numerical values
- **Figures:** `manuscript/figures/` — referenced as `../figures/` from each subdirectory

All statistical results are embedded as LaTeX macros imported from
`generated_macros.tex`, each annotated with the script that produced it.
This ensures the PDF can be compiled without running any code.

## License

© 2026 Seth Tenenbaum.

- **Code:** MIT License (see `LICENSE`)
- **Manuscript text and figures:** CC BY 4.0 (see `LICENSE_MANUSCRIPT.txt`)
- **Data:** UNESCO XML is © UNESCO; redistributed for research under fair use.
  Wikidata extracts are CC0. (see `LICENSE_DATA.txt`)

See `LICENSE_NOTES.md` for the full per-file breakdown.

## Citation

If you use this code or data, please cite:

> Tenenbaum, S. (2026). Domed Monuments Cluster at Babylonian Bēru Harmonics:
> A Longitude Enrichment Test on the UNESCO World Heritage List. *Preprint.*
> https://github.com/SethTenenbaum/gerizim-paper-a

## Contact

Seth Tenenbaum — sethtenenbaum1@gmail.com
