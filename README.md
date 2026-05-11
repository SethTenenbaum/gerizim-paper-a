# A Geometric-Null Binomial Framework for Harmonic Longitude Proximity Testing in Monument Distributions

**Paper A — Primary UNESCO Analysis** · `v1.7.0`

Seth Tenenbaum · Independent Scholar
ORCID: [0009-0008-5797-2498](https://orcid.org/0009-0008-5797-2498)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20123709.svg)](https://doi.org/10.5281/zenodo.20123709)

---

## Abstract

This repository contains the manuscript, analysis code, and data for
"A Geometric-Null Binomial Framework for Harmonic Longitude Proximity Testing
in Monument Distributions" — a statistical methods paper addressing a
test–hypothesis mismatch in spatial periodicity analysis.

The Rayleigh statistic measures mean phase concentration at period *T*, whereas
a repeated-line grid hypothesis predicts thresholded proximity to any harmonic
longitude — two summaries that can dissociate within the same corpus. This paper
develops the matched test: the geometric-null binomial framework. A period sweep
of the UNESCO dome/stupa corpus (*N* = 90) demonstrates the dissociation
empirically. Applied to three hemispherical-monument corpora at the 3° period
and ~35.25°E Levantine corridor phase, the framework finds consistent
corpus-level proximity enrichment. The dome-excluded two-corpus combination
(Wikidata + OSM) gives χ²(4) significant at Tier A. A max-statistic
anchor-sweep correction confirms the frame-conditional results should not be
treated as frame-free results. Following Silva (2020), all results are treated
as statistical description of spatial structure.

## Repository Structure

```
gerizim-paper-a/
├── manuscript/
│   ├── paper_a_primary.tex           # LaTeX source (JAS Reports target)
│   ├── paper_a_primary.pdf           # Compiled PDF
│   ├── shared_content.tex            # Shared body content
│   ├── generated_macros.tex          # All pipeline-emitted LaTeX macros
│   ├── reproduce_all_macros.sh       # Regenerate all macros + figures
│   ├── generate_figures.py           # Generates all manuscript figures
│   └── figures/
│       ├── fig_geo_trail.{pdf,png}        # Eurasian corridor tier map
│       └── fig_period_sweep.{pdf,png}     # Period sweep dissociation figure
│
├── analysis/
│   ├── unesco/                       # Primary UNESCO corpus tests
│   │   ├── period_sweep_dissociation.py           # Rayleigh vs binomial dissociation sweep
│   │   ├── spherical_monument_raw_sweep.py        # Dome enrichment (primary)
│   │   ├── spherical_monument_test.py             # Dome enrichment (context-validated)
│   │   ├── cluster_asymmetry_test.py              # Phase-split cluster asymmetry
│   │   ├── temporal_gradient_test.py              # Temporal gradient
│   │   ├── deep_temporal_analysis.py              # Sequential inscription analysis
│   │   ├── origin_sites_test.py                   # Founding/origin sites
│   │   ├── founding_sites_analysis.py             # Founding sites detailed analysis
│   │   ├── dome_founding_stratification.py        # Dome × founding stratification
│   │   ├── sacred_origin_test.py                  # Sacred origin keyword test
│   │   ├── meta_keyword_test.py                   # Meta-keyword cross-check
│   │   ├── meta_keyword_audit.py                  # Meta-keyword audit output
│   │   ├── harmonic_density_attractor_test.py     # Harmonic density attractor
│   │   ├── unesco_buddhist_heritage_test.py       # Buddhist heritage keyword test
│   │   ├── tumulus_dome_evolution_raw_sweep.py    # Mound→dome morphological test
│   │   ├── tumulus_dome_evolution_test.py         # Mound→dome evolution (validated)
│   │   ├── unit_sweep_fill.py                     # Unit sensitivity sweep
│   │   ├── unit_sweep_montecarlo.py               # Monte Carlo unit sweep
│   │   ├── tier_logsweep_sensitivity.py           # Log-scale threshold sweep
│   │   ├── bonferroni_correction.py               # Bonferroni family correction
│   │   ├── simulation_null_model.py               # Permutation/bootstrap null models
│   │   ├── spatial_independence_test.py           # Spatial independence correction
│   │   ├── regional_temporal_gradient.py          # Regional temporal gradient
│   │   ├── region_conditioned_permutation.py      # Region-conditioned permutation test
│   │   ├── fdr_multiple_comparisons.py            # FDR multiple comparisons
│   │   ├── fine_sweep_audit.py                    # Fine unit sweep audit
│   │   ├── corpus_exclusion_audit.py              # Corpus exclusion audit
│   │   ├── leave_one_out_sensitivity.py           # Leave-one-out sensitivity
│   │   ├── dome_leave_k_out.py                    # Dome leave-k-out robustness
│   │   ├── evo_leave_k_out.py                     # Evolution corpus leave-k-out robustness
│   │   ├── dome_footprint_window_sensitivity.py   # Footprint window sensitivity
│   │   ├── dome_geographic_concentration_test.py  # Geographic concentration test
│   │   ├── stupa_coordinate_perturbation.py       # Stupa coordinate perturbation
│   │   ├── stupa_geographic_concentration_test.py # Stupa geographic concentration
│   │   ├── keyword_sensitivity_stupa.py           # Keyword-stratified sensitivity (leave-stupa-out)
│   │   ├── sensitivity_slope_permutation_test.py  # Sensitivity slope permutation
│   │   ├── sensitivity_slope_specificity_test.py  # Sensitivity slope specificity
│   │   ├── americas_directional_test.py           # Americas directional control
│   │   ├── mound_keyword_context_audit.py         # Mound keyword context audit
│   │   ├── verify_phase_peak_periodicity.py       # Phase-peak periodicity formal tests
│   │   └── wikidata_q180987_stupa_audit.py        # Wikidata Q180987 stupa corpus audit
│   │
│   ├── global/                       # Global robustness checks
│   │   ├── anchor_uniqueness_audit.py             # Global anchor sweep
│   │   ├── anchor_site_comparison.py              # Anchor site comparison
│   │   ├── peak_geography_audit.py                # Peak geography audit
│   │   ├── dome_periodicity_audit.py              # Dome periodicity audit
│   │   ├── geodesic_sensitivity.py                # Geodesic vs. planar sensitivity
│   │   ├── landmark_anchor_ranking.py             # Gerizim vs. other anchors
│   │   ├── emit_constants.py                      # Emit pure constants as LaTeX macros
│   │   ├── emit_sig_macros.py                     # Emit significance label macros
│   │   └── owtrad_route_alignment.py              # OWTRAD vertex tier + cluster asymmetry tests
│   │
│   └── americas/                     # Control comparison
│       └── mounds_geo_control.py                  # N. American mound type control
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
│   │   ├── fetch_wikidata_q180987.py       # Fetch / validate Wikidata Q180987 stupa corpus
│   │   ├── fetch_osm_stupas.py             # Fetch OSM stupa features via Overpass API
│   │   └── fetch_owtrad_silk_road.py       # Download and parse OWTRAD route datasets
│   └── store/
│       ├── results.json                    # Aggregated pipeline results
│       ├── unesco/
│       │   ├── unesco.xml                  # UNESCO WHC XML export (1,248 sites)
│       │   ├── extended_cache.json         # Cached extended descriptions
│       │   ├── meta_keyword_results.json   # Pre-computed keyword results
│       │   ├── wikidata_stupas_q180987.csv # Wikidata Q180987 stupa corpus (229 sites)
│       │   └── osm_stupas.csv              # OSM stupa corpus (259 sites)
│       ├── silk_road/
│       │   ├── owtrad_nodes.csv            # OWTRAD network nodes (1,674 unique cities)
│       │   ├── owtrad_routes.csv           # OWTRAD route edges (1,946 segments)
│       │   ├── owtrad_tier_report.txt      # Tier analysis report
│       │   ├── maritime_routes.csv         # Curated maritime Silk Road polylines
│       │   └── overland_trunk.csv          # Curated overland trunk route
│       └── wikidata/                       # (empty — P1435 corpus removed in v1.2.0)
│
├── supplementary/
│   ├── audit/                        # Keyword-classification and tier audit files
│   │   ├── dome_keyword_audit.txt         # Dome/spherical monument sweep
│   │   ├── dome_mound_keyword_audit.txt   # Dome + mound evolution sweep
│   │   ├── founding_keyword_audit.txt     # Founding/sacred-origin classifier
│   │   ├── aplus_sites_audit.txt          # A/A+/A++ tier site listing
│   │   ├── anchor_sweep_audit.txt         # Global anchor sweep ranking
│   │   ├── corridor_audit.txt             # Levantine corridor site breakdown
│   │   ├── fine_sweep_audit.txt           # Fine unit sweep
│   │   ├── interharmonic_audit.txt        # C-band / inter-harmonic site listing
│   │   ├── religion_keyword_audit.txt     # Religion keyword audit
│   │   ├── site_as_anchor_audit.txt       # Site-as-own-anchor ranking
│   │   ├── stupa_geo_audit.txt            # Stupa geographic concentration audit
│   │   ├── stupa_q180987_geo_audit.txt    # Wikidata Q180987 stupa geographic audit
│   │   ├── owtrad_audit.txt               # OWTRAD Silk Road harmonic alignment audit
│   │   ├── fdr.txt                        # FDR multiple-comparisons output
│   │   └── README_audit.txt              # Audit provenance and reproducibility notes
│   ├── interactive_map/
│   │   └── silkroad_interactive.html  # Interactive Silk Road tier map (Leaflet)
│   └── UNESCO/                       # Archived source materials for anchor citation
│       ├── unesco_5706_rendered.html
│       ├── unesco_5706.pdf
│       └── unesco_5706.png
│
├── tools/                            # Audit generation scripts
│   ├── generate_audit_dome.py
│   ├── generate_audit_dome_mound.py
│   ├── generate_audit_founding.py
│   ├── generate_audit_aplus_sites.py
│   ├── generate_audit_anchor_sweep.py
│   ├── generate_audit_corridor.py
│   ├── generate_audit_fine_sweep.py
│   ├── generate_audit_interharmonic.py
│   ├── generate_audit_religion.py
│   ├── generate_audit_site_as_anchor.py
│   ├── generate_audit_stupa_geo.py
│   ├── generate_audit_stupa_q180987_geo.py
│   ├── generate_audit_fdr.py
│   ├── generate_audit_owtrad.py
│   ├── update_all_audits.sh
│   ├── generate_interactive_map.py
│   ├── owtrad_tier_analysis.py
│   └── md2pdf.sh
│
├── results/                          # Cached permutation null distributions
│   ├── x18_maxperm_null_A.npy
│   ├── x18_maxperm_null_B.npy
│   └── x18_maxperm_perm_max.npy
│
├── guide/
│   ├── statistical_methods_guide.md
│   └── statistical_methods_guide.pdf
│
├── tests/
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
├── LICENSE_DATA.txt          # Data licensing notes
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

# Period sweep dissociation demonstration (Section 2)
python3 analysis/unesco/period_sweep_dissociation.py

# Primary dome enrichment test
python3 analysis/unesco/spherical_monument_raw_sweep.py

# Bonferroni correction (all tests)
python3 analysis/unesco/bonferroni_correction.py

# Global anchor sweep
python3 analysis/global/anchor_uniqueness_audit.py

# Regenerate manuscript figures
python3 manuscript/generate_figures.py

# Compile LaTeX → PDF
cd manuscript && pdflatex paper_a_primary.tex
```

### Regenerating macros

Every number in the PDF is pipeline-driven via `manuscript/generated_macros.tex`.
No values are hand-curated.

```bash
bash manuscript/reproduce_all_macros.sh          # fast: regenerates macros only
bash manuscript/reproduce_all_macros.sh --full   # slow: reruns all analysis scripts first
```

> **Pipeline status:** all values in the manuscript are emitted by analysis scripts · 0 hardcoded values

## Data Sources

| Dataset | Source | *N* |
|---------|--------|-----|
| UNESCO World Heritage List | [whc.unesco.org/en/list/xml](https://whc.unesco.org/en/list/xml) | 1,248 sites |
| Extended descriptions | Scraped from individual UNESCO site pages | 1,248 entries |
| Wikidata Q180987 stupa corpus | [query.wikidata.org](https://query.wikidata.org) (`instance of: stupa`, with P625) | 229 sites |
| OSM stupa layer | OpenStreetMap via Overpass API (`building/historic/ruins=stupa`) | 259 sites |
| OWTRAD Silk Road network | [ciolek.com/owtrad](http://www.ciolek.com/owtrad.html) (Ciolek 2004) | 1,946 edges, 1,674 nodes |

The UNESCO XML is included in `data/store/unesco/unesco.xml`. To refresh
the extended descriptions cache, run `python3 data/scripts/fetch_extended.py`
(requires internet; may take ~20 min due to rate limiting).

The Wikidata corpus can be refreshed with `python3 data/scripts/fetch_wikidata_q180987.py`.

The OSM corpus can be refreshed with `python3 data/scripts/fetch_osm_stupas.py`
(requires internet; queries the Overpass API).

The OWTRAD data is pre-processed in `data/store/silk_road/`. To re-download
and rebuild from source, run `python3 data/scripts/fetch_owtrad_silk_road.py`
(requires internet; fetches 29 MapInfo ZIP archives from ciolek.com).

OSM data © OpenStreetMap contributors, available under the Open Database
Licence (ODbL 1.0).

## Manuscript

The manuscript targets *Journal of Archaeological Science: Reports* (JAS Reports).

- **Source:** `manuscript/paper_a_primary.tex`
- **Compiled PDF:** `manuscript/paper_a_primary.pdf`
- **Macros:** `manuscript/generated_macros.tex` — pipeline-emitted numerical values, annotated with the script that produced each value

All statistical results are embedded as LaTeX macros imported from
`generated_macros.tex`. The PDF can be compiled without running any code,
using the pre-generated macros committed to the repository.

## License

© 2026 Seth Tenenbaum.

- **Code:** MIT License (see `LICENSE`)
- **Manuscript text and figures:** CC BY 4.0 (see `LICENSE_MANUSCRIPT.txt`)
- **Data:** UNESCO XML is © UNESCO; redistributed for research under fair use.
  Wikidata extracts are CC0. OSM data is ODbL 1.0. (see `LICENSE_DATA.txt`)

See `LICENSE_NOTES.md` for the full per-file breakdown.

## Citation

If you use this code or data, please cite:

> Tenenbaum, S. (2026). A Geometric-Null Binomial Framework for Harmonic
> Longitude Proximity Testing in Monument Distributions.
> https://doi.org/10.5281/zenodo.20123709

## Contact

Seth Tenenbaum — sethtenenbaum1@gmail.com
