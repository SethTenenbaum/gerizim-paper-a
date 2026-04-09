# Longitude Quantization in the UNESCO World Heritage Corpus: Domes, Stupas, and the Babylonian Beru

**Paper A — Primary UNESCO Analysis**

Seth Tenenbaum · Independent Scholar  
ORCID: [0009-0008-5797-2498](https://orcid.org/0009-0008-5797-2498)

---

## Abstract

This repository contains the manuscript, analysis code, and data for
"Longitude Quantization in the UNESCO World Heritage Corpus: Domes, Stupas,
and the Babylonian Beru" — a statistical study testing whether domed and
spherical monumental heritage sites in the UNESCO World Heritage List are
non-randomly concentrated near integer-multiple longitudes of the Babylonian
*beru* (30° of arc) measured from Mount Gerizim (35.269°E).

The confirmatory test finds that **domed/spherical UNESCO monuments** cluster
on beru harmonics at rates significantly exceeding the geometric null
(binomial *p* < 0.001), with supporting evidence from a temporal gradient
(pre-2000 vs. post-2000 inscription cohorts), morphological evolution
(hemispherical mound → dome), unit-sensitivity sweep, and global robustness
checks.

## Repository Structure

```
gerizim-paper-a/
├── manuscript/
│   ├── paper_a_primary_unesco.tex    # LaTeX source (main manuscript)
│   ├── paper_a_primary_unesco.pdf    # Compiled PDF
│   ├── generate_figures.py           # Generates all 3 manuscript figures
│   └── figures/
│       ├── fig_devhist.{pdf,png}     # Figure 1: Beru deviation histogram
│       ├── fig_temporal.{pdf,png}    # Figure 2: Temporal gradient
│       └── fig_unitsweep.{pdf,png}   # Figure 3: Unit sensitivity sweep
│
├── analysis/
│   ├── unesco/                       # Primary UNESCO corpus tests
│   │   ├── spherical_monument_raw_sweep.py   # Test 2: Dome enrichment (confirmatory)
│   │   ├── spherical_monument_test.py        # Test 2x: Dome enrichment (context-validated)
│   │   ├── cluster_asymmetry_test.py         # Test 3: Phase-split cluster asymmetry
│   │   ├── temporal_gradient_test.py         # Test 4: Temporal gradient
│   │   ├── deep_temporal_analysis.py         # Sequential inscription analysis
│   │   ├── origin_sites_test.py              # Test E: Founding/origin sites
│   │   ├── founding_sites_analysis.py        # Founding sites detailed analysis
│   │   ├── sacred_origin_test.py             # Sacred origin keyword test
│   │   ├── meta_keyword_test.py              # Meta-keyword cross-check
│   │   ├── harmonic_density_attractor_test.py # Harmonic density attractor
│   │   ├── unesco_buddhist_heritage_test.py  # Buddhist heritage keyword test
│   │   ├── tumulus_dome_evolution_raw_sweep.py # Mound→dome morphological test
│   │   ├── unit_sweep_fill.py                # Unit sensitivity sweep
│   │   ├── bonferroni_correction.py          # Bonferroni family correction
│   │   ├── simulation_null_model.py          # Permutation/bootstrap null models
│   │   ├── spatial_independence_test.py      # Spatial independence correction
│   │   ├── regional_temporal_gradient.py     # Regional temporal gradient
│   │   ├── fdr_multiple_comparisons.py       # FDR multiple comparisons
│   │   └── verify_x18_periodicity.py         # x.18°E artifact check
│   │
│   ├── global/                       # Global robustness checks
│   │   ├── anchor_uniqueness_audit.py        # Global anchor sweep
│   │   ├── peak_geography_audit.py           # Peak geography audit
│   │   ├── dome_periodicity_audit.py         # Dome periodicity / x.18° artifact test
│   │   ├── global_corridor_comparison.py     # Corridor comparison
│   │   ├── corridor_precision_test.py        # Corridor precision test
│   │   └── landmark_anchor_ranking.py        # Gerizim vs. other anchors
│   │
│   └── americas/                     # Control comparison
│       └── americas_harmonic_depletion_audit.py  # Americas P1435 depletion control
│
├── lib/                              # Shared analysis library
│   ├── beru.py               # Beru-unit calculations, tier classification
│   ├── dome_filter.py        # Context-aware dome/stupa keyword matching
│   ├── founding_filter.py    # Founding/origin site classifier
│   ├── stats.py              # Statistical test wrappers
│   ├── units.py              # Unit conversion utilities
│   ├── landmarks.py          # Landmark coordinate definitions
│   ├── reporting.py          # Output formatting, LaTeX macro generation
│   └── sweep.py              # Anchor sweep analysis
│
├── data/
│   ├── unesco_corpus.py              # Shim for backward-compatible imports
│   ├── scripts/
│   │   ├── unesco_corpus.py          # Canonical UNESCO data loader
│   │   ├── fetch_extended.py         # Fetch extended descriptions from UNESCO
│   │   └── fetch_p1435_global.py     # Fetch Wikidata P1435 global control
│   └── store/
│       ├── unesco/
│       │   ├── unesco.xml                  # UNESCO WHC XML export (1,248 sites)
│       │   ├── extended_cache.json         # Cached extended descriptions
│       │   └── meta_keyword_results.json   # Pre-computed keyword results
│       └── wikidata/
│           └── p1435_global_control.csv    # Wikidata P1435 global heritage control
│
├── tests/                            # Unit tests for shared library
│   ├── test_beru.py
│   ├── test_config.py
│   ├── test_dome_filter.py
│   ├── test_stats.py
│   └── test_sweep.py
│
├── config.json               # All parameters, anchors, keywords, thresholds
├── conftest.py               # Pytest path configuration
├── pytest.ini                # Pytest settings
├── requirements.txt          # Python dependencies
├── Makefile                  # Reproduce all results with `make all`
├── LICENSE                   # MIT License (code)
└── README.md                 # This file
```

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

# Primary confirmatory test (Test 2: Dome enrichment)
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

# Compile LaTeX → PDF (requires pdflatex)
cd manuscript && pdflatex paper_a_primary_unesco.tex
```

### Key outputs

Each analysis script prints results to stdout, including LaTeX `\newcommand`
macros that are pasted into the manuscript header. The manuscript `.tex` file
contains all macro values inline — no external macro files are needed to
compile the PDF.

## Data Sources

| Dataset | Source | N |
|---------|--------|---|
| UNESCO World Heritage List | [whc.unesco.org/en/list/xml](https://whc.unesco.org/en/list/xml) | 1,248 sites |
| Extended descriptions | Scraped from individual UNESCO site pages | 1,248 entries |
| Wikidata P1435 global control | SPARQL query on `wdt:P1435` | ~161k monuments |

The UNESCO XML is included in `data/store/unesco/unesco.xml`. To refresh
the extended descriptions cache, run `python3 data/scripts/fetch_extended.py`
(requires internet; may take ~20 min due to rate limiting).

## Manuscript

- **LaTeX source:** `manuscript/paper_a_primary_unesco.tex`
- **Compiled PDF:** `manuscript/paper_a_primary_unesco.pdf`
- **Figures:** `manuscript/figures/`

The manuscript is self-contained: all statistical results are embedded as
LaTeX macros in the document preamble (lines ~75–600), each annotated with
the script that produced it. This ensures the PDF can be compiled without
running any code.

## What This Repository Does *Not* Include

This is a focused preprint package for Paper A only. The following are
**excluded** as they belong to separate, in-progress analyses:

- Americas / Mesoamerica regional analysis (Paper B/C)
- Wikidata angular/arc-distance tests
- Meru 3° grid analysis
- Antinode analysis
- Stupa database construction
- Arxiv exploratory notebooks
- Website / visualization code

These will be released in separate repositories when their respective
manuscripts are complete.

## License

- **Code:** MIT License (see `LICENSE`)
- **Manuscript text and figures:** CC BY 4.0
- **Data:** UNESCO XML is © UNESCO; redistributed for research under fair use.
  Wikidata extracts are CC0.

## Citation

If you use this code or data, please cite:

> Tenenbaum, S. (2026). Longitude Quantization in the UNESCO World Heritage
> Corpus: Domes, Stupas, and the Babylonian Beru. *Preprint.*
> https://github.com/SethTenenbaum/gerizim-paper-a

## Contact

Seth Tenenbaum — sethtenenbaum1@gmail.com
