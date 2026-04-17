# License Notes

This repository uses three licenses. The breakdown by file type:

## Code — MIT License (`LICENSE`)

Applies to all `.py`, `.sh`, `.json`, `.ini`, and configuration files:

- `analysis/**/*.py`
- `lib/*.py`
- `data/scripts/*.py`, `data/unesco_corpus.py`
- `tools/*.py`, `tools/*.sh`
- `tests/*.py`
- `conftest.py`, `config.json`, `keywords.json`, `pytest.ini`
- `Makefile`, `requirements.txt`
- `manuscript/generate_figures.py`, `manuscript/reproduce_all_macros.sh`
- `manuscript/generated_macros.tex` (machine-generated output)

## Manuscript Text and Figures — CC BY 4.0 (`LICENSE_MANUSCRIPT.txt`)

Applies to the authored manuscript content and figures:

- `manuscript/primary/paper_a_primary_unesco.tex` and `.pdf`
- `manuscript/archaeometry/paper_a_archaeometry.tex` and `.pdf`
- `manuscript/figures/*.pdf`, `manuscript/figures/*.png`
- `README.md`, `CHANGELOG.md`, `CITATION.cff`
- `guide/*.md`, `guide/*.pdf`
- `supplementary/audit/*.txt`

## Data — Third-Party Licenses (`LICENSE_DATA.txt`)

- **UNESCO World Heritage XML** (`data/store/unesco/unesco.xml`, `extended_cache.json`):
  Copyright UNESCO. Redistributed for non-commercial research under UNESCO syndication terms.
- **Wikidata extracts** (`data/store/unesco/wikidata_stupas_q180987.csv`):
  CC0 1.0 Universal (public domain).
- **UNESCO supplementary materials** (`supplementary/UNESCO/*`):
  Copyright UNESCO. Archived for citation purposes.

## Wiley Template Files — Not Licensed by This Repository

The following files in `manuscript/archaeometry/` are part of the Wiley USG
template and are subject to Wiley's own terms:

- `USG.cls`, `NJDapacite.sty`, `NJDnatbib.sty`
- `algorithm.sty`, `algorithmicx.sty`, `amssymb.sty`, `appendix.sty`, `lettersp.sty`
- `Wiley_logo.eps`, `allergy.eps`, `images/*`
