# License Notes

## Code (analysis scripts, library, data loaders)

MIT License — see `LICENSE`

All Python files in `analysis/`, `lib/`, `data/scripts/`, `manuscript/generate_figures.py`,
and the root configuration files (`config.json`, `conftest.py`, `pytest.ini`) are
released under the MIT License.

## Manuscript text and figures

© 2026 Seth Tenenbaum.
Released under **Creative Commons Attribution 4.0 International (CC BY 4.0)**.
https://creativecommons.org/licenses/by/4.0/

Applies to: `manuscript/paper_a_primary_unesco.tex`, `manuscript/paper_a_primary_unesco.pdf`,
and all files in `manuscript/figures/`.

You are free to share and adapt this material for any purpose, provided you
give appropriate credit, link to the license, and indicate if changes were made.

## Data

| File | Source | License |
|------|--------|---------|
| `data/store/unesco/unesco.xml` | UNESCO World Heritage Centre — https://whc.unesco.org/en/list/xml | © UNESCO; redistributed here for non-commercial research purposes only. |
| `data/store/unesco/extended_cache.json` | Scraped from individual UNESCO site pages | © UNESCO; as above. |
| `data/store/unesco/meta_keyword_results.json` | Derived from UNESCO data | © UNESCO (derived); as above. |
| `data/store/wikidata/p1435_global_control.csv` | Wikidata SPARQL query on `wdt:P1435` | CC0 1.0 Universal |

Users who wish to redistribute this repository should be aware that the UNESCO
XML data is subject to UNESCO's terms of use. See:
https://whc.unesco.org/en/syndication/
