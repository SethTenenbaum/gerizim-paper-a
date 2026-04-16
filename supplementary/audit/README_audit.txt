Gerizim Paper-A — Supplementary Keyword-Classification Audit Files
==================================================================

Generated : Mon Apr 13 22:00:00 UTC 2026
Version   : v1.0.8 (see CITATION.cff and README.md)
Source    : https://github.com/sethtenenbaum/gerizim-paper-a

All keyword lists are defined in a single source of truth:
  keywords.json  (repo root)

No keyword strings are hardcoded in any analysis script.
Each audit file was generated from a deterministic script that can be
re-run at any time to reproduce the results identically.

────────────────────────────────────────────────────────────────────────────────
FILE INDEX
────────────────────────────────────────────────────────────────────────────────

1. dome_keyword_audit.txt
   Script  : tools/generate_audit_dome.py
   Test    : Spherical / Domed Monument sweep (Test 2 in main paper)
   Keywords: stupa, stupas, tholos (unambiguous)
             dome, domed, domes, spherical (context-validated)
   Method  : Ambiguous keywords require the hit sentence to contain ≥1
             architectural-context term AND no geological/natural-feature
             negative-context term (see keywords.json → dome_forms).
             Full extended UNESCO descriptions are searched (not XML only).
   Results : 83 included, 7 context-rejected (FP rate 7.8%), 11 A+ sites
   Output  : Every Cultural/Mixed UNESCO site that matched ≥1 keyword,
             showing which keywords matched, their validation status, and
             the trigger sentence(s).  Sites where every match was
             context-rejected appear in a separate REJECTED section.

2. dome_mound_keyword_audit.txt
   Script  : tools/generate_audit_dome_mound.py
   Test    : Dome + Hemispherical Mound Evolution sweep (Test 2b)
   Keywords: All dome keywords above PLUS:
             tumulus, tumuli, barrow, barrows, kofun (unambiguous mound)
             mound (context-validated: requires archaeological context)
   Method  : Mound keyword "mound" validated by co-occurrence with
             archaeological-context patterns from
             keywords.json → mound_evolution → mound_positive_context.
   Results : 117 included, 15 context-rejected (FP rate 11.4%), 14 A+ sites
             (strict superset of dome audit: 83 dome + 34 mound-only sites)
   Output  : Same format as (1) but covering the broader mound+dome class.

3. founding_keyword_audit.txt
   Script  : tools/generate_audit_founding.py
   Test    : Founding / Sacred-Origin Site classifier (Test 3 / Table 2)
   Categories (from keywords.json → founding_capital / sacred_origin /
               founding_monument / founding_axis / ancient_landscape):
     F  Founding Capital / Seat of Power
     S  Sacred Origin / Birth of a Tradition
     M  Founding Monument / Prototype
     X  Founding Axis / Imperial Infrastructure
     L  Ancient Continuous Landscape
   Method  : Two-stage matching via lib/founding_filter.py:
             (1) unambiguous keywords → accepted unconditionally
             (2) ambiguous keywords → sentence-level context validation
             (positive-context + negative-context pattern pairs per category)
   Results : 355 classified, 28 A+ sites
             F: 192  S: 56  M: 70  L: 82  X: 26
   Output  : All classified Cultural/Mixed sites, sorted by beru tier,
             showing primary category, all matched categories, matched
             keyword reasons, and trigger sentences.

4. religion_keyword_audit.txt
   Script  : tools/generate_audit_religion.py
   Test    : World Religion keyword sub-population test
   Keywords: Per-religion sets from keywords.json → religion_sets
             (Christianity, Islam, Hinduism, Buddhism, Judaism, Shinto,
              Zoroastrianism, Jainism, Sikhism, Taoism, Confucianism)
   Method  : Substring match against UNESCO site name + short description.
             Binomial test (A+ and A-tier) for each religion sub-corpus.
             Judaism small sample: Fisher exact vs corpus A-tier rate.
   Results : Per-religion N, A+%, A%, enrichment, p-values; union stats.
   Output  : Summary table, Judaism Fisher exact detail, per-religion site
             listings sorted by tier, full union site listing with religion tags.

5. corridor_audit.txt
   Script  : tools/generate_audit_corridor.py
   Test    : Gerizim-Lumbini Corridor Precision Test (Table 5 / Section 3.4)
   Method  : 17 harmonics from 0.0 to 1.6 beru; A+ occupancy per harmonic.
             Binomial (17/17 occupancy), Fisher combined probability,
             A++ enrichment vs non-corridor Fisher exact, Mann-Whitney precision.
             Anchor comparison: Gerizim vs Jerusalem, Megiddo, Bethel,
             plus extended notable-anchor table.
   Results : 17/17 corridor harmonics occupied by A+ UNESCO sites.
   Output  : Per-harmonic site table, full A+ catalogue, statistical test
             results, anchor comparison tables.

6. aplus_sites_audit.txt
   Script  : tools/generate_audit_aplus_sites.py
   Test    : Head-to-head A-tier site comparison: Gerizim vs Jerusalem
   Method  : Classify every Cultural/Mixed UNESCO site under both anchors.
             Lists all A++ / A+ / A sites sorted by deviation for each.
             Identifies sites unique to each anchor vs sites that are A+ under both.
   Results : Gerizim A++=7, A+=56, A=228 | Jerusalem A++=8, A+=57, A=218
             A+ overlap: 39 sites  |  Gerizim-unique: 17  |  Jerusalem-unique: 18
   Output  : Summary comparison table, full A-tier site listings for each
             anchor, overlap section, and unique-to-each-anchor breakdowns.

   NOTE — Jerusalem treatment (two versions):
   (a) Primary tests (Tests 1–4): Jerusalem is a corpus member (N=1011).
       It appears as A+ under the Gerizim anchor (4.1 km) and is counted
       among Gerizim's 56 A+ sites.  This audit reflects that treatment:
       N=1011 for both columns, no exclusions applied.
   (b) Anchor-sensitivity sweeps (§Global Anchor Sweep in main paper):
       Sweep A removes Jerusalem entirely (N=1010, Gerizim → 55 A+);
       Sweep B retains Jerusalem with self-exclusion (N=1010,
       Jerusalem → 56 A+, itself not counted).
       The single-site difference (55 vs. 56 A+ in the sweeps, or
       56 vs. 57 A+ in this raw audit) is the Old City of Jerusalem.

7. fdr.txt
   Script  : tools/generate_audit_fdr.py
   Test    : Multiple comparisons: FDR and family-wise corrections
   Method  : Reads p-values from data/store/results.json (populated by all
             analysis scripts).  Applies Bonferroni, Holm step-down, and
             Benjamini-Hochberg FDR across all confirmatory + exploratory tests.
   Results : Full correction table; counts of tests surviving each threshold.
   Note    : LaTeX macros emitted by fdr_multiple_comparisons.py are NOT
             included in this file; they live in manuscript/generated_macros.tex.

7. unesco_site_by_site_audit.txt
   Legacy audit output (generated by analysis/unesco/spherical_monument_test.py
   in an earlier version; superseded by dome_keyword_audit.txt above).
   Retained for archival continuity.

────────────────────────────────────────────────────────────────────────────────
REPRODUCIBILITY
────────────────────────────────────────────────────────────────────────────────

To regenerate all audit files from scratch (quick — uses existing results.json):

  bash tools/update_all_audits.sh

To regenerate everything including re-running all analysis scripts first:

  bash tools/update_all_audits.sh --full

Individual scripts can also be run directly:

  python3 tools/generate_audit_dome.py
  python3 tools/generate_audit_dome_mound.py
  python3 tools/generate_audit_founding.py
  python3 tools/generate_audit_religion.py
  python3 tools/generate_audit_corridor.py
  python3 tools/generate_audit_aplus_sites.py
  python3 tools/generate_audit_fdr.py

Requires: Python ≥3.9, dependencies in requirements.txt, and the UNESCO
corpus XML + extended-description cache (data/store/).

