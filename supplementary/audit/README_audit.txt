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

5. anchor_sweep_audit.txt
   Script  : tools/generate_audit_anchor_sweep.py
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
   Method  : Classify every Cultural/Mixed UNESCO site under both anchors using
             symmetric self-exclusion.  Mount Gerizim is added as a synthetic
             corpus entry (Tentative List, ref. 5706; not yet inscribed) so that
             both anchors self-exclude exactly one site, leaving N=1011 each.
               Gerizim corpus : N=1011 (1011 inscribed, excl. Gerizim; incl. Jerusalem)
               Jerusalem corpus: N=1011 (1011 inscribed, excl. Jerusalem; incl. Gerizim)
             Lists all A++ / A+ / A sites sorted by deviation for each anchor.
             Identifies sites unique to each anchor vs sites that are A+ under both.
   Results : Gerizim A++=7, A+=56, A=228 | Jerusalem A++=7, A+=57, A=219
             A+ overlap: 38 sites  |  Gerizim-unique: 18  |  Jerusalem-unique: 19
             Gerizim's 56 A+ matches primary Test 1 (N=1011, p=0.0103) exactly.
             Jerusalem's extra A+ hit is Mount Gerizim itself (4.1 km, Tier A+).
   Output  : Summary comparison table, full A-tier site listings for each
             anchor, overlap section, and unique-to-each-anchor breakdowns.

7. fdr.txt
   Script  : tools/generate_audit_fdr.py
   Test    : Multiple comparisons: FDR and family-wise corrections
   Method  : Reads p-values from data/store/results.json (populated by all
             analysis scripts).  Applies Bonferroni, Holm step-down, and
             Benjamini-Hochberg FDR across all confirmatory + exploratory tests.
   Results : Full correction table; counts of tests surviving each threshold.
   Note    : LaTeX macros emitted by fdr_multiple_comparisons.py are NOT
             included in this file; they live in manuscript/generated_macros.tex.

8. fine_sweep_audit.txt
   Script  : tools/generate_audit_fine_sweep.py
   Test    : Fine-resolution unit sensitivity sweep (§5.2)
   Method  : Evaluates A+ hit rates and p-values at fine bandwidth steps
             (0.05–0.20 beru range) around the canonical unit, confirming the
             sensitivity slope reported in the paper.
   Results : p-values, enrichment ratios, and hit counts at each step.
   Output  : Per-step table and summary of the monotonic sensitivity slope.

9. interharmonic_audit.txt
   Script  : tools/generate_audit_interharmonic.py
   Test    : Harmonic-node spacing check (§3.4 Tier Classification)
   Method  : Computes distances between successive harmonic nodes at the
             canonical unit; verifies the 0.1-beru spacing assumption.
   Results : Node coordinates and inter-node distances across the longitude range.
   Output  : Full node listing with spacing checks.

10. owtrad_audit.txt
   Script  : tools/generate_audit_owtrad.py
   Test    : OWTRad trade-route alignment (§7.2)
   Method  : Identifies OWTRad route vertices and edges that fall on A+
             harmonic nodes; computes alignment statistics vs random expectation.
   Results : Vertex/edge counts, A+ alignment rate, enrichment vs null.
   Output  : Matched vertex and edge listings with tier annotations.

11. site_as_anchor_audit.txt
   Script  : tools/generate_audit_site_as_anchor.py
   Test    : Site-as-anchor ranking — each of 1012 sites used as its own anchor
   Method  : For each site in the extended corpus (1011 inscribed + Gerizim),
             count how many of the remaining 1011 sites fall in the A+ band
             relative to that site's longitude.  Rank all 1012 sites by A+ count.
             Directly answers: "How does Gerizim rank compared to every inscribed
             UNESCO site when each is used as its own beru anchor?"
   Results : Gerizim ranks 28/1012 (top 2.7%); Jerusalem ranks 12/1012 (top 1.1%).
             The top sites by A+ are heavily concentrated in the x.18° phase band;
             Gerizim and Jerusalem are the highest-ranked Levantine sites.
   Output  : Distribution summary, focal anchor table, full ranked listing of all
             1012 sites with A++, A+, A counts and longitude.

12. stupa_geo_audit.txt
   Script  : tools/generate_audit_stupa_geo.py
   Test    : Stupa geographic distribution (§5.6)
   Method  : Summarises the geographic spread of the UNESCO stupa sub-corpus;
             checks for regional concentration that could explain harmonic hits.
   Results : Regional counts, A+ rates by region, concentration tests.
   Output  : Region-by-region breakdown with site listings.

13. stupa_q180987_geo_audit.txt
   Script  : tools/generate_audit_stupa_q180987_geo.py
   Test    : Wikidata Q180987 stupa corpus geographic distribution (§5.6)
   Method  : Same as (12) but for the independent Wikidata Q180987 stupa corpus;
             verifies geographic distribution is not a selection artefact.
   Results : Regional counts, A+ rates by region, concentration tests.
   Output  : Region-by-region breakdown with site listings.

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
  python3 tools/generate_audit_anchor_sweep.py
  python3 tools/generate_audit_site_as_anchor.py
  python3 tools/generate_audit_fine_sweep.py
  python3 tools/generate_audit_interharmonic.py
  python3 tools/generate_audit_owtrad.py
  python3 tools/generate_audit_stupa_geo.py
  python3 tools/generate_audit_stupa_q180987_geo.py

Requires: Python ≥3.9, dependencies in requirements.txt, and the UNESCO
corpus XML + extended-description cache (data/store/).

