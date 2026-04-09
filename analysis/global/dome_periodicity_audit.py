"""
dome_periodicity_audit.py
=========================
Tests whether the dome sub-population's longitude distribution carries
the same inherent 3° periodicity (x.18°E artifact) as the full UNESCO
corpus, and whether that artifact alone can account for the dome A+
enrichment result.

MOTIVATION
----------
The x.18°E periodic artifact (see anchor_uniqueness_audit.py,
peak_geography_audit.py) shows that the full UNESCO corpus already
has 3° periodicity baked into its longitude distribution: every anchor
at phase x.18°E scores ~60 A+, regardless of historical motivation.
A reviewer may ask: does the dome sub-population (N=90) have the same
inherent periodicity? If so, any anchor near the x.18° phase would
produce elevated A+ counts in the dome population too, meaning the
x.18° artifact—not the beru hypothesis—could explain Test 2.

This script answers that question with four explicit tests:

  PART 1 — Full anchor sweep of the dome sub-population.
      Sweep all 36,000 trial anchors against only the N=90 dome sites.
      Record Gerizim's percentile in this dome-specific sweep.
      The x.18° artifact predicts: every x.18°E anchor scores as well
      as Gerizim. The beru hypothesis predicts: Gerizim is in the top tier
      AND the x.18° anchors show no consistent advantage over Gerizim
      (because the dome population is small and the artifact is driven by
      full-corpus density, not dome-specific concentration).

  PART 2 — Phase distribution of the dome population.
      Compute each dome site's longitude modulo 3.0° and test whether
      the phase distribution is uniform (χ² goodness-of-fit, 10 bins of
      0.3° each). Compare the dome phase distribution to the full corpus.
      If the dome sub-population has the same 3° structure as the full
      corpus, a KS test between the two phase distributions will be non-
      significant. If dome sites are independently concentrated at beru
      harmonics, the phase distribution will peak near 0 (mod 3°) and
      the two distributions will differ.

  PART 3 — Correlation of dome A+ count with x.18° proximity.
      For anchors within ±5° of the optimal Gerizim-corridor (30–40°E),
      plot dome A+ as a function of phase (lon mod 3.0°). If the dome
      result is driven by the x.18° artifact, dome A+ should peak at
      2.18° phase just like the full corpus. If it is driven by genuine
      beru alignment, dome A+ should be anchored to Gerizim's specific
      phase (2.269°) and insensitive to the global x.18° optimum.

  PART 4 — Residual enrichment after removing x.18°-phase sites.
      Exclude all dome sites whose longitude (mod 3.0°) is within 0.15°
      of the optimal 2.18° phase (i.e., sites that would be A+ under any
      x.18°E anchor). Re-run the binomial A+ test on the residual.
      If the dome enrichment is entirely explained by the artifact, the
      residual will be non-significant. If genuine beru structure remains,
      the residual will retain significance or a directional elevation.

USAGE
-----
    cd /path/to/gerizim-analysis
    python3 analysis/global/dome_periodicity_audit.py

OUTPUT
------
    Console report + summary table suitable for copying into the
    manuscript's x.18° artifact subsection (§5.2 / §6).

MANUSCRIPT MACROS PRODUCED (to add to paper_a_primary_unesco.tex)
-----------------------------------------------------------------
    \\domePhaseChiP        — chi-sq p for dome phase uniformity
    \\domePhaseKSP         — KS p comparing dome vs full-corpus phases
    \\domeGerizimPctile    — Gerizim's percentile in dome-specific sweep
    \\domeX18AdvantageN    — number of x.18° anchors beating Gerizim in dome sweep
    \\domeResidualN        — dome sites after removing x.18°-proximate sites
    \\domeResidualAp       — A+ count in residual
    \\domeResidualP        — binomial p in residual
    \\domeResidualEnrich   — enrichment in residual
"""

import sys
import numpy as np
from pathlib import Path
from scipy.stats import chi2_contingency, ks_2samp, binomtest, chisquare

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus
from lib.beru import (
    GERIZIM, BERU, TIER_APLUS, TIER_A_MAX,
    P_NULL_AP, P_NULL_A, CONFIG,
)
from lib.dome_filter import FORM_KEYWORDS, FORM_KEYWORD_RES
from lib.stats import significance_label as sig

# ── Constants ─────────────────────────────────────────────────────────────────
OPTIMAL_X18_PHASE = 2.18      # the artifact phase (mod 3.0°)
PHASE_ARTIFACT_WINDOW = 0.15  # ±0.15° around x.18 phase (matches artifact half-width)
HARMONIC_STEP_DEG = 3.0       # 0.1 beru = 3°
GERIZIM_PHASE = GERIZIM % HARMONIC_STEP_DEG  # 35.269 mod 3.0 = 2.269°

# ── Load corpus and build dome sub-population ─────────────────────────────────
corpus = load_corpus()

all_dome_sites = []
all_cultural_lons = []

for site_obj in corpus:
    if site_obj.category == "Natural":
        continue
    if not site_obj.has_coords:
        continue

    all_cultural_lons.append(site_obj.longitude)
    full_text = site_obj.full_text
    matched_kws = [kw for kw in FORM_KEYWORDS if FORM_KEYWORD_RES[kw].search(full_text)]
    if matched_kws:
        arc = abs(site_obj.longitude - GERIZIM)
        beru_val = arc / BERU
        nearest = round(beru_val * 10) / 10
        dev = abs(beru_val - nearest)
        all_dome_sites.append({
            "name":    site_obj.site,
            "lon":     site_obj.longitude,
            "beru":    beru_val,
            "nearest": nearest,
            "dev":     dev,
            "dev_km":  dev * BERU * 111.0,
            "phase":   site_obj.longitude % HARMONIC_STEP_DEG,
            "keywords": matched_kws,
        })

dome_lons = np.array([s["lon"] for s in all_dome_sites])
full_lons  = np.array(all_cultural_lons)
N_dome = len(all_dome_sites)
N_full = len(full_lons)

n_dome_ap = sum(1 for s in all_dome_sites if s["dev"] <= TIER_APLUS)
dome_ap_rate = n_dome_ap / N_dome

# ══════════════════════════════════════════════════════════════════════════════
# PART 1 — Full anchor sweep of the dome sub-population (36,000 trial anchors)
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 80)
print("  DOME PERIODICITY AUDIT")
print(f"  Dome sub-population: N = {N_dome}  (raw keyword sweep, Test 2)")
print(f"  Full corpus:         N = {N_full}")
print(f"  Gerizim: {GERIZIM}°E  (phase mod 3.0° = {GERIZIM_PHASE:.3f}°)")
print(f"  Optimal x.18° phase: {OPTIMAL_X18_PHASE}°")
print("=" * 80)

print("\n── PART 1: Full anchor sweep of dome population ──────────────────────────")
sweep_anchors = np.arange(0.00, 360.00, 0.01)

def count_dome_ap(anchor_lon):
    arcs  = np.abs(dome_lons - anchor_lon)
    arcs  = np.minimum(arcs, 360.0 - arcs)
    bvs   = arcs / BERU
    nears = np.round(bvs * 10) / 10
    devs  = np.abs(bvs - nears)
    return int(np.sum(devs <= TIER_APLUS))

print("  Running dome-specific anchor sweep (36,000 anchors)…")
dome_ap_sweep = np.array([count_dome_ap(a) for a in sweep_anchors])

# Gerizim's score and percentile in this dome-specific sweep
ger_dome_ap    = count_dome_ap(GERIZIM)
dome_pctile    = 100.0 * float(np.mean(dome_ap_sweep <= ger_dome_ap))
pct_above_ger  = 100.0 * float(np.mean(dome_ap_sweep >= ger_dome_ap))

# x.18° anchors in this dome sweep
x18_mask_sw   = np.abs((sweep_anchors % 3.0) - OPTIMAL_X18_PHASE) < 0.006
x18_dome_ap   = dome_ap_sweep[x18_mask_sw]
x18_dome_mean = float(np.mean(x18_dome_ap))
x18_dome_max  = int(np.max(x18_dome_ap))
n_x18_beating_ger = int(np.sum(x18_dome_ap > ger_dome_ap))
n_x18_total   = int(np.sum(x18_mask_sw))

# Global max of dome sweep
dome_global_max = int(np.max(dome_ap_sweep))
n_at_dome_max   = int(np.sum(dome_ap_sweep == dome_global_max))

print(f"\n  Gerizim dome A+:          {ger_dome_ap}  (of {N_dome}; {100*ger_dome_ap/N_dome:.1f}%)")
print(f"  Gerizim percentile:       {dome_pctile:.1f}th (dome-specific sweep)")
print(f"  Anchors scoring ≥ Gerizim: {pct_above_ger:.1f}%")
print(f"  Dome sweep global max:    {dome_global_max}  (at {n_at_dome_max} anchors)")
print(f"\n  x.18°E anchors (N={n_x18_total}):   mean dome A+ = {x18_dome_mean:.2f},  max = {x18_dome_max}")
print(f"  x.18°E anchors beating Gerizim: {n_x18_beating_ger} of {n_x18_total}")
print(f"\n  INTERPRETATION:")
if n_x18_beating_ger == 0:
    print("  → NO x.18° anchor outperforms Gerizim in the dome sub-population.")
    print("    The x.18° artifact does not explain the dome result.")
elif n_x18_beating_ger < n_x18_total / 2:
    print(f"  → Only {n_x18_beating_ger}/{n_x18_total} x.18° anchors beat Gerizim in the dome sweep.")
    print("    The x.18° artifact provides at most partial explanation.")
else:
    print(f"  → {n_x18_beating_ger}/{n_x18_total} x.18° anchors beat Gerizim in the dome sweep.")
    print("    The x.18° artifact may contribute to the dome result — investigate further.")


# ══════════════════════════════════════════════════════════════════════════════
# PART 2 — Phase distribution of dome vs full corpus (mod 3.0°)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── PART 2: Phase distribution (longitude mod 3.0°) ───────────────────────")

N_BINS = 10
phase_edges = np.linspace(0.0, 3.0, N_BINS + 1)
phase_centers = 0.5 * (phase_edges[:-1] + phase_edges[1:])

dome_phases = dome_lons % HARMONIC_STEP_DEG
full_phases  = full_lons % HARMONIC_STEP_DEG

dome_phase_counts, _ = np.histogram(dome_phases, bins=phase_edges)
full_phase_counts, _ = np.histogram(full_phases,  bins=phase_edges)

# χ² uniformity test for dome phases
dome_expected = np.full(N_BINS, N_dome / N_BINS)
chi2_dome_uniform, p_chi2_dome_uniform = chisquare(dome_phase_counts, dome_expected)

# χ² uniformity test for full corpus phases
full_expected = np.full(N_BINS, N_full / N_BINS)
chi2_full_uniform, p_chi2_full_uniform = chisquare(full_phase_counts, full_expected)

# KS test: dome phase distribution vs full corpus phase distribution
ks_stat, p_ks = ks_2samp(dome_phases, full_phases)

# Identify which phase bin contains Gerizim's phase and x.18° phase
ger_phase_bin  = int(GERIZIM_PHASE // (3.0 / N_BINS))
x18_phase_bin  = int(OPTIMAL_X18_PHASE // (3.0 / N_BINS))

print(f"\n  Phase bins (0–3°, {N_BINS} bins of {3.0/N_BINS:.2f}° each):")
print(f"  {'Bin center':>10}  {'Dome N':>7}  {'Dome %':>7}  {'Full N':>7}  {'Full %':>7}  Note")
for i, (ctr, dc, fc) in enumerate(zip(phase_centers, dome_phase_counts, full_phase_counts)):
    note = ""
    if abs(ctr - GERIZIM_PHASE) < 3.0 / N_BINS / 2:
        note += " ← Gerizim phase"
    if abs(ctr - OPTIMAL_X18_PHASE) < 3.0 / N_BINS / 2:
        note += " ← x.18° optimal"
    print(f"  {ctr:>10.2f}°  {dc:>7}  {100*dc/N_dome:>6.1f}%  {fc:>7}  {100*fc/N_full:>6.1f}%  {note}")

print(f"\n  χ² uniformity — dome phases:        χ²={chi2_dome_uniform:.2f}, df={N_BINS-1}, p={p_chi2_dome_uniform:.4f} {sig(p_chi2_dome_uniform)}")
print(f"  χ² uniformity — full corpus phases: χ²={chi2_full_uniform:.2f}, df={N_BINS-1}, p={p_chi2_full_uniform:.4f} {sig(p_chi2_full_uniform)}")
print(f"  KS test (dome vs full):             D={ks_stat:.4f}, p={p_ks:.4f} {sig(p_ks)}")

print(f"\n  INTERPRETATION:")
if p_chi2_dome_uniform > 0.10:
    print("  → Dome phases are UNIFORM (no 3° periodicity in dome longitudes).")
    print("    The x.18° artifact is a property of the full corpus, not the dome sub-population.")
    print("    This directly refutes the hypothesis that dome enrichment is driven by corpus periodicity.")
else:
    print(f"  → Dome phases are non-uniform (χ² p={p_chi2_dome_uniform:.4f}).")
    print("    Investigate which bins drive the excess (see table above).")

if p_ks > 0.10:
    print("  → Dome and full-corpus phase distributions are INDISTINGUISHABLE (KS ns).")
    print("    Both are influenced by the same underlying geographic distribution.")
else:
    print(f"  → Dome and full-corpus phase distributions DIFFER (KS p={p_ks:.4f}).")
    print("    The dome sub-population has a distinct longitude structure.")


# ══════════════════════════════════════════════════════════════════════════════
# PART 3 — Dome A+ vs. phase in the Levant corridor (30–40°E)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── PART 3: Dome A+ vs. anchor phase in the Levant corridor (30–40°E) ──────")

lev_mask   = (sweep_anchors >= 30.0) & (sweep_anchors <= 40.0)
lev_anchors = sweep_anchors[lev_mask]
lev_dome_ap = dome_ap_sweep[lev_mask]
lev_phases  = lev_anchors % HARMONIC_STEP_DEG

# Anchor closest to x.18° optimal in Levant corridor
x18_lev_idx = np.argmin(np.abs(lev_phases - OPTIMAL_X18_PHASE))
x18_lev_lon = float(lev_anchors[x18_lev_idx])
x18_lev_ap  = int(lev_dome_ap[x18_lev_idx])

# Anchor closest to Gerizim's actual phase
ger_lev_idx = np.argmin(np.abs(lev_phases - GERIZIM_PHASE))
ger_lev_ap_check = int(lev_dome_ap[ger_lev_idx])

# Peak dome A+ in Levant corridor and its phase
lev_peak_idx = int(np.argmax(lev_dome_ap))
lev_peak_lon = float(lev_anchors[lev_peak_idx])
lev_peak_ap  = int(lev_dome_ap[lev_peak_idx])
lev_peak_phase = lev_peak_lon % HARMONIC_STEP_DEG

# Correlation between phase and dome A+ in Levant
phase_corr = np.corrcoef(lev_phases, lev_dome_ap)[0, 1]

print(f"\n  Levant corridor anchors ({len(lev_anchors)} trial anchors, 30–40°E):")
print(f"  Gerizim ({GERIZIM}°E, phase={GERIZIM_PHASE:.3f}°): dome A+ = {ger_dome_ap}")
print(f"  x.18° anchor in corridor ({x18_lev_lon:.2f}°E, phase={x18_lev_lon % 3.0:.3f}°): dome A+ = {x18_lev_ap}")
print(f"  Levant peak dome A+: {lev_peak_ap} at {lev_peak_lon:.3f}°E (phase={lev_peak_phase:.3f}°)")
print(f"  Phase–dome-A+ correlation in corridor: r={phase_corr:.3f}")
print(f"\n  INTERPRETATION:")
if x18_lev_ap <= ger_dome_ap:
    print(f"  → The x.18° anchor in the Levant corridor scores ≤ Gerizim in the dome sweep.")
    print(f"    Gerizim's dome result is not explained by phase proximity to the artifact optimum.")
else:
    print(f"  → The x.18° anchor scores {x18_lev_ap - ger_dome_ap} more than Gerizim in the dome sweep.")
    print(f"    Partial x.18° contribution to dome result; see residual test (Part 4).")


# ══════════════════════════════════════════════════════════════════════════════
# PART 4 — Residual enrichment after removing x.18°-proximate dome sites
# ══════════════════════════════════════════════════════════════════════════════
print("\n── PART 4: Residual dome enrichment after removing x.18°-proximate sites ──")

# A dome site is "x.18°-proximate" if its longitude (mod 3.0°) is within
# ±PHASE_ARTIFACT_WINDOW of the optimal artifact phase (2.18°).
# These are sites that would be A+ for ANY x.18°E anchor; they represent
# the artifact's maximum possible contribution to the dome enrichment.

residual_sites = [
    s for s in all_dome_sites
    if abs(s["phase"] - OPTIMAL_X18_PHASE) > PHASE_ARTIFACT_WINDOW
]
artifact_sites = [
    s for s in all_dome_sites
    if abs(s["phase"] - OPTIMAL_X18_PHASE) <= PHASE_ARTIFACT_WINDOW
]

N_residual   = len(residual_sites)
N_artifact   = len(artifact_sites)
res_ap_sites = [s for s in residual_sites if s["dev"] <= TIER_APLUS]
N_res_ap     = len(res_ap_sites)
res_ap_rate  = N_res_ap / N_residual if N_residual > 0 else 0.0

# Binomial test on residual
res_binom    = binomtest(N_res_ap, N_residual, P_NULL_AP, alternative="greater")
p_res        = res_binom.pvalue
res_enrich   = res_ap_rate / P_NULL_AP

# Artifact-proximate sites: are they disproportionately A+?
art_ap_sites  = [s for s in artifact_sites if s["dev"] <= TIER_APLUS]
N_art_ap      = len(art_ap_sites)
art_ap_rate   = N_art_ap / N_artifact if N_artifact > 0 else 0.0
art_binom     = binomtest(N_art_ap, N_artifact, P_NULL_AP, alternative="greater")
p_art         = art_binom.pvalue

# 2×2 Fisher: artifact-proximate vs residual, A+ vs not
contingency = np.array([
    [N_art_ap,            N_artifact - N_art_ap],
    [N_res_ap,            N_residual - N_res_ap],
])
from scipy.stats import fisher_exact
fisher_or, p_fisher_res_vs_art = fisher_exact(contingency, alternative="two-sided")

print(f"\n  Artifact-window definition: phase ± {PHASE_ARTIFACT_WINDOW}° of {OPTIMAL_X18_PHASE}° (x.18°)")
print(f"  Dome sites within artifact window: {N_artifact} ({100*N_artifact/N_dome:.1f}%)")
print(f"  Dome sites outside artifact window (residual): {N_residual} ({100*N_residual/N_dome:.1f}%)")

print(f"\n  A+ rates:")
print(f"  {'Group':<35}  {'N':>5}  {'A+':>5}  {'Rate':>7}  {'Enrich':>7}  {'p (binom)':>12}")
print(f"  {'─'*35}  {'─'*5}  {'─'*5}  {'─'*7}  {'─'*7}  {'─'*12}")
print(f"  {'Full dome (Test 2)':<35}  {N_dome:>5}  {n_dome_ap:>5}  {100*dome_ap_rate:>6.1f}%  {dome_ap_rate/P_NULL_AP:>6.2f}×  p={binomtest(n_dome_ap, N_dome, P_NULL_AP, alternative='greater').pvalue:.4f} {sig(binomtest(n_dome_ap, N_dome, P_NULL_AP, alternative='greater').pvalue)}")
print(f"  {'Artifact-proximate sites':<35}  {N_artifact:>5}  {N_art_ap:>5}  {100*art_ap_rate:>6.1f}%  {art_ap_rate/P_NULL_AP:>6.2f}×  p={p_art:.4f} {sig(p_art)}")
print(f"  {'Residual (non-artifact-proximate)':<35}  {N_residual:>5}  {N_res_ap:>5}  {100*res_ap_rate:>6.1f}%  {res_enrich:>6.2f}×  p={p_res:.4f} {sig(p_res)}")
print(f"\n  Fisher exact (artifact vs residual A+ rates): OR={fisher_or:.2f}, p={p_fisher_res_vs_art:.4f} {sig(p_fisher_res_vs_art)}")

print(f"\n  Artifact-proximate A+ dome sites (if any):")
for s in art_ap_sites:
    print(f"    {s['name']:<45}  lon={s['lon']:.3f}°  phase={s['phase']:.3f}°  δ={s['dev_km']:.1f} km")

print(f"\n  Residual A+ dome sites:")
for s in res_ap_sites:
    print(f"    {s['name']:<45}  lon={s['lon']:.3f}°  phase={s['phase']:.3f}°  δ={s['dev_km']:.1f} km")

print(f"\n  INTERPRETATION:")
if p_res < 0.05:
    print(f"  → Residual dome enrichment is SIGNIFICANT (p={p_res:.4f}).")
    print("    Even after removing all dome sites within the x.18° artifact window,")
    print("    the dome population retains statistically significant A+ enrichment.")
    print("    The x.18° artifact does NOT account for the Test 2 result.")
elif p_res < 0.10:
    print(f"  → Residual dome enrichment is TREND-LEVEL (p={p_res:.4f}).")
    print("    The x.18° artifact accounts for part but not all of the dome result.")
else:
    print(f"  → Residual dome enrichment is NON-SIGNIFICANT (p={p_res:.4f}).")
    if N_res_ap == 0 and N_art_ap > 0:
        print(f"""
  !! IMPORTANT DIAGNOSTIC — ALL DOME A+ SITES ARE IN THE ARTIFACT WINDOW !!
  All {N_art_ap} dome A+ sites have longitudes (mod 3.0°) within ±{PHASE_ARTIFACT_WINDOW}°
  of the x.18° artifact phase ({OPTIMAL_X18_PHASE}°).

  This means the dome A+ sites are concentrated near longitudes of the form
  x.{OPTIMAL_X18_PHASE:.2f}°E (±{PHASE_ARTIFACT_WINDOW}°), i.e., multiples of 3° shifted by ~2.25°.
  This is EXACTLY the longitude pattern that produces the x.18° artifact in
  the full corpus. The dome sub-population's A+ enrichment appears to be
  explained by the same underlying geographic structure as the full corpus —
  dome sites globally happen to be concentrated near this phase.

  CRUCIAL DISTINCTION (two-sided interpretation):
  ─────────────────────────────────────────────────────────────────────────────
  INTERPRETATION A (artifact explains dome result):
    The dome A+ concentration is driven by the x.18° geographic structure,
    not by beru-specific placement. Domed monuments are built in regions
    (South Asia, the Near East, Central Asia, Italy, Japan) that happen to
    cluster near longitudes of the form x.22–2.33°E (mod 3°).
    The Test 2 result would then reflect the same structural periodicity
    as the full corpus, not a morphology-specific signal.

  INTERPRETATION B (beru alignment causes the phase concentration):
    If ancient builders placed domed monuments at beru harmonics from
    ~35°E, then BOTH the dome sites AND many other UNESCO sites would
    cluster at the same phase (x.22–2.33°E mod 3°), which is precisely
    what the x.18° artifact measures. Under this interpretation, the
    artifact is not a confound but a CONSEQUENCE of widespread beru
    alignment. The fact that Gerizim's phase (2.269°) coincides with
    the enriched dome-site phase band is then evidence FOR the hypothesis,
    not against it.

  HOW TO DISTINGUISH THESE:
    1. Compare dome-site phase concentration to non-dome UNESCO sites at
       the same longitudes: if dome enrichment is morphology-specific,
       non-dome sites at the same longitudes should NOT show comparable
       A+ rates.
    2. Check whether the enriched phase (2.23–2.33° mod 3.0°) matches
       Gerizim's phase (2.269°) more than the nominal x.18° optimum (2.18°).
    3. Independent replication: does a non-UNESCO dome corpus show the
       same phase concentration?

  MANUSCRIPT IMPLICATION:
    The original claim that "the x.18° artifact does not explain the dome
    result" must be REVISED. The dome A+ sites are all phase-proximate to
    the artifact window. This is a stronger confound than previously stated.
    However, Interpretation B remains viable and should be presented as
    an open question requiring independent replication.
    The unit-sensitivity result (§5.1) and the morphological specificity
    (dome 3.06× vs full corpus 1.38×) are the remaining independent
    lines of evidence.
  ─────────────────────────────────────────────────────────────────────────────
        """)
    else:
        print("    The x.18° artifact may account for the dome result — investigate.")


# ══════════════════════════════════════════════════════════════════════════════
# PART 5 — Non-dome sites at same longitudes: is dome enrichment morphology-specific?
# This distinguishes Interpretation A (artifact explains dome) from
# Interpretation B (beru alignment causes phase concentration).
# If dome enrichment is morphology-specific, non-dome sites at the same
# longitudes should NOT show the same A+ rate.
# ══════════════════════════════════════════════════════════════════════════════
print("\n── PART 5: Morphological specificity within the artifact-proximate band ───")

dome_lons_set = set(round(s["lon"], 3) for s in all_dome_sites)

# All Cultural/Mixed sites, partitioned into dome vs non-dome
all_sites_data = []
for site_obj in corpus:
    if site_obj.category == "Natural":
        continue
    if not site_obj.has_coords:
        continue
    arc = abs(site_obj.longitude - GERIZIM)
    beru_val = arc / BERU
    nearest = round(beru_val * 10) / 10
    dev = abs(beru_val - nearest)
    phase = site_obj.longitude % HARMONIC_STEP_DEG
    all_sites_data.append({
        "name":    site_obj.site,
        "lon":     site_obj.longitude,
        "dev":     dev,
        "phase":   phase,
        "is_dome": round(site_obj.longitude, 3) in dome_lons_set,
    })

# Artifact-proximate band: phase within ±PHASE_ARTIFACT_WINDOW of OPTIMAL_X18_PHASE
art_band_all   = [s for s in all_sites_data if abs(s["phase"] - OPTIMAL_X18_PHASE) <= PHASE_ARTIFACT_WINDOW]
art_band_dome  = [s for s in art_band_all if s["is_dome"]]
art_band_nondome = [s for s in art_band_all if not s["is_dome"]]

N_band_all     = len(art_band_all)
N_band_dome    = len(art_band_dome)
N_band_nondome = len(art_band_nondome)

band_dome_ap     = sum(1 for s in art_band_dome    if s["dev"] <= TIER_APLUS)
band_nondome_ap  = sum(1 for s in art_band_nondome if s["dev"] <= TIER_APLUS)

band_dome_rate    = band_dome_ap    / N_band_dome    if N_band_dome    > 0 else 0.0
band_nondome_rate = band_nondome_ap / N_band_nondome if N_band_nondome > 0 else 0.0

# Fisher exact: dome vs non-dome within artifact band
cont_band = np.array([
    [band_dome_ap,    N_band_dome    - band_dome_ap],
    [band_nondome_ap, N_band_nondome - band_nondome_ap],
])
fisher_band_or, p_fisher_band = fisher_exact(cont_band, alternative="two-sided")

# Binomial tests
p_band_dome_binom    = binomtest(band_dome_ap,    N_band_dome,    P_NULL_AP, alternative="greater").pvalue
p_band_nondome_binom = binomtest(band_nondome_ap, N_band_nondome, P_NULL_AP, alternative="greater").pvalue

print(f"\n  Within the artifact-proximate band (phase {OPTIMAL_X18_PHASE}° ± {PHASE_ARTIFACT_WINDOW}°):")
print(f"  {'Group':<35}  {'N':>5}  {'A+':>5}  {'Rate':>7}  {'Enrich':>7}  {'p (binom)':>12}")
print(f"  {'─'*35}  {'─'*5}  {'─'*5}  {'─'*7}  {'─'*7}  {'─'*12}")
print(f"  {'Dome sites in band':<35}  {N_band_dome:>5}  {band_dome_ap:>5}  {100*band_dome_rate:>6.1f}%  {band_dome_rate/P_NULL_AP:>6.2f}×  p={p_band_dome_binom:.4f} {sig(p_band_dome_binom)}")
print(f"  {'Non-dome sites in band':<35}  {N_band_nondome:>5}  {band_nondome_ap:>5}  {100*band_nondome_rate:>6.1f}%  {band_nondome_rate/P_NULL_AP:>6.2f}×  p={p_band_nondome_binom:.4f} {sig(p_band_nondome_binom)}")
print(f"\n  Fisher (dome vs non-dome within band): OR={fisher_band_or:.2f}, p={p_fisher_band:.4f} {sig(p_fisher_band)}")

# Also compare outside the band
nonband_dome    = [s for s in all_sites_data if s["is_dome"]     and abs(s["phase"] - OPTIMAL_X18_PHASE) > PHASE_ARTIFACT_WINDOW]
nonband_nondome = [s for s in all_sites_data if not s["is_dome"] and abs(s["phase"] - OPTIMAL_X18_PHASE) > PHASE_ARTIFACT_WINDOW]
nonband_dome_ap    = sum(1 for s in nonband_dome    if s["dev"] <= TIER_APLUS)
nonband_nondome_ap = sum(1 for s in nonband_nondome if s["dev"] <= TIER_APLUS)
nonband_dome_rate    = nonband_dome_ap    / len(nonband_dome)    if nonband_dome    else 0.0
nonband_nondome_rate = nonband_nondome_ap / len(nonband_nondome) if nonband_nondome else 0.0

print(f"\n  Outside the artifact-proximate band:")
print(f"  {'Group':<35}  {'N':>5}  {'A+':>5}  {'Rate':>7}  {'Enrich':>7}")
print(f"  {'─'*35}  {'─'*5}  {'─'*5}  {'─'*7}  {'─'*7}")
print(f"  {'Dome sites':<35}  {len(nonband_dome):>5}  {nonband_dome_ap:>5}  {100*nonband_dome_rate:>6.1f}%  {nonband_dome_rate/P_NULL_AP:>6.2f}×")
print(f"  {'Non-dome sites':<35}  {len(nonband_nondome):>5}  {nonband_nondome_ap:>5}  {100*nonband_nondome_rate:>6.1f}%  {nonband_nondome_rate/P_NULL_AP:>6.2f}×")

print(f"\n  INTERPRETATION:")
if p_fisher_band < 0.05:
    print(f"  → Within the artifact band, dome sites are significantly MORE enriched")
    print(f"    than non-dome sites (OR={fisher_band_or:.2f}, p={p_fisher_band:.4f}).")
    print(f"    This supports MORPHOLOGICAL SPECIFICITY: the dome concentration")
    print(f"    at the x.18° phase is stronger than the general corpus concentration,")
    print(f"    supporting the beru hypothesis over a pure artifact explanation.")
elif band_dome_rate > band_nondome_rate:
    print(f"  → Within the artifact band, dome sites trend more enriched than non-dome")
    print(f"    (OR={fisher_band_or:.2f}, p={p_fisher_band:.4f}, ns) but the sample is small.")
    print(f"    Direction supports morphological specificity; insufficient power to confirm.")
else:
    print(f"  → Within the artifact band, dome and non-dome sites show similar A+ rates.")
    print(f"    Dome enrichment within the band may not be morphology-specific.")

# Summary values for this part
print(f"\n  Key values for manuscript:")
print(f"    Band dome A+ rate:    {100*band_dome_rate:.1f}%  (N={N_band_dome}, A+={band_dome_ap})")
print(f"    Band non-dome A+ rate:{100*band_nondome_rate:.1f}%  (N={N_band_nondome}, A+={band_nondome_ap})")
print(f"    Fisher OR: {fisher_band_or:.2f}, p={p_fisher_band:.4f}")


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY — Manuscript macro values
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("  SUMMARY — COPY THESE MACRO VALUES INTO paper_a_primary_unesco.tex")
print("=" * 80)

print(f"""
% Dome Periodicity Audit (analysis/global/dome_periodicity_audit.py)
% Tests whether x.18°E artifact explains the dome A+ enrichment (Test 2).
\\newcommand{{\\domePhaseChiSq}}{{{chi2_dome_uniform:.2f}}}     % chi-sq for dome phase uniformity
\\newcommand{{\\domePhaseChiP}}{{{p_chi2_dome_uniform:.4f}}}     % p-value, dome phase chi-sq
\\newcommand{{\\domePhaseChiSig}}{{{sig(p_chi2_dome_uniform)}}}        % significance label
\\newcommand{{\\fullPhaseChiP}}{{{p_chi2_full_uniform:.4f}}}     % p-value, full-corpus phase chi-sq
\\newcommand{{\\domePhaseKSD}}{{{ks_stat:.4f}}}      % KS statistic, dome vs full phases
\\newcommand{{\\domePhaseKSP}}{{{p_ks:.4f}}}      % KS p-value, dome vs full phases
\\newcommand{{\\domeGerizimPctile}}{{{dome_pctile:.1f}}}    % Gerizim's percentile in dome sweep
\\newcommand{{\\domeArtMeanAp}}{{{x18_dome_mean:.2f}}}      % mean dome A+ at x.18°E anchors
\\newcommand{{\\domeArtMaxAp}}{{{x18_dome_max}}}         % max dome A+ at any x.18°E anchor
\\newcommand{{\\domeArtBeatGer}}{{{n_x18_beating_ger}}}         % x.18° anchors beating Gerizim in dome sweep
\\newcommand{{\\domeArtTotal}}{{{n_x18_total}}}          % total x.18°E anchors in sweep
\\newcommand{{\\domeArtifactN}}{{{N_artifact}}}          % dome sites in x.18° artifact window
\\newcommand{{\\domeResidualN}}{{{N_residual}}}          % dome sites outside artifact window
\\newcommand{{\\domeResidualAp}}{{{N_res_ap}}}          % A+ in residual
\\newcommand{{\\domeResidualRate}}{{{100*res_ap_rate:.1f}}}     % A+ rate in residual (%)
\\newcommand{{\\domeResidualEnrich}}{{{res_enrich:.2f}}}    % enrichment in residual
\\newcommand{{\\domeResidualP}}{{{p_res:.4f}}}      % binomial p, residual
\\newcommand{{\\domeResidualSig}}{{{sig(p_res)}}}        % significance label, residual
% Part 5 — morphological specificity within artifact band
\\newcommand{{\\domeBandN}}{{{N_band_dome}}}           % dome sites in artifact-proximate band
\\newcommand{{\\domeBandAp}}{{{band_dome_ap}}}          % dome A+ in band
\\newcommand{{\\domeBandRate}}{{{100*band_dome_rate:.1f}}}      % dome A+ rate in band (%)
\\newcommand{{\\nonDomeBandN}}{{{N_band_nondome}}}       % non-dome sites in band
\\newcommand{{\\nonDomeBandAp}}{{{band_nondome_ap}}}     % non-dome A+ in band
\\newcommand{{\\nonDomeBandRate}}{{{100*band_nondome_rate:.1f}}} % non-dome A+ rate in band (%)
\\newcommand{{\\domeBandFisherOR}}{{{fisher_band_or:.2f}}}  % Fisher OR, dome vs non-dome in band
\\newcommand{{\\domeBandFisherP}}{{{p_fisher_band:.4f}}}  % Fisher p, dome vs non-dome in band
\\newcommand{{\\domeBandFisherSig}}{{{sig(p_fisher_band)}}}  % significance label
""")

print("  Done.")
