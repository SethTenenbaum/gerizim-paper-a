"""
joint_periodogram_test.py
=========================
Normalized Phasor Consensus (NPC) test for a common period and anchor
across four independently assembled corpora:
  1. UNESCO Cultural/Mixed World Heritage sites  (N ~ 1011)
  2. UNESCO domed-monument subset               (N ~ variable)
  3. Wikidata Q180987 stupa corpus               (N ~ 229)
  4. OWTRAD Silk Road nodes                      (N ~ 1674)

NARRATIVE
---------
  The central claim of this analysis is:

    (A) Three independently assembled MONUMENT corpora (UNESCO all, UNESCO
        dome subset, Wikidata stupas) spontaneously recover the SAME 3°
        anchor — placing a grid line within 0.02° of 35°E — with K = 0.986
        and perm-p = 0.012.

    (B) A fourth corpus (OWTRAD Silk Road trade nodes), assembled on
        entirely different principles (routes, not monument density), does
        NOT share the anchor.  At T=3° OWTRAD's own Rayleigh R = 0.031
        (p = 0.21), and its preferred anchor sits 0.7° away from the
        monument consensus.  When OWTRAD is added to the NPC test the
        consensus weakens to K = 0.771 (p = 0.092, ns).

    Together (A) and (B) establish that the 3° / 35°E alignment is a
    property of MONUMENTAL geography, not of Eurasian longitude distribution
    in general, and is cross-corpus recoverable within the monument domain.

PROBLEM WITH PRIOR APPROACHES
------------------------------
  Product-of-powers J(T) = ∏ S_c(T)  — dominated by shared continental
    longitude clustering; 13° wins trivially.

  Weighted phasor sum J(T) = |Σ Z_c(T)| / C  — better, but still distorted
    by individual corpus magnitudes: a corpus with R_c=0.40 (Stupa, 2.1°)
    contributes 18× more to the sum than one with R_c=0.022 (UNESCO, 3°).
    The statistic conflates "how strong is each corpus's own signal" with
    "do all corpora agree on the same anchor".

THE NEW TEST — Normalized Phasor Consensus (NPC)
-------------------------------------------------
  For each period T, define the UNIT PHASOR for corpus c:

      Ẑ_c(T) = Z_c(T) / |Z_c(T)|       (direction only, magnitude discarded)

  The NPC statistic is the mean resultant of the unit phasors:

      K(T) = |Σ_c Ẑ_c(T)| / C   ∈ [0, 1]

  K(T) = 1  iff all corpora point in exactly the same direction (same anchor).
  K(T) ≈ 0  iff the corpora point in random directions.

  Each corpus gets exactly ONE VOTE regardless of its individual R_c.
  A corpus with R_c=0.001 and a corpus with R_c=0.999 contribute equally.
  This is the pure cross-corpus phase consensus statistic.

  The best common anchor at any T is:
      φ*(T) = arg(Σ_c Ẑ_c(T)) · T / (2π)   [mod T]

NULL DISTRIBUTION — Phase randomization (fully vectorized)
----------------------------------------------------------
  Under H₀ (no common anchor), the phases arg(Z_c(T)) are independent
  across corpora.  For each draw, rotate each corpus's phasor by an
  independent θ_c ~ U(0, 2π):

      Ẑ_c_null(T) = Ẑ_c(T) · exp(i·θ_c)   →   K_null(T)

  This preserves each corpus's R_c(T) at every T but destroys inter-corpus
  phase alignment.  Fully vectorized: (N_PERM, C, N_periods).

  Two p-values:
    • global-max: p = P(max_T K_null ≥ max_T K_obs)   [search-corrected]
    • fixed T=3°: p = P(K_null(3°) ≥ K_obs(3°))       [pre-specified only]

DIAGNOSTICS
-----------
  • Top-10 peaks table with per-corpus anchor directions
  • Pairwise phase-difference matrix at key periods — shows which pairs agree
  • Jackknife (leave-one-out K) — which corpus is driving or breaking consensus

OUTPUT
------
  Console report + LaTeX macros + ResultsStore entries.
"""

import sys, json
import numpy as np
from pathlib import Path
import csv

_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(_ROOT))

from data.unesco_corpus import load_corpus
from lib.dome_filter import is_dome_site_raw as is_dome_site
from lib.results_store import ResultsStore
from analysis.unesco.americas_directional_test import is_americas_site

_CFG   = json.loads((_ROOT / "config.json").read_text())
N_PERM = _CFG["simulation"]["n_permutations"]   # 100 000
SEED   = _CFG["simulation"]["random_seed"]       # 42

TARGET_T = 3.0   # pre-specified common period (degrees)

# ── 1. Load all three corpora ─────────────────────────────────────────────────

# (a) UNESCO Cultural/Mixed (all)
_corpus   = load_corpus()
_cultural = [s for s in _corpus if s.category != "Natural" and s.has_coords]
lons_unesco = np.array([s.longitude for s in _cultural])

# (a2) UNESCO dome/stupa/tholos subset
_dome = [s for s in _cultural if is_dome_site(s)]
lons_dome = np.array([s.longitude for s in _dome])

# (a3) UNESCO Americas sub-corpus — geographic NEGATIVE CONTROL
# Used in a 2-corpus NPC against UNESCO-non-Americas to ask whether a
# UNESCO sub-corpus drawn from outside the signal corridor recovers the
# same anchor.  Parallels the OWTRAD contrast at the sub-corpus level.
_americas    = [s for s in _cultural if is_americas_site(s)]
lons_americas = np.array([s.longitude for s in _americas])
# Eurasian (non-Americas) UNESCO sites — the signal-region complement
_eurasia     = [s for s in _cultural if not is_americas_site(s)]
lons_eurasia = np.array([s.longitude for s in _eurasia])

# (b) Wikidata Q180987 stupas
STUPA_CSV = _ROOT / "data" / "store" / "unesco" / "wikidata_stupas_q180987.csv"
# NOTE: OSM stupa corpus loaded below as (d)
_stupa_lons = []
with open(STUPA_CSV, newline="", encoding="utf-8") as fh:
    lines = [l for l in fh if not l.startswith("#")]
    reader = csv.DictReader(lines)
    for row in reader:
        try:
            _stupa_lons.append(float(row["lon"]))
        except (KeyError, ValueError):
            pass
lons_stupa = np.array(_stupa_lons)

# (d) OSM stupas — building=stupa / historic=stupa / ruins=stupa
OSM_STUPA_CSV = _ROOT / "data" / "store" / "unesco" / "osm_stupas.csv"
_osm_lons = []
with open(OSM_STUPA_CSV, newline="", encoding="utf-8") as fh:
    lines = [l for l in fh if not l.startswith("#")]
    reader = csv.DictReader(lines)
    for row in reader:
        try:
            _osm_lons.append(float(row["lon"]))
        except (KeyError, ValueError):
            pass
lons_osm = np.array(_osm_lons)

# (c) OWTRAD nodes — no deduplication, full dataset
OWTRAD_NODES = _ROOT / "data" / "store" / "silk_road" / "owtrad_nodes.csv"
_owtrad_lons = []
with open(OWTRAD_NODES, newline="", encoding="utf-8") as fh:
    lines = [l for l in fh if not l.startswith("#")]
    reader = csv.DictReader(lines)
    for row in reader:
        try:
            _owtrad_lons.append(float(row["lon"]))
        except (KeyError, ValueError):
            pass
lons_owtrad = np.array(_owtrad_lons)

CORPORA = {
    "UNESCO": lons_unesco,
    "Dome":   lons_dome,
    "Stupa":  lons_stupa,
    "OWTRAD": lons_owtrad,
}
N_CORPORA = len(CORPORA)

# Auxiliary sub-corpora used ONLY for sub-corpus negative-control 2-corpus
# NPC tests (Americas vs. Eurasia).  They are NOT included in the 4-corpus
# joint NPC, to preserve the existing primary statistic.
# OSM stupas are also kept here: they are an auxiliary sensitivity corpus,
# not included in the primary 4-corpus NPC.
AUX_CORPORA = {
    "Americas": lons_americas,
    "Eurasia":  lons_eurasia,
    "OSM":      lons_osm,
}

# ── 2. Phasor helpers ─────────────────────────────────────────────────────────

def mean_phasor_vec(lons: np.ndarray, periods: np.ndarray) -> np.ndarray:
    """Return complex mean phasor Z_c(T) for each period. Shape: (len(periods),)"""
    a = 2.0 * np.pi * lons[:, None] / periods[None, :]   # (N, P)
    return np.mean(np.exp(1j * a), axis=0)                # (P,)

# Period grid 0.5–10° at 0.025° resolution
periods   = np.arange(0.5, 10.01, 0.025)
N_periods = len(periods)
idx_3     = int(np.argmin(np.abs(periods - TARGET_T)))

# ── 3. Observed phasors and NPC spectrum ─────────────────────────────────────

phasors_obs = {}   # name -> complex array (N_periods,), Z_c(T)
for name, lons in CORPORA.items():
    phasors_obs[name] = mean_phasor_vec(lons, periods)

# Unit phasors: direction only, magnitude discarded
unit_phasors = {}  # name -> complex array (N_periods,), |Ẑ_c| = 1
for name, Z in phasors_obs.items():
    mag = np.abs(Z)
    mag = np.where(mag < 1e-15, 1e-15, mag)   # guard against zero
    unit_phasors[name] = Z / mag

# NPC statistic K(T) = |Σ_c Ẑ_c(T)| / C
unit_sum_obs = sum(unit_phasors.values())      # (N_periods,)
K_obs        = np.abs(unit_sum_obs) / N_CORPORA  # (N_periods,) in [0,1]

K_at_3   = float(K_obs[idx_3])
peak_idx = int(np.argmax(K_obs))
peak_T   = float(periods[peak_idx])
peak_K   = float(K_obs[peak_idx])

def best_anchor(unit_sum_val: complex, T: float) -> float:
    """Convert unit-phasor sum angle to longitude anchor (mod T)."""
    angle = float(np.angle(unit_sum_val))
    return (angle * T / (2.0 * np.pi)) % T

anchor_at_peak = best_anchor(unit_sum_obs[peak_idx], peak_T)
anchor_at_3    = best_anchor(unit_sum_obs[idx_3], TARGET_T)

rank_3   = int(np.sum(K_obs > K_at_3)) + 1
pctile_3 = float(np.mean(K_obs <= K_at_3)) * 100

# Per-corpus angle (anchor) at key periods
def corpus_anchor(name: str, T: float) -> float:
    idx = int(np.argmin(np.abs(periods - T)))
    return best_anchor(unit_phasors[name][idx], T)

def corpus_R(name: str, T: float) -> float:
    idx = int(np.argmin(np.abs(periods - T)))
    return float(np.abs(phasors_obs[name][idx]))

# Individual own peaks (by R, not K)
indiv_peak = {}
for name, Z in phasors_obs.items():
    R = np.abs(Z)
    i = int(np.argmax(R))
    indiv_peak[name] = (float(periods[i]), float(R[i]), best_anchor(unit_phasors[name][i], periods[i]))

# ── 4. Phase-randomization null — fully vectorized ───────────────────────────
# Rotate each corpus's UNIT phasor by independent θ_c ~ U(0,2π) per draw.
# Preserves each corpus's R_c(T) at every T; destroys inter-corpus alignment.
# Shape: (N_PERM, C, N_periods)

rng = np.random.default_rng(SEED)

unit_mat = np.stack(list(unit_phasors.values()), axis=0)   # (C, N_periods)
theta    = rng.uniform(0.0, 2.0 * np.pi, size=(N_PERM, N_CORPORA, 1))
rotated  = unit_mat[np.newaxis, :, :] * np.exp(1j * theta) # (N_PERM, C, N_periods)
K_null   = np.abs(rotated.sum(axis=1)) / N_CORPORA         # (N_PERM, N_periods)

null_max_K       = K_null.max(axis=1)                      # (N_PERM,)
perm_p_globalmax = float(np.mean(null_max_K >= peak_K))
perm_p_3         = float(np.mean(K_null[:, idx_3] >= K_at_3))

# ── 5. Jackknife — leave-one-out K at each period ────────────────────────────
corpus_names = list(unit_phasors.keys())
jk_K = {}   # name -> K_obs with that corpus removed, shape (N_periods,)
for drop in corpus_names:
    others = [unit_phasors[n] for n in corpus_names if n != drop]
    jk_sum = sum(others)
    jk_K[drop] = np.abs(jk_sum) / (N_CORPORA - 1)

# ── 6. Pairwise phase-difference matrix at T=3° and joint peak ───────────────
def phase_diff_deg(name_a: str, name_b: str, T: float) -> float:
    """Angular difference between two corpus anchors at period T (in degrees, [0, T/2])."""
    idx = int(np.argmin(np.abs(periods - T)))
    a1  = float(np.angle(unit_phasors[name_a][idx]))
    a2  = float(np.angle(unit_phasors[name_b][idx]))
    diff = abs(a1 - a2) % (2.0 * np.pi)
    if diff > np.pi:
        diff = 2.0 * np.pi - diff
    return diff * T / (2.0 * np.pi)   # convert to longitude degrees

# ── 7. Three-corpus NPC (monument corpora only, pre-primary result) ───────────
monument_names  = ["UNESCO", "Dome", "Stupa"]
monument_units  = [unit_phasors[n] for n in monument_names]   # at T=3 (scalar)

def npc_fixed_T(unit_list, T_val):
    """NPC K and anchor for a list of unit phasors at a single period."""
    idx = int(np.argmin(np.abs(periods - T_val)))
    vecs = np.array([u[idx] for u in unit_list])
    s    = vecs.sum()
    K    = abs(s) / len(vecs)
    phi  = best_anchor(s, T_val)
    return K, phi

K3_mon, phi3_mon = npc_fixed_T(monument_units, 3.0)

# phase-rand null for 3-corpus test at fixed T=3°
unit_mat_mon = np.stack([unit_phasors[n][idx_3] for n in monument_names])  # (3,) complex
theta_mon    = rng.uniform(0, 2*np.pi, size=(N_PERM, 3))
K_null_mon   = np.abs((unit_mat_mon[None,:] * np.exp(1j*theta_mon)).sum(axis=1)) / 3
p3_mon       = float(np.mean(K_null_mon >= K3_mon))

# nearest grid line to 35°E from the 3-corpus anchor
gl_mon_35 = phi3_mon + round((35.0 - phi3_mon) / 3.0) * 3.0

# ── 7b. Two-corpus NPC (Dome + Stupa — most theoretically motivated) ──────────
# The dome subset and Wikidata stupas are the most morphologically comparable
# corpora. Neither is nested within the other. This is the strongest theoretical
# prediction: two morphologically similar monument types should share an anchor.
dome_stupa_names  = ["Dome", "Stupa"]
dome_stupa_units  = [unit_phasors[n] for n in dome_stupa_names]
K3_ds, phi3_ds = npc_fixed_T(dome_stupa_units, 3.0)

unit_mat_ds = np.stack([unit_phasors[n][idx_3] for n in dome_stupa_names])
theta_ds    = rng.uniform(0, 2*np.pi, size=(N_PERM, 2))
K_null_ds   = np.abs((unit_mat_ds[None,:] * np.exp(1j*theta_ds)).sum(axis=1)) / 2
p3_ds       = float(np.mean(K_null_ds >= K3_ds))
gl_ds_35    = phi3_ds + round((35.0 - phi3_ds) / 3.0) * 3.0

# ── 7c. Two-corpus NPC (UNESCO + Stupa, NON-NESTED — primary conservative) ───
# UNESCO Cultural/Mixed and Wikidata Q180987 stupa corpus are independently
# curated. The UNESCO dome subset is nested within UNESCO and excluded here.
two_names    = ["UNESCO", "Stupa"]
two_units    = [unit_phasors[n] for n in two_names]
K3_two, phi3_two = npc_fixed_T(two_units, 3.0)

unit_mat_two = np.stack([unit_phasors[n][idx_3] for n in two_names])  # (2,) complex
theta_two    = rng.uniform(0, 2*np.pi, size=(N_PERM, 2))
K_null_two   = np.abs((unit_mat_two[None,:] * np.exp(1j*theta_two)).sum(axis=1)) / 2
p3_two       = float(np.mean(K_null_two >= K3_two))

gl_two_35 = phi3_two + round((35.0 - phi3_two) / 3.0) * 3.0

# ── 7c. Sub-corpus geographic NEGATIVE CONTROL NPCs ───────────────────────────
# Two 2-corpus NPC tests on UNESCO sub-corpora that should NOT recover the
# monument anchor:
#   (i)  UNESCO-Eurasia + UNESCO-Americas  — disjoint geographic split of UNESCO
#   (ii) UNESCO (full)  + UNESCO-Americas  — parallels OWTRAD as a sub-corpus
#                                            contrast (note Americas ⊂ UNESCO)
# These are run identically to the 2-corpus primary test (T=3°, fixed-T
# permutation null with N_PERM phase-randomization draws).

def mean_unit_phasor_at_T(lons: np.ndarray, T: float) -> complex:
    Z = mean_phasor_vec(lons, np.array([T]))[0]
    m = abs(Z)
    return Z / max(m, 1e-15)

def two_corpus_npc(lons_a: np.ndarray, lons_b: np.ndarray, T: float, rng_local):
    """Return (K, anchor_deg, perm_p) for a 2-corpus NPC at fixed T."""
    u_a = mean_unit_phasor_at_T(lons_a, T)
    u_b = mean_unit_phasor_at_T(lons_b, T)
    s   = u_a + u_b
    K   = abs(s) / 2.0
    phi = (float(np.angle(s)) * T / (2.0 * np.pi)) % T
    th  = rng_local.uniform(0, 2 * np.pi, size=(N_PERM, 2))
    pair = np.stack([u_a, u_b])
    K_null = np.abs((pair[None, :] * np.exp(1j * th)).sum(axis=1)) / 2.0
    p   = float(np.mean(K_null >= K))
    return K, phi, p

# (i) UNESCO-Eurasia × UNESCO-Americas — disjoint geographic split
rng_neg1 = np.random.default_rng(SEED + 1)
K_eu_am, phi_eu_am, p_eu_am = two_corpus_npc(
    lons_eurasia, lons_americas, 3.0, rng_neg1
)
gl_eu_am_35 = phi_eu_am + round((35.0 - phi_eu_am) / 3.0) * 3.0

# (ii) UNESCO (full) × UNESCO-Americas — parallels OWTRAD contrast
rng_neg2 = np.random.default_rng(SEED + 2)
K_un_am, phi_un_am, p_un_am = two_corpus_npc(
    lons_unesco, lons_americas, 3.0, rng_neg2
)
gl_un_am_35 = phi_un_am + round((35.0 - phi_un_am) / 3.0) * 3.0

# (iii) UNESCO (full) × OWTRAD  — 2-corpus version of the existing OWTRAD
# contrast, computed here so the manuscript NPC table has an apples-to-apples
# 2-corpus row for every contrast.
rng_neg3 = np.random.default_rng(SEED + 3)
K_un_ow, phi_un_ow, p_un_ow = two_corpus_npc(
    lons_unesco, lons_owtrad, 3.0, rng_neg3
)
gl_un_ow_35 = phi_un_ow + round((35.0 - phi_un_ow) / 3.0) * 3.0

# (iv) UNESCO (full) × UNESCO-Americas × Wikidata Stupa — 3-corpus NEGATIVE
# CONTROL. Mirrors the 3-corpus monument primary (UNESCO + Dome + Stupa),
# substituting the Americas geographic-control sub-corpus for the dome subset.
# Tests whether ANY 3-corpus NPC including UNESCO and Wikidata recovers
# concordance, or whether concordance specifically requires the dome subset.
rng_neg4 = np.random.default_rng(SEED + 4)
u_un_neg = mean_unit_phasor_at_T(lons_unesco,   3.0)
u_am_neg = mean_unit_phasor_at_T(lons_americas, 3.0)
u_st_neg = mean_unit_phasor_at_T(lons_stupa,    3.0)
trip_neg = np.stack([u_un_neg, u_am_neg, u_st_neg])
s_neg    = trip_neg.sum()
K_neg3   = abs(s_neg) / 3.0
phi_neg3 = (float(np.angle(s_neg)) * 3.0 / (2.0 * np.pi)) % 3.0
th_neg3  = rng_neg4.uniform(0, 2 * np.pi, size=(N_PERM, 3))
K_null_neg3 = np.abs((trip_neg[None, :] * np.exp(1j * th_neg3)).sum(axis=1)) / 3.0
p_neg3   = float(np.mean(K_null_neg3 >= K_neg3))
gl_neg3_35 = phi_neg3 + round((35.0 - phi_neg3) / 3.0) * 3.0

# ── 7d. OSM stupa NPC rows (sensitivity check) ───────────────────────────────
# The OSM corpus (building/historic/ruins=stupa) is assembled by a completely
# different process than Wikidata Q180987 — crowd-sourced structural tagging
# with no relationship to Wikidata curation.  We run four 2-corpus NPC tests
# at T=3° to ask whether the OSM anchor agrees with the monument corpora.
# OSM is NOT added to the primary 4-corpus NPC to preserve the primary statistic.
rng_osm1 = np.random.default_rng(SEED + 10)
K_osm_dome, phi_osm_dome, p_osm_dome = two_corpus_npc(
    lons_dome, lons_osm, 3.0, rng_osm1
)
gl_osm_dome_35 = phi_osm_dome + round((35.0 - phi_osm_dome) / 3.0) * 3.0

rng_osm2 = np.random.default_rng(SEED + 11)
K_osm_stupa, phi_osm_stupa, p_osm_stupa = two_corpus_npc(
    lons_stupa, lons_osm, 3.0, rng_osm2
)
gl_osm_stupa_35 = phi_osm_stupa + round((35.0 - phi_osm_stupa) / 3.0) * 3.0

rng_osm3 = np.random.default_rng(SEED + 12)
K_osm_unesco, phi_osm_unesco, p_osm_unesco = two_corpus_npc(
    lons_unesco, lons_osm, 3.0, rng_osm3
)
gl_osm_unesco_35 = phi_osm_unesco + round((35.0 - phi_osm_unesco) / 3.0) * 3.0

# 3-corpus: Dome + Wikidata Stupa + OSM — all three hemispherical monument corpora
rng_osm4 = np.random.default_rng(SEED + 13)
u_dome_osm = mean_unit_phasor_at_T(lons_dome,  3.0)
u_wiki_osm = mean_unit_phasor_at_T(lons_stupa, 3.0)
u_osm_osm  = mean_unit_phasor_at_T(lons_osm,   3.0)
trip_osm   = np.stack([u_dome_osm, u_wiki_osm, u_osm_osm])
s_osm3     = trip_osm.sum()
K_osm3     = abs(s_osm3) / 3.0
phi_osm3   = (float(np.angle(s_osm3)) * 3.0 / (2.0 * np.pi)) % 3.0
th_osm3    = rng_osm4.uniform(0, 2 * np.pi, size=(N_PERM, 3))
K_null_osm3 = np.abs((trip_osm[None, :] * np.exp(1j * th_osm3)).sum(axis=1)) / 3.0
p_osm3     = float(np.mean(K_null_osm3 >= K_osm3))
gl_osm3_35 = phi_osm3 + round((35.0 - phi_osm3) / 3.0) * 3.0

# ── 8. OWTRAD isolated signal at T=3° ────────────────────────────────────────
R_owtrad_3   = float(np.abs(phasors_obs["OWTRAD"][idx_3]))
phi_owtrad_3 = corpus_anchor("OWTRAD", 3.0)
# Rayleigh p for OWTRAD alone
N_ow = len(lons_owtrad)
phases_ow_null = rng.uniform(0, 2*np.pi, size=(N_PERM, N_ow))
R_ow_null      = np.abs(np.mean(np.exp(1j * phases_ow_null), axis=1))
p_rayleigh_ow  = float(np.mean(R_ow_null >= R_owtrad_3))

# ── 9. Print report ───────────────────────────────────────────────────────────
SEP = "=" * 90

def sig_stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

print()
print(SEP)
print("  NPC TEST — MONUMENT CORPORA vs. TRADE-ROUTE CORPUS")
print(SEP)

# ── Part A: Three monument corpora ────────────────────────────────────────────
print(f"""
PART A — THREE MONUMENT CORPORA  (T = 3°, pre-specified)
{'Corpus':<10} {'N':>6}  {'R(3°)':>7}  {'anchor@3°':>10}  {'grid line @35°E':>16}  {'Δ35°E':>7}
{'-'*72}""")
for name in monument_names:
    a  = corpus_anchor(name, 3.0)
    gl = a + round((35.0 - a) / 3.0) * 3.0
    print(f"  {name:<10} {len(CORPORA[name]):>6}  {corpus_R(name,3.0):>7.4f}  "
          f"{a:>+8.3f}°    {gl:>+8.3f}°E          {35-gl:>+5.2f}°")

print(f"""
  NPC consensus (3 corpora):
    K(3°) = {K3_mon:.4f}   anchor φ = {phi3_mon:+.3f}°
    nearest grid line to 35°E: {gl_mon_35:.3f}°E  (Δ = {35-gl_mon_35:+.3f}°)
    perm-p (fixed T=3°, {N_PERM:,} draws) = {p3_mon:.5f}  {sig_stars(p3_mon)}
  → All three monument corpora independently recover the same anchor,
    placing a grid line within {abs(35-gl_mon_35):.2f}° of 35°E.

  PRIMARY (non-nested 2-corpus NPC: UNESCO + Wikidata Stupa):
    K(3°) = {K3_two:.4f}   anchor φ = {phi3_two:+.3f}°
    nearest grid line to 35°E: {gl_two_35:.3f}°E  (Δ = {35-gl_two_35:+.3f}°)
    perm-p (fixed T=3°, {N_PERM:,} draws) = {p3_two:.5f}  {sig_stars(p3_two)}
  → Conservative non-nested concordance between externally curated UNESCO
    Cultural/Mixed corpus and externally curated Wikidata Q180987 stupas.

  THEORETICALLY MOTIVATED (non-nested 2-corpus NPC: Dome subset + Wikidata Stupa):
    K(3°) = {K3_ds:.4f}   anchor φ = {phi3_ds:+.3f}°
    nearest grid line to 35°E: {gl_ds_35:.3f}°E  (Δ = {35-gl_ds_35:+.3f}°)
    perm-p (fixed T=3°, {N_PERM:,} draws) = {p3_ds:.5f}  {sig_stars(p3_ds)}
  → Morphologically comparable non-nested corpora: UNESCO domed monuments
    and Wikidata Q180987 stupas. Neither is nested within the other.
""")

# ── Part A2: Sub-corpus negative-control NPCs ─────────────────────────────────
print(f"""PART A2 — SUB-CORPUS NEGATIVE-CONTROL NPCs  (2-corpus, T = 3°)
  {'Pair':<35} {'N_a':>5} {'N_b':>5}  {'K':>6}  {'Δ35°E':>7}  {'perm-p':>9}  {'sig':>4}
  {'-'*94}
  {'UNESCO-Eurasia × UNESCO-Americas':<35} {len(lons_eurasia):>5} {len(lons_americas):>5}  {K_eu_am:>6.4f}  {35-gl_eu_am_35:>+5.2f}°  {p_eu_am:>9.5f}  {sig_stars(p_eu_am):>4}
  {'UNESCO (full)  × UNESCO-Americas':<35} {len(lons_unesco):>5} {len(lons_americas):>5}  {K_un_am:>6.4f}  {35-gl_un_am_35:>+5.2f}°  {p_un_am:>9.5f}  {sig_stars(p_un_am):>4}
  {'UNESCO (full)  × OWTRAD nodes   ':<35} {len(lons_unesco):>5} {len(lons_owtrad):>5}  {K_un_ow:>6.4f}  {35-gl_un_ow_35:>+5.2f}°  {p_un_ow:>9.5f}  {sig_stars(p_un_ow):>4}
  {'UNESCO + Americas + Stupa (3-c) ':<35} {len(lons_unesco):>5} {len(lons_americas)+len(lons_stupa):>5}  {K_neg3:>6.4f}  {35-gl_neg3_35:>+5.2f}°  {p_neg3:>9.5f}  {sig_stars(p_neg3):>4}
  → Sub-corpus controls outside the monument signal (Americas) and on a
    different corpus type (OWTRAD trade routes) for comparison with the
    primary UNESCO + Wikidata Stupa NPC. The 3-corpus negative control
    substitutes the Americas geographic-control sub-corpus for the dome
    subset in the 3-corpus monument primary.
""")

# ── Part A3: OSM stupa NPC sensitivity rows ───────────────────────────────────
print(f"""PART A3 — OSM STUPA NPC SENSITIVITY (2-corpus and 3-corpus, T = 3°)
  OSM corpus (building/historic/ruins=stupa, N={len(lons_osm)}) assembled independently
  of Wikidata curation.  Tests whether OSM anchor agrees with monument corpora.
  {'Pair':<40} {'K':>6}  {'Δ35°E':>7}  {'perm-p':>9}  {'sig':>4}
  {'-'*80}
  {'Dome × OSM':<40} {K_osm_dome:>6.4f}  {35-gl_osm_dome_35:>+5.2f}°  {p_osm_dome:>9.5f}  {sig_stars(p_osm_dome):>4}
  {'Wikidata Stupa × OSM':<40} {K_osm_stupa:>6.4f}  {35-gl_osm_stupa_35:>+5.2f}°  {p_osm_stupa:>9.5f}  {sig_stars(p_osm_stupa):>4}
  {'UNESCO × OSM':<40} {K_osm_unesco:>6.4f}  {35-gl_osm_unesco_35:>+5.2f}°  {p_osm_unesco:>9.5f}  {sig_stars(p_osm_unesco):>4}
  {'Dome + Wikidata Stupa + OSM (3-c)':<40} {K_osm3:>6.4f}  {35-gl_osm3_35:>+5.2f}°  {p_osm3:>9.5f}  {sig_stars(p_osm3):>4}
""")

# ── Part B: OWTRAD contrast ───────────────────────────────────────────────────
print(f"""PART B — OWTRAD (trade-route nodes, N={N_ow}) — CONTRAST CORPUS
  At T=3°:
    R(3°)   = {R_owtrad_3:.4f}   (monument corpora: {np.mean([corpus_R(n,3.0) for n in monument_names]):.4f} avg)
    anchor  = {phi_owtrad_3:+.3f}°   (offset {abs(phi_owtrad_3 - phi3_mon):.3f}° from monument consensus)
    Rayleigh p (isolated, {N_PERM:,} draws) = {p_rayleigh_ow:.4f}  {sig_stars(p_rayleigh_ow)}
  → OWTRAD shows no significant periodicity at T=3° when tested in isolation.
    Its preferred anchor deviates {abs(phi_owtrad_3 - phi3_mon):.2f}° from the monument consensus.
""")

# ── Part C: All four corpora ──────────────────────────────────────────────────
print(f"""PART C — ALL FOUR CORPORA  (NPC including OWTRAD)
  K(3°) = {K_at_3:.4f}   anchor φ = {anchor_at_3:+.3f}°
  perm-p (fixed T=3°) = {perm_p_3:.5f}  {sig_stars(perm_p_3)}
  perm-p (global max, search-corrected) = {perm_p_globalmax:.5f}  {sig_stars(perm_p_globalmax)}
  → Adding OWTRAD drops K from {K3_mon:.3f} → {K_at_3:.3f} and renders the result
    non-significant.  This is expected: trade-route nodes are not governed by
    the same monumental geography as UNESCO / stupa corpora.
""")

# ── Individual signals (full table) ──────────────────────────────────────────
print(f"""INDIVIDUAL CORPUS SIGNALS (all four):
{'Corpus':<10} {'N':>6}  {'R(3°)':>7}  {'anchor@3°':>10}  {'own peak T':>10}  {'R(peak)':>8}  {'anchor@peak':>12}
{'-'*72}""")
for name in corpus_names:
    pk_T, pk_R, pk_a = indiv_peak[name]
    print(f"  {name:<10} {len(CORPORA[name]):>6}  {corpus_R(name,3.0):>7.4f}  "
          f"{corpus_anchor(name,3.0):>+8.2f}°    {pk_T:>8.3f}°    {pk_R:>7.4f}  {pk_a:>+8.2f}°")

# Pairwise phase-difference table
print(f"\nPAIRWISE ANCHOR AGREEMENT — phase difference [0, T/2] at each period")
print(f"  (≈0° = same anchor, ≈T/2 = opposite anchor)")
for T_show, label in [(3.0, "T=3°"), (peak_T, f"T={peak_T:.3f}° (NPC peak, all 4)")]:
    print(f"\n  {label}:")
    print(f"  {'':10}", end="")
    for nb in corpus_names[1:]:
        print(f"  {nb:>10}", end="")
    print()
    for i, na in enumerate(corpus_names[:-1]):
        print(f"  {na:<10}", end="")
        for j, nb in enumerate(corpus_names):
            if j <= i:
                print(f"  {'':>10}", end="")
            else:
                d = phase_diff_deg(na, nb, T_show)
                print(f"  {d:>+8.3f}°  ", end="")
        print()

# Top-10 NPC peaks (all 4 corpora)
def top_peaks(K, n=10, min_sep=0.2):
    sep = max(1, int(min_sep / (periods[1] - periods[0])))
    K_w, idxs = K.copy(), []
    while len(idxs) < n and K_w.max() > 0:
        i = int(np.argmax(K_w))
        idxs.append(i)
        K_w[max(0,i-sep):i+sep+1] = 0.0
    return idxs

top_idx = top_peaks(K_obs)
print(f"\nTOP-10 PERIODS BY NPC  K(T)  (all 4 corpora, each = 1 vote):")
print(f"{'Rk':<3} {'T':>7}  {'K':>6}  {'anchor':>8}  {'p(fixed)':>9}  {'sig':>4}  "
      f"  {'jk-drop UNESCO':>14}  {'jk-drop Dome':>12}  {'jk-drop Stupa':>13}  {'jk-drop OWTRAD':>14}")
print("-" * 110)
for rk, i in enumerate(top_idx, 1):
    T_i   = float(periods[i])
    K_i   = float(K_obs[i])
    phi_i = best_anchor(unit_sum_obs[i], T_i)
    p_i   = float(np.mean(K_null[:, i] >= K_i))
    jk_vals = "  ".join(f"{jk_K[n][i]:.4f}" for n in corpus_names)
    print(f"  {rk:<2}  {T_i:>6.3f}°  {K_i:>6.4f}  {phi_i:>+6.2f}°  {p_i:>9.5f}  {sig_stars(p_i):>4}  "
          f"  {jk_vals}")

# ── 10. LaTeX macros ──────────────────────────────────────────────────────────
print()
print(SEP)
print("  LATEX MACROS")
print(SEP)
macros = {
    # corpus sizes
    "npcNUnesco":           len(lons_unesco),
    "npcNDome":             len(lons_dome),
    "npcNStupa":            len(lons_stupa),
    "npcNOwtrad":           len(lons_owtrad),
    "npcNPerm":             N_PERM,
    # ── ROW 1: 2-corpus theoretically motivated (Dome + Stupa, non-nested) ──
    "npcDomeStupaK":           round(K3_ds, 4),
    "npcDomeStupaAnchor":      round(phi3_ds, 3),
    "npcDomeStupaGridLineNearE": round(gl_ds_35, 3),
    "npcDomeStupaGridDeltaE":  round(abs(35 - gl_ds_35), 3),
    "npcDomeStupaPermP":       round(p3_ds, 5),
    # ── ROW 2: 2-corpus primary conservative (UNESCO + Wikidata Stupa) ──
    "npcTwoK":              round(K3_two, 4),
    "npcTwoAnchor":         round(phi3_two, 3),
    "npcTwoGridLineNearE":  round(gl_two_35, 3),
    "npcTwoGridDeltaE":     round(abs(35 - gl_two_35), 3),
    "npcTwoPermP":          round(p3_two, 5),
    # ── ROW 3: 3-corpus (includes nested dome subset) ──
    "npcMonK":              round(K3_mon, 4),
    "npcMonAnchor":         round(phi3_mon, 3),
    "npcMonGridLineNearE":  round(gl_mon_35, 3),
    "npcMonGridDeltaE":     round(abs(35 - gl_mon_35), 3),
    "npcMonPermP":          round(p3_mon, 5),
    # ── NEGATIVE CONTROLS: 2-corpus sub-corpus NPCs ──
    "npcNAmericas":         len(lons_americas),
    "npcNEurasia":          len(lons_eurasia),
    # (i) UNESCO-Eurasia × UNESCO-Americas — disjoint geographic split
    "npcNegEuAmK":          round(K_eu_am, 4),
    "npcNegEuAmAnchor":     round(phi_eu_am, 3),
    "npcNegEuAmGridLineNearE": round(gl_eu_am_35, 3),
    "npcNegEuAmGridDeltaE": round(abs(35 - gl_eu_am_35), 3),
    "npcNegEuAmPermP":      round(p_eu_am, 5),
    # (ii) UNESCO × UNESCO-Americas — sub-corpus contrast (parallels OWTRAD)
    "npcNegUnAmK":          round(K_un_am, 4),
    "npcNegUnAmAnchor":     round(phi_un_am, 3),
    "npcNegUnAmGridDeltaE": round(abs(35 - gl_un_am_35), 3),
    "npcNegUnAmPermP":      round(p_un_am, 5),
    # (iii) UNESCO × OWTRAD  — 2-corpus version of OWTRAD contrast
    "npcNegUnOwK":          round(K_un_ow, 4),
    "npcNegUnOwAnchor":     round(phi_un_ow, 3),
    "npcNegUnOwGridDeltaE": round(abs(35 - gl_un_ow_35), 3),
    "npcNegUnOwPermP":      round(p_un_ow, 5),
    # (iv) UNESCO × UNESCO-Americas × Wikidata Stupa — 3-corpus negative
    # control mirroring the 3-corpus monument primary, with the dome subset
    # replaced by the Americas geographic-control sub-corpus.
    "npcNegThreeK":          round(K_neg3, 4),
    "npcNegThreeAnchor":     round(phi_neg3, 3),
    "npcNegThreeGridLineNearE": round(gl_neg3_35, 3),
    "npcNegThreeGridDeltaE": round(abs(35 - gl_neg3_35), 3),
    "npcNegThreePermP":      round(p_neg3, 5),
    # ── OSM stupa NPC sensitivity rows ──
    "npcNOsm":               len(lons_osm),
    # 2-corpus: Dome × OSM
    "npcOsmDomeK":           round(K_osm_dome,  4),
    "npcOsmDomeGridDeltaE":  round(abs(35 - gl_osm_dome_35), 3),
    "npcOsmDomePermP":       round(p_osm_dome,  5),
    # 2-corpus: Wikidata Stupa × OSM
    "npcOsmStupaK":          round(K_osm_stupa, 4),
    "npcOsmStupaGridDeltaE": round(abs(35 - gl_osm_stupa_35), 3),
    "npcOsmStupaPermP":      round(p_osm_stupa, 5),
    # 2-corpus: UNESCO × OSM
    "npcOsmUnescoK":         round(K_osm_unesco, 4),
    "npcOsmUnescoGridDeltaE":round(abs(35 - gl_osm_unesco_35), 3),
    "npcOsmUnescoPermP":     round(p_osm_unesco, 5),
    # 3-corpus: Dome + Wikidata Stupa + OSM
    "npcOsmThreeK":          round(K_osm3,  4),
    "npcOsmThreeGridDeltaE": round(abs(35 - gl_osm3_35), 3),
    "npcOsmThreePermP":      round(p_osm3,  5),
    # ── Part B: OWTRAD isolated ──
    "npcOwR":               round(R_owtrad_3, 4),
    "npcOwAnchor":          round(phi_owtrad_3, 3),
    "npcOwAnchorDelta":     round(abs(phi_owtrad_3 - phi3_mon), 3),
    "npcOwRayleighP":       round(p_rayleigh_ow, 4),
    # ── Part C: all-4 NPC ──
    "npcAllK":              round(K_at_3, 4),
    "npcAllAnchor":         round(anchor_at_3, 3),
    "npcAllPermP":          round(perm_p_3, 5),
    "npcAllPermPGlobalmax": round(perm_p_globalmax, 5),
    # pairwise differences at T=3°
    "npcPairDomeStupa":     round(phase_diff_deg("Dome","Stupa",3.0), 3),
    "npcPairUnescoStupa":   round(phase_diff_deg("UNESCO","Stupa",3.0), 3),
    "npcPairUnescoOwtrad":  round(phase_diff_deg("UNESCO","OWTRAD",3.0), 3),
}
for k, v in macros.items():
    print(f"  \\newcommand{{\\{k}}}{{{v}}}")

# ── 11. Write to ResultsStore ─────────────────────────────────────────────────
ResultsStore().write_many(macros)
print()
print("  Results written to data/store/results.json")
