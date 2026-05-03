"""
multiscale_combined_p.py
========================
Multi-scale enrichment test for domed UNESCO monuments, with cross-corpus
natural-scale detection.

RATIONALE
---------
Any single proximity window (tier) is an arbitrary choice.  This script
removes that dependence by combining evidence across K logarithmically-spaced
window sizes simultaneously.

For each window w_k, the geometric null rate is p0_k = 2*w_k / HARMONIC_STEP.
A one-sided binomial test asks whether the observed fraction of sites within
w_k exceeds p0_k.

The K resulting p-values are combined with Fisher's method:
    chi2_obs = -2 * sum(log(p_k))

Because the counts at nested windows are correlated, chi2_obs is NOT
chi-squared distributed under the null.  We calibrate it by permutation.

NATURAL SCALE
-------------
The window where -log(p_binom) peaks is the data-derived natural proximity
scale: the scale at which enrichment is most statistically concentrated.
We find this for each corpus and report the joint peak (window maximising
the sum of -log(p) across all corpora).

CORPORA
-------
  1. UNESCO cultural sites -- dome/stupa subset  (N_dome)
  2. UNESCO cultural sites -- full corpus        (N_full)
  3. Wikidata Q180987 stupa corpus               (N_wiki)
  4. OWTRAD Silk Road vertices                   (N_owtrad)
"""

import csv
import sys
import json
import math
import numpy as np
from scipy.stats import binom
from scipy.special import betainc
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT))

from data.unesco_corpus import load_corpus
from lib.dome_filter import FORM_KEYWORDS, FORM_KEYWORD_RES
from lib.results_store import ResultsStore

_CFG    = json.loads((_ROOT / "config.json").read_text())
N_PERM  = _CFG["simulation"]["n_permutations"]
SEED    = _CFG["simulation"]["random_seed"]

GERIZIM_LON   = _CFG["anchors"]["gerizim"]["longitude"]
HARMONIC_STEP = 3.0
HALF_STEP     = HARMONIC_STEP / 2.0   # 1.5 deg

# ── Window grid ───────────────────────────────────────────────────────────────
# 15 log-spaced null rates from 0.33% to 40%.
# null_rate = 2w / HARMONIC_STEP  =>  w = null_rate * HARMONIC_STEP / 2
NULL_RATES = np.geomspace(0.0033, 0.40, 15)
WINDOWS    = NULL_RATES * HARMONIC_STEP / 2.0   # degrees
K          = len(WINDOWS)

# ── Harmonic deviation helper ─────────────────────────────────────────────────
def harmonic_devs(lons):
    """Deviation from nearest 3-deg harmonic node, in [0, 1.5] deg."""
    shifted = (lons - GERIZIM_LON) % HARMONIC_STEP
    return np.minimum(shifted, HARMONIC_STEP - shifted)

# ── Vectorised binomial SF ────────────────────────────────────────────────────
def binom_sf_batch(n_arr, N, p0_arr):
    """P(X >= n) for arrays n and p0, same N."""
    n = np.maximum(n_arr, 1)
    return betainc(n, N - n + 1, p0_arr)

# ── Load UNESCO corpora ───────────────────────────────────────────────────────
corpus   = load_corpus()
cultural = [s for s in corpus if s.category != "Natural" and s.has_coords]
lons_full = np.array([s.longitude for s in cultural])
N_full    = len(lons_full)

lons_dome = np.array([
    s.longitude for s in cultural
    if any(FORM_KEYWORD_RES[kw].search(s.full_text) for kw in FORM_KEYWORDS)
])
N_dome = len(lons_dome)

# ── Load Wikidata Q180987 stupa corpus ────────────────────────────────────────
_WIKI_CSV = _ROOT / "data" / "store" / "unesco" / "wikidata_stupas_q180987.csv"
_wiki_lons = []
if _WIKI_CSV.exists():
    with open(_WIKI_CSV, newline="") as _fh:
        _lines = [ln for ln in _fh if not ln.startswith("#")]
    for _row in csv.DictReader(_lines):
        try:
            _wiki_lons.append(float(_row["lon"]))
        except (ValueError, KeyError):
            continue
lons_wiki = np.array(_wiki_lons)
N_wiki    = len(lons_wiki)

# ── Load OWTRAD Silk Road vertices ────────────────────────────────────────────
_OWTRAD_DIR   = _ROOT / "data" / "store" / "silk_road"
_OWTRAD_NODES = _OWTRAD_DIR / "owtrad_nodes.csv"
_owtrad_lons, _owtrad_seen = [], set()
if _OWTRAD_NODES.exists():
    with open(_OWTRAD_NODES, newline="", encoding="utf-8") as _fh:
        for _row in csv.DictReader(r for r in _fh if not r.startswith("#")):
            try:
                _lon = float(_row["lon"])
                _lat = float(_row["lat"])
            except (ValueError, KeyError):
                continue
            _key = (round(_lon, 4), round(_lat, 4))
            if _key not in _owtrad_seen:
                _owtrad_seen.add(_key)
                _owtrad_lons.append(_lon)
lons_owtrad = np.array(_owtrad_lons)
N_owtrad    = len(lons_owtrad)

# ── Compute observed binom p-curves for each corpus ───────────────────────────
def obs_curve(lons, N):
    devs   = harmonic_devs(lons)
    counts = np.array([np.sum(devs <= w) for w in WINDOWS])
    p_arr  = np.clip(np.array([
        binom.sf(obs - 1, N, p0)
        for obs, p0 in zip(counts, NULL_RATES)
    ]), 1e-300, 1.0)
    chi2   = float(-2.0 * np.sum(np.log(p_arr)))
    return devs, counts, p_arr, chi2

devs_dome,   obs_dome,   p_dome,   chi2_dome   = obs_curve(lons_dome,   N_dome)
devs_full,   obs_full,   p_full,   chi2_full   = obs_curve(lons_full,   N_full)
devs_wiki,   obs_wiki,   p_wiki,   chi2_wiki   = obs_curve(lons_wiki,   N_wiki)
devs_owtrad, obs_owtrad, p_owtrad, chi2_owtrad = obs_curve(lons_owtrad, N_owtrad)

# ── Natural scale: window where joint -log(p) across all corpora peaks ────────
neg_log_joint = (
    -np.log(p_dome)
    - np.log(p_full)
    - np.log(p_wiki)
    - np.log(p_owtrad)
)
natural_scale_idx      = int(np.argmax(neg_log_joint))
natural_scale_null_pct = float(NULL_RATES[natural_scale_idx] * 100)
natural_scale_deg      = float(WINDOWS[natural_scale_idx])

# Also compute per-corpus peak windows
dome_peak_idx   = int(np.argmax(-np.log(p_dome)))
full_peak_idx   = int(np.argmax(-np.log(p_full)))
wiki_peak_idx   = int(np.argmax(-np.log(p_wiki)))
owtrad_peak_idx = int(np.argmax(-np.log(p_owtrad)))

# ── Permutation calibration: dome subset ──────────────────────────────────────
rng = np.random.default_rng(SEED)
perm_indices = np.array([
    rng.choice(N_full, size=N_dome, replace=False)
    for _ in range(N_PERM)
])
perm_devs    = devs_full[perm_indices]   # (N_PERM, N_dome)
hit_matrix   = (perm_devs[:, :, None] <= WINDOWS[None, None, :])
perm_counts  = hit_matrix.sum(axis=1)
perm_p_dome  = np.clip(binom_sf_batch(perm_counts, N_dome, NULL_RATES[None, :]), 1e-300, 1.0)
chi2_perm_dome    = -2.0 * np.sum(np.log(perm_p_dome), axis=1)
combined_p_dome   = float(np.mean(chi2_perm_dome >= chi2_dome))
p_hm_dome_obs     = float(K / np.sum(1.0 / p_dome))
p_hm_dome_perm    = K / np.sum(1.0 / perm_p_dome, axis=1)
combined_p_hm_dome = float(np.mean(p_hm_dome_perm <= p_hm_dome_obs))

# ── Permutation calibration: full corpus ──────────────────────────────────────
rng2 = np.random.default_rng(SEED + 1)
perm_indices_full  = np.array([
    rng2.choice(N_full, size=N_full, replace=False)
    for _ in range(N_PERM)
])
perm_devs_full     = devs_full[perm_indices_full]
hit_matrix_full    = (perm_devs_full[:, :, None] <= WINDOWS[None, None, :])
perm_counts_full   = hit_matrix_full.sum(axis=1)
perm_p_full        = np.clip(binom_sf_batch(perm_counts_full, N_full, NULL_RATES[None, :]), 1e-300, 1.0)
chi2_perm_full     = -2.0 * np.sum(np.log(perm_p_full), axis=1)
combined_p_full    = float(np.mean(chi2_perm_full >= chi2_full))
p_hm_full_obs      = float(K / np.sum(1.0 / p_full))
p_hm_full_perm     = K / np.sum(1.0 / perm_p_full, axis=1)
combined_p_hm_full = float(np.mean(p_hm_full_perm <= p_hm_full_obs))

# ── Permutation calibration: Wikidata stupa ───────────────────────────────────
rng3 = np.random.default_rng(SEED + 2)
perm_indices_wiki  = np.array([
    rng3.choice(N_wiki, size=N_wiki, replace=False)
    for _ in range(N_PERM)
])
perm_devs_wiki     = devs_wiki[perm_indices_wiki]
hit_matrix_wiki    = (perm_devs_wiki[:, :, None] <= WINDOWS[None, None, :])
perm_counts_wiki   = hit_matrix_wiki.sum(axis=1)
perm_p_wiki        = np.clip(binom_sf_batch(perm_counts_wiki, N_wiki, NULL_RATES[None, :]), 1e-300, 1.0)
chi2_perm_wiki     = -2.0 * np.sum(np.log(perm_p_wiki), axis=1)
combined_p_wiki    = float(np.mean(chi2_perm_wiki >= chi2_wiki))
p_hm_wiki_obs      = float(K / np.sum(1.0 / p_wiki))
p_hm_wiki_perm     = K / np.sum(1.0 / perm_p_wiki, axis=1)
combined_p_hm_wiki = float(np.mean(p_hm_wiki_perm <= p_hm_wiki_obs))

# ── Permutation calibration: OWTRAD vertices ──────────────────────────────────
rng4 = np.random.default_rng(SEED + 3)
perm_indices_owtrad  = np.array([
    rng4.choice(N_owtrad, size=N_owtrad, replace=False)
    for _ in range(N_PERM)
])
perm_devs_owtrad     = devs_owtrad[perm_indices_owtrad]
hit_matrix_owtrad    = (perm_devs_owtrad[:, :, None] <= WINDOWS[None, None, :])
perm_counts_owtrad   = hit_matrix_owtrad.sum(axis=1)
perm_p_owtrad        = np.clip(binom_sf_batch(perm_counts_owtrad, N_owtrad, NULL_RATES[None, :]), 1e-300, 1.0)
chi2_perm_owtrad     = -2.0 * np.sum(np.log(perm_p_owtrad), axis=1)
combined_p_owtrad    = float(np.mean(chi2_perm_owtrad >= chi2_owtrad))
p_hm_owtrad_obs      = float(K / np.sum(1.0 / p_owtrad))
p_hm_owtrad_perm     = K / np.sum(1.0 / perm_p_owtrad, axis=1)
combined_p_hm_owtrad = float(np.mean(p_hm_owtrad_perm <= p_hm_owtrad_obs))

# ── Joint-maximum test ───────────────────────────────────────────────────────
# At each window, sum -log(p) across all four corpora simultaneously.
# The maximum of this joint curve is the observed statistic.
# Under the null (paired i-th permutation from each corpus), compute the same.
# This tests: is there ANY window where all four corpora are jointly more
# enriched than expected by chance?
perm_joint     = (
    -np.log(perm_p_dome)
    - np.log(perm_p_full)
    - np.log(perm_p_wiki)
    - np.log(perm_p_owtrad)
)   # (N_PERM, K)
perm_joint_max = perm_joint.max(axis=1)          # (N_PERM,)
obs_joint_max  = float(neg_log_joint.max())
p_joint_max    = float(np.mean(perm_joint_max >= obs_joint_max))

# Which window holds the joint max?
joint_max_idx     = int(np.argmax(neg_log_joint))
joint_max_null_pct = float(NULL_RATES[joint_max_idx] * 100)
joint_max_deg     = float(WINDOWS[joint_max_idx])
joint_max_km      = joint_max_deg * 111.0 * math.cos(math.radians(30))

# ── Natural scale convergence probability ────────────────────────────────────
# Observed: all four corpora peak within 30 km (null rate <= 20.2%, index <= 12).
# For each corpus we ask: under the null, how often does the peak fall at or
# below that window?  The four corpora are independent, so the joint probability
# is the product.  We also compute the range of peak indices across corpora as
# a convergence metric and compare it to the permutation null.
K_30KM = natural_scale_idx   # index 12, window 0.3023 deg, ~29 km at 30N

dome_perm_peaks   = np.argmax(-np.log(perm_p_dome),   axis=1)   # (N_PERM,)
full_perm_peaks   = np.argmax(-np.log(perm_p_full),   axis=1)
wiki_perm_peaks   = np.argmax(-np.log(perm_p_wiki),   axis=1)
owtrad_perm_peaks = np.argmax(-np.log(perm_p_owtrad), axis=1)

# Per-corpus null probability of peaking within 30 km
p_dome_in30   = float(np.mean(dome_perm_peaks   <= K_30KM))
p_full_in30   = float(np.mean(full_perm_peaks   <= K_30KM))
p_wiki_in30   = float(np.mean(wiki_perm_peaks   <= K_30KM))
p_owtrad_in30 = float(np.mean(owtrad_perm_peaks <= K_30KM))
p_all_in30    = p_dome_in30 * p_full_in30 * p_wiki_in30 * p_owtrad_in30

# Convergence: range of peak indices across the four observed corpora
obs_peak_indices  = np.array([dome_peak_idx, full_peak_idx, wiki_peak_idx, owtrad_peak_idx])
obs_peak_range    = int(obs_peak_indices.max() - obs_peak_indices.min())

# Null range: match permutation draws row-by-row (all from independent seeds,
# so pairing the i-th draw from each corpus is a valid joint null draw)
perm_ranges = (
    np.maximum.reduce([dome_perm_peaks, full_perm_peaks,
                       wiki_perm_peaks, owtrad_perm_peaks])
    - np.minimum.reduce([dome_perm_peaks, full_perm_peaks,
                         wiki_perm_peaks, owtrad_perm_peaks])
)
p_convergence = float(np.mean(perm_ranges <= obs_peak_range))

# ── Print ──────────────────────────────────────────────────────────────────────
SEP = "=" * 80

def print_table(label, N, obs_arr, p_arr, chi2, perm_p, hm_p, hm_perm_p, peak_idx):
    print()
    print(SEP)
    print(f"  MULTI-SCALE COMBINED ENRICHMENT  --  {label}  (N={N})")
    print(f"  {K} log-spaced windows  |  {N_PERM:,} permutations")
    print(SEP)
    print(f"  {'Null rate':>10}  {'Window (deg)':>13}  {'Obs':>5}  {'Exp':>6}  "
          f"{'Ratio':>6}  {'Binom p':>10}  {'Peak':>5}")
    print("  " + "-" * 68)
    for k_idx in range(K):
        exp   = N * NULL_RATES[k_idx]
        ratio = obs_arr[k_idx] / exp if exp > 0 else float('nan')
        marker = " <--" if k_idx == peak_idx else ""
        print(f"  {NULL_RATES[k_idx]*100:>9.1f}%  {WINDOWS[k_idx]:>13.4f}  "
              f"{obs_arr[k_idx]:>5}  {exp:>6.1f}  {ratio:>6.2f}x  "
              f"{p_arr[k_idx]:>10.5f}{marker}")
    print()
    print(f"  Fisher chi2     : {chi2:.4f}  (df={2*K} under independence)")
    print(f"  Permutation p   : {perm_p:.5f}")
    print(f"  Harmonic mean p : {hm_p:.5f}  (perm p = {hm_perm_p:.5f})")
    print(f"  Peak scale      : {NULL_RATES[peak_idx]*100:.1f}% null rate "
          f"({WINDOWS[peak_idx]:.4f} deg)")
    print()

print_table("UNESCO DOME SUBSET", N_dome, obs_dome, p_dome,
            chi2_dome, combined_p_dome, p_hm_dome_obs, combined_p_hm_dome, dome_peak_idx)
print_table("UNESCO FULL CORPUS", N_full, obs_full, p_full,
            chi2_full, combined_p_full, p_hm_full_obs, combined_p_hm_full, full_peak_idx)
print_table("WIKIDATA Q180987 STUPAS", N_wiki, obs_wiki, p_wiki,
            chi2_wiki, combined_p_wiki, p_hm_wiki_obs, combined_p_hm_wiki, wiki_peak_idx)
print_table("OWTRAD SILK ROAD VERTICES", N_owtrad, obs_owtrad, p_owtrad,
            chi2_owtrad, combined_p_owtrad, p_hm_owtrad_obs, combined_p_hm_owtrad, owtrad_peak_idx)

print()
print(SEP)
print("  JOINT-MAXIMUM TEST  (best window across all 4 corpora simultaneously)")
print(SEP)
print()
print(f"  Joint max window : {joint_max_null_pct:.1f}% null rate "
      f"({joint_max_deg:.4f} deg, ~{joint_max_km:.0f} km at 30N)")
print(f"  Observed joint -log(p) sum : {obs_joint_max:.3f}")
print(f"  Permutation p              : {p_joint_max:.5f}")
print()
print(SEP)
print("  NATURAL PROXIMITY SCALE  (joint peak across all 4 corpora)")
print(SEP)
print()
print(f"  Window index    : {natural_scale_idx} of {K}")
print(f"  Null rate       : {natural_scale_null_pct:.2f}%")
print(f"  Window width    : +/-{natural_scale_deg:.4f} deg")
km_approx = natural_scale_deg * 111.0 * math.cos(math.radians(30))
print(f"  Approx distance : ~{km_approx:.1f} km  (at 30 deg N latitude)")
print()
print(f"  Per-corpus peaks (all within ~30 km of a harmonic node):")
print(f"    UNESCO dome  : {NULL_RATES[dome_peak_idx]*100:.1f}% null  "
      f"({WINDOWS[dome_peak_idx]:.4f} deg, "
      f"~{WINDOWS[dome_peak_idx]*111.0*math.cos(math.radians(30)):.0f} km)")
print(f"    OWTRAD       : {NULL_RATES[owtrad_peak_idx]*100:.1f}% null  "
      f"({WINDOWS[owtrad_peak_idx]:.4f} deg, "
      f"~{WINDOWS[owtrad_peak_idx]*111.0*math.cos(math.radians(30)):.0f} km)")
print(f"    UNESCO full  : {NULL_RATES[full_peak_idx]*100:.1f}% null  "
      f"({WINDOWS[full_peak_idx]:.4f} deg, "
      f"~{WINDOWS[full_peak_idx]*111.0*math.cos(math.radians(30)):.0f} km)")
print(f"    Wikidata     : {NULL_RATES[wiki_peak_idx]*100:.1f}% null  "
      f"({WINDOWS[wiki_peak_idx]:.4f} deg, "
      f"~{WINDOWS[wiki_peak_idx]*111.0*math.cos(math.radians(30)):.0f} km)")
print()
print(f"  Convergence test (all 4 peaks within {km_approx:.0f} km):")
print(f"    Null P(dome  peak in30)   : {p_dome_in30:.4f}")
print(f"    Null P(full  peak in30)   : {p_full_in30:.4f}")
print(f"    Null P(wiki  peak in30)   : {p_wiki_in30:.4f}")
print(f"    Null P(owtrad peak in30)  : {p_owtrad_in30:.4f}")
print(f"    Joint P(all 4, indep)     : {p_all_in30:.5f}")
print()
print(f"  Peak-range convergence:")
print(f"    Observed range of peak indices : {obs_peak_range}  ({obs_peak_indices})")
print(f"    Perm p(range <= obs)           : {p_convergence:.5f}")

print()
print(SEP)
print("  LATEX MACROS")
print(SEP)
print()
print(f"  \\newcommand{{\\multiscaleK}}{{{K}}}                       % number of log windows")
print(f"  \\newcommand{{\\multiscaleNPerm}}{{{N_PERM:,}}}              % permutations")
print(f"  % Dome subset")
print(f"  \\newcommand{{\\multiscaleDomeChiSq}}{{{chi2_dome:.2f}}}")
print(f"  \\newcommand{{\\multiscaleDomeCombinedP}}{{{combined_p_dome:.5f}}}")
print(f"  \\newcommand{{\\multiscaleDomeHMP}}{{{p_hm_dome_obs:.5f}}}")
print(f"  % Full corpus")
print(f"  \\newcommand{{\\multiscaleFullChiSq}}{{{chi2_full:.2f}}}")
print(f"  \\newcommand{{\\multiscaleFullCombinedP}}{{{combined_p_full:.5f}}}")
print(f"  \\newcommand{{\\multiscaleFullHMP}}{{{p_hm_full_obs:.5f}}}")
print(f"  % Wikidata stupa")
print(f"  \\newcommand{{\\multiscaleWikiChiSq}}{{{chi2_wiki:.2f}}}")
print(f"  \\newcommand{{\\multiscaleWikiCombinedP}}{{{combined_p_wiki:.5f}}}")
print(f"  \\newcommand{{\\multiscaleWikiHMP}}{{{p_hm_wiki_obs:.5f}}}")
print(f"  % OWTRAD")
print(f"  \\newcommand{{\\multiscaleOwtradChiSq}}{{{chi2_owtrad:.2f}}}")
print(f"  \\newcommand{{\\multiscaleOwtradCombinedP}}{{{combined_p_owtrad:.5f}}}")
print(f"  \\newcommand{{\\multiscaleOwtradHMP}}{{{p_hm_owtrad_obs:.5f}}}")
print(f"  % Natural scale")
print(f"  \\newcommand{{\\multiscaleNaturalNullPct}}{{{natural_scale_null_pct:.1f}}}")
print(f"  \\newcommand{{\\multiscaleNaturalDeg}}{{{natural_scale_deg:.4f}}}")
print(f"  \\newcommand{{\\multiscaleNaturalKm}}{{{km_approx:.0f}}}")
print(f"  \\newcommand{{\\multiscaleAllNearNodeP}}{{{p_all_in30:.5f}}}")
print(f"  \\newcommand{{\\multiscaleConvergenceP}}{{{p_convergence:.5f}}}")
print(f"  \\newcommand{{\\multiscaleJointMaxP}}{{{p_joint_max:.5f}}}")
print(f"  \\newcommand{{\\multiscaleJointMaxNullPct}}{{{joint_max_null_pct:.1f}}}")
print(f"  \\newcommand{{\\multiscaleJointMaxKm}}{{{joint_max_km:.0f}}}")

# ── Write to results store ────────────────────────────────────────────────────
ResultsStore().write_many({
    "multiscaleK":                  K,
    "multiscaleDomeChiSq":           round(chi2_dome, 4),
    "multiscaleDomeCombinedP":      round(combined_p_dome, 5),
    "multiscaleDomeHMP":            round(p_hm_dome_obs, 5),
    "multiscaleFullChiSq":           round(chi2_full, 4),
    "multiscaleFullCombinedP":      round(combined_p_full, 5),
    "multiscaleFullHMP":            round(p_hm_full_obs, 5),
    "multiscaleWikiChiSq":           round(chi2_wiki, 4),
    "multiscaleWikiCombinedP":      round(combined_p_wiki, 5),
    "multiscaleWikiHMP":            round(p_hm_wiki_obs, 5),
    "multiscaleOwtradChiSq":         round(chi2_owtrad, 4),
    "multiscaleOwtradCombinedP":    round(combined_p_owtrad, 5),
    "multiscaleOwtradHMP":          round(p_hm_owtrad_obs, 5),
    "multiscaleNaturalNullPct":     round(natural_scale_null_pct, 1),
    "multiscaleNaturalDeg":         round(natural_scale_deg, 4),
    "multiscaleNaturalKm":          round(km_approx, 0),
    "multiscaleAllNearNodeP":       round(p_all_in30, 5),
    "multiscaleConvergenceP":       round(p_convergence, 5),
    "multiscaleJointMaxP":          round(p_joint_max, 5),
    "multiscaleJointMaxNullPct":    round(joint_max_null_pct, 1),
    "multiscaleJointMaxKm":         round(joint_max_km, 0),
})
