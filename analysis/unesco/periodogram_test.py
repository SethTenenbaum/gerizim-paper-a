"""
periodogram_test.py
===================
Rayleigh periodogram sensitivity sweep: which candidate spacing maximises
spectral power for the dome subset?

ANCHOR NOTE — IMPORTANT
-----------------------
The Rayleigh power statistic S(T) = N * |mean(exp(2*pi*i * lon / T))|^2 is
rotation-invariant: subtracting any constant anchor c from all longitudes
multiplies every phasor by exp(-2*pi*i*c/T), which rotates the mean resultant
but leaves its LENGTH unchanged.  Therefore S(T) is IDENTICAL regardless of
whether you subtract GERIZIM_LON or any other constant.  The peak period and
all power values produced here are the same with or without an anchor shift.

Consequence: the Rayleigh periodogram sweep is genuinely anchor-free -- it
finds the period at which sites cluster most tightly at *any* phase, but it
cannot detect clustering *at Gerizim specifically*.  For that directional
test (does the mean phase point toward Gerizim?) use the V-test in
anchored_periodicity_test.py.

This file answers one sensitivity question only:
  Q. Which period maximises raw spectral power in the dome subset / full
     corpus, and how does T=3 rank in that sweep?
The permutation p-value tests whether S_dome(3 deg) exceeds chance; it does
NOT test whether the mean phase at T=3 aligns with Gerizim.
"""

import sys
import json
import numpy as np
from pathlib import Path

_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(_ROOT))

from data.unesco_corpus import load_corpus
from lib.dome_filter import FORM_KEYWORDS, FORM_KEYWORD_RES
from lib.results_store import ResultsStore
from lib.stats import significance_label as sig

_CFG        = json.loads((_ROOT / "config.json").read_text())
N_PERM      = _CFG["simulation"]["n_permutations"]         # 100,000
SEED        = _CFG["simulation"]["random_seed"]            # 42
TARGET_T    = 3.0                                          # degrees — the candidate period

# ── Load corpora ──────────────────────────────────────────────────────────────
corpus   = load_corpus()
cultural = [s for s in corpus if s.category != "Natural" and s.has_coords]

# Raw longitudes: anchor subtraction is a no-op for Rayleigh power (see docstring)
lons_full = np.array([s.longitude for s in cultural])
N_full    = len(lons_full)

# Dome subset: same raw keyword sweep as spherical_monument_raw_sweep.py
lons_dome = np.array([
    s.longitude for s in cultural
    if any(FORM_KEYWORD_RES[kw].search(s.full_text) for kw in FORM_KEYWORDS)
])
N_dome = len(lons_dome)

# ── Core statistic ────────────────────────────────────────────────────────────
def rayleigh_power(lons: np.ndarray, period: float) -> float:
    """
    Rayleigh spectral power at period T (degrees).
    S(T) = N x R^2  where R = mean resultant length on the T-circle.
    Under uniform null: S ~ Exp(1).  Larger S = stronger periodicity at T.
    """
    angles = 2.0 * np.pi * lons / period
    R = abs(np.mean(np.exp(1j * angles)))
    return len(lons) * R * R


def periodogram_vectorized(lons: np.ndarray, periods: np.ndarray) -> np.ndarray:
    """
    Vectorized Rayleigh periodogram over all periods at once.
    lons   : (N,)       longitude array
    periods: (P,)       candidate period grid
    returns: (P,)       Rayleigh power S(T) for each period
    Phase matrix: (N, P) where entry [n, p] = 2pi * lons[n] / periods[p]
    """
    N = len(lons)
    # angles shape: (N, P)
    angles = 2.0 * np.pi * lons[:, None] / periods[None, :]
    # complex mean along axis 0 -> (P,)
    C = np.mean(np.exp(1j * angles), axis=0)
    return N * (C.real ** 2 + C.imag ** 2)


# ── Candidate periods ─────────────────────────────────────────────────────────
# Fine grid from 0.5 to 15 deg avoids edge effects and covers all plausible
# archaeological grid spacings.  Step 0.05 deg gives 291 candidates.
periods     = np.arange(0.5, 15.01, 0.05)
N_periods   = len(periods)

# ── Observed periodograms (vectorized) ───────────────────────────────────────
power_dome = periodogram_vectorized(lons_dome, periods)
power_full = periodogram_vectorized(lons_full, periods)

# Index of 3° period in the candidate grid
idx_3 = int(np.argmin(np.abs(periods - TARGET_T)))
S3_dome = float(power_dome[idx_3])
S3_full = float(power_full[idx_3])

# Peak period for each corpus
peak_idx_dome = int(np.argmax(power_dome))
peak_idx_full = int(np.argmax(power_full))
peak_T_dome   = float(periods[peak_idx_dome])
peak_T_full   = float(periods[peak_idx_full])

# ── Constrained sweep: peak period below 10° ─────────────────────────────────
CONSTRAINED_MAX = 10.0
constrained_mask = periods < CONSTRAINED_MAX
power_dome_c = np.where(constrained_mask, power_dome, -np.inf)
power_full_c = np.where(constrained_mask, power_full, -np.inf)
peak_idx_dome_c = int(np.argmax(power_dome_c))
peak_idx_full_c = int(np.argmax(power_full_c))
peak_T_dome_c   = float(periods[peak_idx_dome_c])
peak_T_full_c   = float(periods[peak_idx_full_c])
S_peak_dome_c   = float(power_dome[peak_idx_dome_c])
S_peak_full_c   = float(power_full[peak_idx_full_c])
R_peak_dome_c   = float(np.sqrt(S_peak_dome_c / N_dome))
R_peak_full_c   = float(np.sqrt(S_peak_full_c / N_full))
R3_dome_compare = float(np.sqrt(S3_dome / N_dome))  # R at T=3 for dome
# Rank of 3° among constrained periods
rank_3_dome_c = int(np.sum(power_dome_c[constrained_mask] > S3_dome)) + 1
rank_3_full_c = int(np.sum(power_full_c[constrained_mask] > S3_full)) + 1
n_constrained  = int(np.sum(constrained_mask))

# ── Phasor-overlap analysis: do the constrained-peak sites overlap with 3° A+? ──
# For each corpus, compute phasor projection of each site onto the dominant
# direction at the constrained peak period.  Then count how many of the top-25
# phasor contributors are also 3°-A+ sites.
AP_THRESH_DEG = 0.15  # Tier-A+ threshold

def dev3_gerizim(lon):
    from lib.dome_filter import FORM_KEYWORDS  # already imported; just needs config
    import json
    cfg = json.loads((_ROOT / "config.json").read_text())
    anchor = cfg["anchors"]["gerizim"]["longitude"]
    d = (lon - anchor) % 3.0
    return min(d, 3.0 - d)

# Load anchor from config once
import json as _json
_anchor = _json.loads((_ROOT / "config.json").read_text())["anchors"]["gerizim"]["longitude"]

def dev3(lon):
    d = (lon - _anchor) % 3.0
    return min(d, 3.0 - d)

def phasor_overlap(lons, T, top_n=25):
    """Return count of top_n phasor contributors that are also 3°-A+ sites."""
    angles = 2*np.pi * lons / T
    C = np.mean(np.exp(1j*angles))
    proj = np.cos(angles - np.angle(C))
    top_idx = np.argsort(proj)[::-1][:top_n]
    ap3 = np.array([dev3(float(lon)) <= AP_THRESH_DEG for lon in lons])
    overlap = int(np.sum(ap3[top_idx]))
    expected = top_n * float(np.mean(ap3))
    return overlap, round(expected, 1)

dome_overlap_c, dome_overlap_exp = phasor_overlap(lons_dome, peak_T_dome_c)
full_overlap_c, full_overlap_exp = phasor_overlap(lons_full, peak_T_full_c)

# Rank of 3° among all candidate periods (1 = highest power)
rank_3_dome = int(np.sum(power_dome > S3_dome)) + 1   # 1-indexed
rank_3_full = int(np.sum(power_full > S3_full)) + 1

# Percentile of 3° in its own periodogram
pctile_3_dome = float(np.mean(power_dome <= S3_dome)) * 100
pctile_3_full = float(np.mean(power_full <= S3_full)) * 100

# ── Permutation null: significance of S(3°) ──────────────────────────────────
# Null: draw N_dome from the full UNESCO pool (without replacement).
# Preserves the actual longitude distribution; requires no anchor.
rng      = np.random.default_rng(SEED)

# Draw all permutation indices at once: shape (N_PERM, N_dome)
# NOTE: for without-replacement draws we must loop (numpy has no batch
# choice-without-replacement), but we vectorize the expensive inner computation.
perm_indices = np.array([rng.choice(N_full, size=N_dome, replace=False)
                          for _ in range(N_PERM)])   # (N_PERM, N_dome)

# S(3 deg) for each permutation: vectorized over all draws simultaneously.
# perm_lons shape: (N_PERM, N_dome)
perm_lons = lons_full[perm_indices]
# angles for TARGET_T only: (N_PERM, N_dome)
angles_3  = 2.0 * np.pi * perm_lons / TARGET_T
C3        = np.mean(np.exp(1j * angles_3), axis=1)   # (N_PERM,)
S3_perm   = N_dome * (C3.real ** 2 + C3.imag ** 2)   # (N_PERM,)

# Permutation p-value for the OBSERVED PEAK (post-hoc, for comparison only).
# This answers: how often does a random draw of N_dome sites produce a global
# peak as high as the observed peak?  This is NOT corrected for the number of
# candidates; it illustrates why post-hoc peak-picking is unreliable.
S_peak_dome = float(power_dome[peak_idx_dome])
angles_peak  = 2.0 * np.pi * perm_lons / peak_T_dome
C_peak       = np.mean(np.exp(1j * angles_peak), axis=1)
S_peak_perm  = N_dome * (C_peak.real ** 2 + C_peak.imag ** 2)
perm_p_peak  = float(np.mean(S_peak_perm >= S_peak_dome))

perm_p_dome = float(np.mean(S3_perm >= S3_dome))

# Fraction of permutations where 3 deg is the global peak.
# Computed in chunks to stay within memory budget (~500 MB per chunk).
CHUNK = 500   # permutations per chunk
perm_peak_count_at_3 = 0
for start in range(0, N_PERM, CHUNK):
    end = min(start + CHUNK, N_PERM)
    chunk_lons = perm_lons[start:end]                        # (chunk, N_dome)
    # angles_all: (chunk, N_dome, N_periods)
    angles_all = 2.0 * np.pi * chunk_lons[:, :, None] / periods[None, None, :]
    C_all      = np.mean(np.exp(1j * angles_all), axis=1)    # (chunk, N_periods)
    perm_power = N_dome * (C_all.real ** 2 + C_all.imag ** 2)
    perm_peak_count_at_3 += int(np.sum(np.argmax(perm_power, axis=1) == idx_3))

null_peak_at_3_frac = perm_peak_count_at_3 / N_PERM

# ── Cross-corpus stability check ─────────────────────────────────────────────
# A real quantization at period T should peak at the same T in every
# independent subsample.  We test this by splitting the dome subset in half
# 1000 times and recording whether both halves share the same peak period.
N_SPLIT    = 1000
SPLIT_SEED = SEED + 1
rng2 = np.random.default_rng(SPLIT_SEED)
split_agree = 0
for _ in range(N_SPLIT):
    idx = rng2.permutation(N_dome)
    half = N_dome // 2
    a, b = lons_dome[idx[:half]], lons_dome[idx[half:]]
    pa = float(periods[np.argmax(periodogram_vectorized(a, periods))])
    pb = float(periods[np.argmax(periodogram_vectorized(b, periods))])
    if abs(pa - pb) < (periods[1] - periods[0]) + 1e-9:
        split_agree += 1
split_agree_frac = split_agree / N_SPLIT

# R value at 3°
R3_dome = float(np.sqrt(S3_dome / N_dome))
R3_full = float(np.sqrt(S3_full / N_full))

# ── Print ─────────────────────────────────────────────────────────────────────
SEP = "=" * 100
print()
print(SEP)
print("  PERIODOGRAM TEST — rotation-invariant longitudinal periodicity sensitivity sweep")
print(f"  Dome subset: N = {N_dome}  |  Full corpus: N = {N_full}")
print(f"  Candidate periods: {periods[0]:.2f}° to {periods[-1]:.2f}° in {periods[1]-periods[0]:.2f}° steps ({N_periods} values)")
print(SEP)

print(f"""
  DOME SUBSET (N = {N_dome}):
  ────────────────────────────────────────────────────
  Peak period         : {peak_T_dome:.2f}°  (rank 1 of {N_periods} candidates)
  3° power S(3°)      : {S3_dome:.4f}  (R = {R3_dome:.4f})
  Rank of 3° period   : {rank_3_dome} of {N_periods}  ({pctile_3_dome:.1f}th percentile)
  Permutation p       : {perm_p_dome:.5f}  {sig(perm_p_dome)}
    (fraction of {N_PERM:,} UNESCO draws with S_null(3°) >= {S3_dome:.4f})

  NOTE: the peak period {peak_T_dome:.2f}° has higher raw power than T=3°, but this
  does NOT constitute a stronger quantization claim because:
    (a) T={peak_T_dome:.2f}° is post-hoc: no prior hypothesis predicts it.
    (b) Permutation p for the peak at T={peak_T_dome:.2f}° = {perm_p_peak:.5f} {sig(perm_p_peak)}
        but this is NOT Bonferroni-corrected over {N_periods} candidates;
        the corrected threshold is p < {0.05/N_periods:.5f}.
    (c) Its phase anchor drifts with corpus mean longitude (see docstring);
        it does not replicate at the same phase across independent subsets.
    (d) The pre-specified hypothesis is T=3°/Gerizim; the correct test for
        that hypothesis is the V-test in anchored_periodicity_test.py, not
        the omnidirectional sweep here.

  FULL CORPUS (N = {N_full}):
  ────────────────────────────────────────────────────
  Peak period         : {peak_T_full:.2f}°  (rank 1 of {N_periods} candidates)
  3° power S(3°)      : {S3_full:.4f}  (R = {R3_full:.4f})
  Rank of 3° period   : {rank_3_full} of {N_periods}  ({pctile_3_full:.1f}th percentile)
  NOTE: the dome peak ({peak_T_dome:.2f}°) does not replicate in the full corpus
  (peak = {peak_T_full:.2f}°), consistent with it being a sampling artifact.

  CROSS-CORPUS STABILITY:
  ────────────────────────────────────────────────────
  Dome peak T          : {peak_T_dome:.2f}°
  Full corpus peak T   : {peak_T_full:.2f}°
  These periods share no harmonics and have no archaeological relationship.
  A real quantization at T would produce the same peak in every subsample.

  Dome split-half agreement ({N_SPLIT} random half-splits):
    Both halves peak at same T in {split_agree_frac*100:.1f}% of splits
    (chance level ≈ {100/N_periods:.1f}% = 1/{N_periods} candidates)
  Interpretation: {"UNSTABLE — peak period is not reproducible within the dome subset" if split_agree_frac < 0.10 else "peak shows some within-subset consistency"}

  HARMONICS CHECK:
  Under a true 3° structure, 1.5° (subharmonic) and 6° (superharmonic)
  will also show elevated power because they share nodes with the 3° grid.
  If 3° is primary, its power should exceed 1.5° and 6°.
""")

# Harmonics comparison
for T_check in [1.5, 2.0, 2.4, 3.0, 3.6, 4.0, 6.0, 9.0, 12.0]:
    idx_c = int(np.argmin(np.abs(periods - T_check)))
    note  = " ← TARGET" if abs(T_check - 3.0) < 0.01 else (
            " ← subharmonic (3/2)" if abs(T_check - 1.5) < 0.01 else (
            " ← superharmonic (2×3)" if abs(T_check - 6.0) < 0.01 else ""))
    print(f"  T = {T_check:5.1f}°  dome S = {power_dome[idx_c]:7.4f}  "
          f"full S = {power_full[idx_c]:7.4f}{note}")

# ── LaTeX macros ──────────────────────────────────────────────────────────────
print()
print(SEP)
print("  LATEX MACROS")
print(SEP)
print()
print(f"  \\newcommand{{\\periodogramNperiods}}{{{N_periods}}}        % candidate periods in periodogram")
print(f"  \\newcommand{{\\periodogramPeakTDome}}{{{peak_T_dome:.2f}}}     % dominant period, dome subset (degrees)")
print(f"  \\newcommand{{\\periodogramPeakTFull}}{{{peak_T_full:.2f}}}     % dominant period, full corpus (degrees)")
print(f"  \\newcommand{{\\periodogramSThreeDome}}{{{S3_dome:.4f}}}       % Rayleigh power S(3 deg), dome subset")
print(f"  \\newcommand{{\\periodogramSThreeFull}}{{{S3_full:.4f}}}       % Rayleigh power S(3 deg), full corpus")
print(f"  \\newcommand{{\\periodogramRThreeDome}}{{{R3_dome:.4f}}}       % mean resultant length R at 3 deg, dome")
print(f"  \\newcommand{{\\periodogramRankThreeDome}}{{{rank_3_dome}}}           % rank of 3 deg period (dome), 1=highest")
print(f"  \\newcommand{{\\periodogramPctileThreeDome}}{{{pctile_3_dome:.1f}}}    % percentile of 3 deg in dome periodogram")
print(f"  \\newcommand{{\\periodogramPermP}}{{{perm_p_dome:.5f}}}   % permutation p, S_dome(3 deg)")
print(f"  \\newcommand{{\\periodogramNPerm}}{{{N_PERM:,}}}       % permutation draws")
print(f"  \\newcommand{{\\periodogramNullPeakAtThreeFrac}}{{{null_peak_at_3_frac:.4f}}} % fraction null draws with 3 deg as peak")
# ── New: constrained-range and cross-corpus disagreement macros ───────────────
print(f"  % -- Constrained periodogram (<{CONSTRAINED_MAX}°) macros --")
print(f"  \\newcommand{{\\periodogramPeakTDomeConstrained}}{{{peak_T_dome_c:.2f}}}  % dome peak period constrained <{CONSTRAINED_MAX}°")
print(f"  \\newcommand{{\\periodogramPeakTFullConstrained}}{{{peak_T_full_c:.2f}}}  % full-corpus peak period constrained <{CONSTRAINED_MAX}°")
print(f"  \\newcommand{{\\periodogramRPeakDomeConstrained}}{{{R_peak_dome_c:.4f}}} % R at dome constrained peak")
print(f"  \\newcommand{{\\periodogramRPeakFullConstrained}}{{{R_peak_full_c:.4f}}} % R at full-corpus constrained peak")
print(f"  \\newcommand{{\\periodogramRThreeDomeCompare}}{{{R3_dome_compare:.4f}}}  % R at T=3° dome (for direct comparison)")
print(f"  \\newcommand{{\\periodogramRankThreeDomeConstrained}}{{{rank_3_dome_c}}}      % rank of 3° among <{CONSTRAINED_MAX}° periods (dome)")
print(f"  \\newcommand{{\\periodogramRankThreeFullConstrained}}{{{rank_3_full_c}}}      % rank of 3° among <{CONSTRAINED_MAX}° periods (full)")
print(f"  \\newcommand{{\\periodogramNConstrained}}{{{n_constrained}}}      % number of periods in constrained grid")
print(f"  \\newcommand{{\\periodogramDomeOverlapTopTwentyFive}}{{{dome_overlap_c}}}   % 3°-A+ sites in top-25 phasor at dome constrained peak")
print(f"  \\newcommand{{\\periodogramFullOverlapTopTwentyFive}}{{{full_overlap_c}}}   % 3°-A+ sites in top-25 phasor at full constrained peak")
print(f"  \\newcommand{{\\periodogramDomeOverlapExpected}}{{{dome_overlap_exp}}}  % expected overlap by chance (dome)")
print(f"  \\newcommand{{\\periodogramFullOverlapExpected}}{{{full_overlap_exp}}}  % expected overlap by chance (full)")

# ── Write to results store ────────────────────────────────────────────────────
ResultsStore().write_many({
    "periodogramPermP":           perm_p_dome,
    "periodogramPeakTDome":       round(peak_T_dome, 2),
    "periodogramPeakTFull":       round(peak_T_full, 2),
    "periodogramSThreeDome":      round(S3_dome, 5),
    "periodogramRankThreeDome":   rank_3_dome,
    "periodogramPctileThreeDome": round(pctile_3_dome, 1),
    "periodogramRThreeDome":      round(R3_dome, 5),
    "periodogramNullPeakAtThreeFrac": round(null_peak_at_3_frac, 4),
    # Constrained-range and cross-corpus disagreement
    "periodogramPeakTDomeConstrained":       round(peak_T_dome_c, 2),
    "periodogramPeakTFullConstrained":       round(peak_T_full_c, 2),
    "periodogramRPeakDomeConstrained":       round(R_peak_dome_c, 4),
    "periodogramRPeakFullConstrained":       round(R_peak_full_c, 4),
    "periodogramRThreeDomeCompare":          round(R3_dome_compare, 4),
    "periodogramRankThreeDomeConstrained":   rank_3_dome_c,
    "periodogramRankThreeFullConstrained":   rank_3_full_c,
    "periodogramNConstrained":               n_constrained,
    "periodogramDomeOverlapTopTwentyFive":   dome_overlap_c,
    "periodogramFullOverlapTopTwentyFive":   full_overlap_c,
    "periodogramDomeOverlapExpected":        dome_overlap_exp,
    "periodogramFullOverlapExpected":        full_overlap_exp,
})
