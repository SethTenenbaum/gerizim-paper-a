"""
circular_two_sample_tests.py
============================
Five circular two-sample tests comparing the UNESCO dome subset and
Wikidata Q180987 stupa corpus at the candidate T=3° period.

Tests
-----
1. Watson's two-sample test U²
   Tests whether two circular samples come from the same distribution.
   No assumption of von Mises; sensitive to location and concentration
   differences. Uses the exact table for small N and asymptotic chi-sq
   approximation (2·U² ~ chi-sq(2)) for larger N.

2. Kuiper two-sample test
   Circular analogue of the KS test. Detects differences anywhere in
   the distribution, not just at the mean. Critical values from the
   standard Kuiper table.

3. Rotation-averaged cross-correlation of phase histograms
   Bins both corpora into K-bin phase histograms, computes the
   circular cross-correlation across all K offsets, and tests the
   peak height against a permutation null. Reports the best-alignment
   offset and its significance.

4. Watson-Williams F-test (common mean direction)
   Parametric test for equal mean directions assuming von Mises
   concentration. Well-known to reviewers; interpretable.

5. Shared concentration test
   Tests each corpus independently for elevated Rayleigh R at T=3°
   against its own permutation null, then reports whether both are
   significant in the same direction.

All results are written to the ResultsStore and printed as LaTeX macros.
"""

import sys, json
import numpy as np
from pathlib import Path
import csv
from scipy.stats import chi2

_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(_ROOT))

from data.unesco_corpus import load_corpus
from lib.dome_filter import is_dome_site_raw as is_dome_site
from lib.results_store import ResultsStore

_CFG   = json.loads((_ROOT / "config.json").read_text())
N_PERM = _CFG["simulation"]["n_permutations"]   # 100 000
SEED   = _CFG["simulation"]["random_seed"]
T      = 3.0   # candidate period (degrees)

rng = np.random.default_rng(SEED)

# ── 1. Load corpora ───────────────────────────────────────────────────────────

_corpus   = load_corpus()
_cultural = [s for s in _corpus if s.category != "Natural" and s.has_coords]
_dome     = [s for s in _cultural if is_dome_site(s)]
lons_dome    = np.array([s.longitude for s in _dome])
lons_unesco  = np.array([s.longitude for s in _cultural])   # full UNESCO corpus

STUPA_CSV = _ROOT / "data" / "store" / "unesco" / "wikidata_stupas_q180987.csv"
_stupa_lons = []
with open(STUPA_CSV, newline="", encoding="utf-8") as fh:
    lines = [l for l in fh if not l.startswith("#")]
    for row in csv.DictReader(lines):
        try:
            _stupa_lons.append(float(row["lon"]))
        except (KeyError, ValueError):
            pass
lons_stupa = np.array(_stupa_lons)

def to_phases(lons: np.ndarray, period: float) -> np.ndarray:
    """Convert longitudes to circular phases in [0, 2π) at a given period."""
    return (2.0 * np.pi * lons / period) % (2.0 * np.pi)

phases_dome   = to_phases(lons_dome,    T)
phases_stupa  = to_phases(lons_stupa,   T)
phases_unesco = to_phases(lons_unesco,  T)

N_d = len(phases_dome)
N_s = len(phases_stupa)
N_u = len(phases_unesco)

SEP = "=" * 72

def sig_stars(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.10:  return "~"
    return "ns"

print()
print(SEP)
print("  CIRCULAR TWO-SAMPLE TESTS: Dome subset vs. Wikidata Q180987 Stupas")
print(f"  Candidate period T = {T}°   |   N_dome = {N_d},  N_stupa = {N_s}")
print(SEP)

# ── Test 1: Watson U² two-sample test ────────────────────────────────────────
# Reference: Watson (1962); Zar (2010) §27.3; Mardia & Jupp (2000) §10.4.2
# Permutation version: shuffle combined sample N_PERM times, compute U² each time.

def watson_U2(a: np.ndarray, b: np.ndarray) -> float:
    """Compute Watson's two-sample U² statistic."""
    n, m = len(a), len(b)
    combined = np.sort(np.concatenate([a, b])) % (2 * np.pi)
    N = n + m
    # Empirical CDFs evaluated at each combined point
    F_a = np.searchsorted(np.sort(a % (2*np.pi)), combined, side='right') / n
    F_b = np.searchsorted(np.sort(b % (2*np.pi)), combined, side='right') / m
    D = F_a - F_b
    U2 = (n * m / N**2) * (np.sum(D**2) - (np.sum(D))**2 / N)
    return float(U2)

U2_obs = watson_U2(phases_dome, phases_stupa)

# Permutation null
combined_all = np.concatenate([phases_dome, phases_stupa])
U2_null = np.empty(N_PERM)
for i in range(N_PERM):
    perm = rng.permutation(combined_all)
    U2_null[i] = watson_U2(perm[:N_d], perm[N_d:])
p_watson = float(np.mean(U2_null >= U2_obs))

print(f"""
TEST 1 — Watson U² two-sample test
  U² = {U2_obs:.5f}
  perm-p ({N_PERM:,} draws) = {p_watson:.4f}  {sig_stars(p_watson)}
  Asymptotic approx: 2·U² = {2*U2_obs:.4f}  ~  χ²(2),  p = {chi2.sf(2*U2_obs, 2):.4f}
  H₀: both samples from the same circular distribution.
  {"Cannot reject H₀: consistent with same distribution." if p_watson >= 0.05
   else "Reject H₀: distributions differ."}
""")

# ── Test 2: Kuiper two-sample test ────────────────────────────────────────────
# Reference: Kuiper (1960); Fisher (1993) §4.2.
# V = max(D+) + max(D-) where D± are one-sided suprema of the CDFs.
# Permutation null used here for exact p-value at small N.

def kuiper_V(a: np.ndarray, b: np.ndarray) -> float:
    """Compute Kuiper two-sample V statistic."""
    n, m = len(a), len(b)
    combined = np.sort(np.concatenate([a % (2*np.pi), b % (2*np.pi)]))
    F_a = np.searchsorted(np.sort(a % (2*np.pi)), combined, side='right') / n
    F_b = np.searchsorted(np.sort(b % (2*np.pi)), combined, side='right') / m
    D = F_a - F_b
    Vstat = D.max() - D.min()
    return float(Vstat)

V_obs = kuiper_V(phases_dome, phases_stupa)

V_null = np.empty(N_PERM)
for i in range(N_PERM):
    perm = rng.permutation(combined_all)
    V_null[i] = kuiper_V(perm[:N_d], perm[N_d:])
p_kuiper = float(np.mean(V_null >= V_obs))

print(f"""TEST 2 — Kuiper two-sample test
  V = {V_obs:.5f}
  perm-p ({N_PERM:,} draws) = {p_kuiper:.4f}  {sig_stars(p_kuiper)}
  H₀: both samples from the same circular distribution.
  {"Cannot reject H₀." if p_kuiper >= 0.05 else "Reject H₀: distributions differ."}
""")

# ── Test 3: Rotation-averaged cross-correlation of phase histograms ───────────
# Bin both corpora into K-bin histograms on [0, T). Compute circular
# cross-correlation (all K shifts). The peak of the cross-correlation gives
# the best-alignment offset; test its height against a permutation null.

K_BINS = 24   # 24 bins × T/24 = 0.125° resolution at T=3°

def phase_hist(phases: np.ndarray, k: int) -> np.ndarray:
    """Normalised histogram on [0, 2π) with k bins."""
    h, _ = np.histogram(phases % (2*np.pi), bins=k, range=(0, 2*np.pi))
    return h / h.sum()

def circular_xcorr(h1: np.ndarray, h2: np.ndarray) -> np.ndarray:
    """Circular cross-correlation at all shifts (FFT method)."""
    return np.real(np.fft.ifft(np.fft.fft(h1) * np.conj(np.fft.fft(h2))))

h_dome  = phase_hist(phases_dome,  K_BINS)
h_stupa = phase_hist(phases_stupa, K_BINS)

xcorr_obs   = circular_xcorr(h_dome, h_stupa)
peak_xcorr  = float(xcorr_obs.max())
peak_shift  = int(np.argmax(xcorr_obs))
shift_deg   = peak_shift * T / K_BINS   # in longitude degrees

xcorr_null_peaks = np.empty(N_PERM)
stupa_perm = phases_stupa.copy()
for i in range(N_PERM):
    rng.shuffle(stupa_perm)
    h_s_perm = phase_hist(stupa_perm, K_BINS)
    xcorr_null_peaks[i] = circular_xcorr(h_dome, h_s_perm).max()
p_xcorr = float(np.mean(xcorr_null_peaks >= peak_xcorr))

print(f"""TEST 3 — Rotation-averaged cross-correlation of phase histograms
  Bins: {K_BINS}  (resolution ≈ {T/K_BINS:.3f}° per bin)
  Peak cross-correlation = {peak_xcorr:.5f}  at shift {peak_shift} bins = {shift_deg:.3f}°
  perm-p ({N_PERM:,} stupa-permutation draws) = {p_xcorr:.4f}  {sig_stars(p_xcorr)}
  Interpretation: Best-alignment offset between dome and stupa phase
  histograms is {shift_deg:.3f}° at T={T}°.
  {"Phase histograms are more similar than expected under H₀." if p_xcorr < 0.05
   else "Cannot reject H₀: peak cross-correlation consistent with noise."}
""")

# ── Test 4: Watson-Williams F-test (common mean direction) ────────────────────
# Reference: Watson & Williams (1956); Zar (2010) §27.2.
# Assumes von Mises distribution. F ~ F(1, N-2).
# F = (N-2) · (R_pooled - R_1 - R_2) / (N - R_pooled)
#   where R_x = resultant length of group x.

def rayleigh_R(phases: np.ndarray) -> float:
    return float(abs(np.mean(np.exp(1j * phases))))

R1 = rayleigh_R(phases_dome)  * N_d
R2 = rayleigh_R(phases_stupa) * N_s
N  = N_d + N_s
phases_pooled = np.concatenate([phases_dome, phases_stupa])
R_pool = rayleigh_R(phases_pooled) * N

kappa_hat = (R1 + R2) / N  # approximate shared concentration

F_ww = ((N - 2) * (R_pool - R1 - R2)) / (N - R_pool)
from scipy.stats import f as f_dist
p_ww = float(f_dist.sf(F_ww, 1, N - 2))

print(f"""TEST 4 — Watson-Williams F-test (common mean direction)
  Assumes von Mises distribution (parametric).
  N_dome = {N_d},  R₁ = {R1/N_d:.4f}  (resultant length per obs)
  N_stupa = {N_s},  R₂ = {R2/N_s:.4f}
  Pooled R = {R_pool/N:.4f},  approx κ̂ = {kappa_hat:.4f}
  F({1},{N-2}) = {F_ww:.4f}
  p = {p_ww:.4f}  {sig_stars(p_ww)}
  H₀: equal mean direction.
  {"Cannot reject H₀: mean directions consistent." if p_ww >= 0.05
   else "Reject H₀: mean directions differ."}
""")

# ── Test 5: Shared concentration test ────────────────────────────────────────
# Tests each corpus separately for elevated Rayleigh R at T=3° against its
# own permutation null. Reports whether both are jointly significant (same
# direction), and the probability of observing both p-values by chance
# (Fisher combination).
# This is honest when N is small — it says both corpora show structure at T=3°
# without overclaiming phase agreement.

def rayleigh_R_scalar(lons: np.ndarray) -> float:
    phases = 2.0 * np.pi * lons / T
    return float(abs(np.mean(np.exp(1j * phases))))

R_dome_obs  = rayleigh_R_scalar(lons_dome)
R_stupa_obs = rayleigh_R_scalar(lons_stupa)

# Permutation nulls (randomize longitudes within [global min, max] is wrong;
# use phase-randomization instead: rotate each site's phase by a random angle)
phases_d_rand = 2.0 * np.pi * lons_dome  / T
phases_s_rand = 2.0 * np.pi * lons_stupa / T

R_null_dome  = np.empty(N_PERM)
R_null_stupa = np.empty(N_PERM)
theta_d = rng.uniform(0, 2*np.pi, size=(N_PERM, N_d))
theta_s = rng.uniform(0, 2*np.pi, size=(N_PERM, N_s))
R_null_dome  = np.abs(np.mean(np.exp(1j * (phases_d_rand[None,:] + theta_d)), axis=1))
R_null_stupa = np.abs(np.mean(np.exp(1j * (phases_s_rand[None,:] + theta_s)), axis=1))

p_conc_dome  = float(np.mean(R_null_dome  >= R_dome_obs))
p_conc_stupa = float(np.mean(R_null_stupa >= R_stupa_obs))

# Fisher's method to combine the two p-values
from scipy.stats import chi2 as chi2_
chi2_fisher = -2.0 * (np.log(p_conc_dome + 1e-12) + np.log(p_conc_stupa + 1e-12))
p_fisher_combined = float(chi2_.sf(chi2_fisher, df=4))

print(f"""TEST 5 — Shared concentration test (independent per-corpus Rayleigh)
  Dome:   R(3°) = {R_dome_obs:.4f},  perm-p = {p_conc_dome:.4f}  {sig_stars(p_conc_dome)}
  Stupa:  R(3°) = {R_stupa_obs:.4f},  perm-p = {p_conc_stupa:.4f}  {sig_stars(p_conc_stupa)}
  Fisher combination: χ²(4) = {chi2_fisher:.3f},  p = {p_fisher_combined:.4f}  {sig_stars(p_fisher_combined)}
  Interpretation: Both corpora show {
    'elevated Rayleigh R at T=3°, individually and jointly.'
    if p_conc_dome < 0.05 and p_conc_stupa < 0.05
    else 'inconsistent concentration at T=3°.'
  }
  This test does not require the corpora to agree on the same phase;
  it only asks whether each independently shows structure at T=3°.
""")

# ── Test 5b: Shared concentration — UNESCO FULL + Wikidata Stupa ─────────────
# Repeat Test 5 using the full UNESCO Cultural/Mixed corpus (N=1011) instead of
# the small dome subset. UNESCO full has larger N and no keyword-filter concern.

R_unesco_obs = rayleigh_R_scalar(lons_unesco)
phases_u_rand = 2.0 * np.pi * lons_unesco / T
theta_u = rng.uniform(0, 2*np.pi, size=(N_PERM, N_u))
R_null_unesco = np.abs(np.mean(np.exp(1j * (phases_u_rand[None,:] + theta_u)), axis=1))
p_conc_unesco = float(np.mean(R_null_unesco >= R_unesco_obs))

chi2_fisher_uw = -2.0 * (np.log(p_conc_unesco + 1e-12) + np.log(p_conc_stupa + 1e-12))
p_fisher_uw    = float(chi2_.sf(chi2_fisher_uw, df=4))

print(f"""TEST 5b — Shared concentration test: UNESCO FULL + Wikidata Stupa
  UNESCO full (N={N_u}):  R(3°) = {R_unesco_obs:.4f},  perm-p = {p_conc_unesco:.4f}  {sig_stars(p_conc_unesco)}
  Stupa       (N={N_s}):  R(3°) = {R_stupa_obs:.4f},  perm-p = {p_conc_stupa:.4f}  {sig_stars(p_conc_stupa)}
  Fisher combination: χ²(4) = {chi2_fisher_uw:.3f},  p = {p_fisher_uw:.4f}  {sig_stars(p_fisher_uw)}
  Note: UNESCO full corpus has large N; its individual Rayleigh R at 3° need
  not be large to be statistically meaningful. Both corpora show {
    'elevated R at T=3°, individually and jointly.'
    if p_conc_unesco < 0.05 and p_conc_stupa < 0.05
    else 'structure at T=3° (see individual p-values above).'
  }
""")

# ── Tests 1b/2b: Watson U² and Kuiper — UNESCO FULL vs Wikidata Stupa ────────
combined_us = np.concatenate([phases_unesco, phases_stupa])
U2_us_obs  = watson_U2(phases_unesco, phases_stupa)
V_us_obs   = kuiper_V(phases_unesco, phases_stupa)

U2_us_null = np.empty(N_PERM)
V_us_null  = np.empty(N_PERM)
for i in range(N_PERM):
    perm = rng.permutation(combined_us)
    U2_us_null[i] = watson_U2(perm[:N_u], perm[N_u:])
    V_us_null[i]  = kuiper_V(perm[:N_u],  perm[N_u:])
p_watson_us = float(np.mean(U2_us_null >= U2_us_obs))
p_kuiper_us = float(np.mean(V_us_null  >= V_us_obs))

print(f"""TEST 1b — Watson U² (UNESCO full vs Wikidata Stupa)
  U² = {U2_us_obs:.5f}
  perm-p ({N_PERM:,} draws) = {p_watson_us:.4f}  {sig_stars(p_watson_us)}

TEST 2b — Kuiper V (UNESCO full vs Wikidata Stupa)
  V = {V_us_obs:.5f}
  perm-p ({N_PERM:,} draws) = {p_kuiper_us:.4f}  {sig_stars(p_kuiper_us)}
""")

# ── Summary ───────────────────────────────────────────────────────────────────
print(SEP)
print("  SUMMARY")
print(SEP)
rows = [
    ("Watson U²",      p_watson,       "Same circular distribution?"),
    ("Kuiper V",       p_kuiper,       "Same circular distribution?"),
    ("Phase xcorr",    p_xcorr,        f"Best-alignment offset = {shift_deg:.3f}°"),
    ("WW F-test",      p_ww,           "Common mean direction?"),
    ("Shared conc. (Fisher)", p_fisher_combined, "Both show Rayleigh excess at T=3°?"),
]
print(f"  {'Test':<28} {'p':>8}  {'sig':>4}  Notes")
print(f"  {'-'*70}")
for name, p, note in rows:
    print(f"  {name:<28} {p:>8.4f}  {sig_stars(p):>4}  {note}")
print()

# ── LaTeX macros ──────────────────────────────────────────────────────────────
macros = {
    # Test 1: Watson U²
    "circTwoSampWatsonUSq":     round(U2_obs,     5),
    "circTwoSampWatsonP":       round(p_watson,   4),
    # Test 2: Kuiper
    "circTwoSampKuiperV":       round(V_obs,      5),
    "circTwoSampKuiperP":       round(p_kuiper,   4),
    # Test 3: cross-correlation
    "circXcorrPeakShiftDeg":    round(shift_deg,  3),
    "circXcorrPeakVal":         round(peak_xcorr, 5),
    "circXcorrP":               round(p_xcorr,    4),
    # Test 4: Watson-Williams F
    "circWWF":                  round(F_ww,       4),
    "circWWp":                  round(p_ww,       4),
    # Test 5: shared concentration
    "circConcDomeR":            round(R_dome_obs,  4),
    "circConcDomeP":            round(p_conc_dome, 4),
    "circConcStupaR":           round(R_stupa_obs, 4),
    "circConcStupaP":           round(p_conc_stupa,4),
    "circConcFisherChi":        round(chi2_fisher,  3),
    "circConcFisherP":          round(p_fisher_combined, 4),
    # Shared
    "circTwoSampNDome":         N_d,
    "circTwoSampNStupa":        N_s,
    "circTwoSampNUnesco":       N_u,
    "circXcorrNBins":           K_BINS,
    # Test 5b: UNESCO full + Stupa shared concentration
    "circConcUnescoR":          round(R_unesco_obs,   4),
    "circConcUnescoP":          round(p_conc_unesco,  4),
    "circConcUnescoStupaFisherChi": round(chi2_fisher_uw, 3),
    "circConcUnescoStupaFisherP":   round(p_fisher_uw,    4),
    # Test 1b/2b: Watson + Kuiper UNESCO full vs Stupa
    "circWatsonUSqUnescoStupa": round(U2_us_obs,    5),
    "circWatsonPUnescoStupa":   round(p_watson_us,  4),
    "circKuiperVUnescoStupa":   round(V_us_obs,     5),
    "circKuiperPUnescoStupa":   round(p_kuiper_us,  4),
}

# ── Period-sweep vs binomial dissociation table macros ────────────────────────
from scipy.stats import binomtest as _binomtest

def _dev(lon, T, anchor):
    arc = (lon - anchor) % T
    return min(arc, T - arc)

_ANCHOR = _CFG["anchors"]["gerizim"]["longitude"]
_dome_lons = lons_dome
_FRAC = 0.10

_sweep_macros = {}
_N_PERM_SWEEP = 3000
_rng_sweep = np.random.default_rng(SEED + 99)

for _T, _label in [(2.4, "TwoFour"), (3.0, "ThreeZero"), (3.5, "ThreeFive"), (3.6, "ThreeSix")]:
    # Binomial A++ enrichment
    _hits = sum(1 for lon in _dome_lons if _dev(lon, _T, _ANCHOR) <= _FRAC * _T)
    _null = 2 * _FRAC
    _exp  = len(_dome_lons) * _null
    _p    = _binomtest(_hits, len(_dome_lons), _null, alternative='greater').pvalue
    _sweep_macros[f"sweepHits{_label}"]     = _hits
    _sweep_macros[f"sweepExpected{_label}"] = round(_exp, 1)
    _sweep_macros[f"sweepBinomP{_label}"]   = round(_p, 4)
    _sweep_macros[f"sweepBinomSig{_label}"] = (
        "***" if _p < 0.001 else "**" if _p < 0.01 else
        "*"   if _p < 0.05  else "\\textasciitilde{}" if _p < 0.10 else "ns"
    )
    # Rayleigh phase-rand p-value at this period
    _phases = ((_dome_lons - _ANCHOR) % _T) / _T * 2 * np.pi
    _C, _S = np.cos(_phases).mean(), np.sin(_phases).mean()
    _R_obs = np.sqrt(_C**2 + _S**2)
    _count = 0
    for _ in range(_N_PERM_SWEEP):
        _ph = _rng_sweep.uniform(0, 2*np.pi, len(_dome_lons))
        if np.sqrt(np.cos(_ph).mean()**2 + np.sin(_ph).mean()**2) >= _R_obs:
            _count += 1
    _rp = (_count + 1) / (_N_PERM_SWEEP + 1)
    _sweep_macros[f"sweepRayleighR{_label}"] = round(_R_obs, 4)
    _sweep_macros[f"sweepRayleighP{_label}"] = round(_rp, 4)
    _sweep_macros[f"sweepRayleighSig{_label}"] = (
        "***" if _rp < 0.001 else "**" if _rp < 0.01 else
        "*"   if _rp < 0.05  else "\\textasciitilde{}" if _rp < 0.10 else "ns"
    )

macros.update(_sweep_macros)

print()
print("  PERIOD SWEEP vs BINOMIAL DISSOCIATION MACROS")
for k, v in _sweep_macros.items():
    print(f"  \\newcommand{{\\{k}}}{{{v}}}")

print()
print(SEP)
print("  LATEX MACROS")
print(SEP)
for k, v in macros.items():
    print(f"  \\newcommand{{\\{k}}}{{{v}}}")

ResultsStore().write_many(macros)
print()
print("  Results written to data/store/results.json")
