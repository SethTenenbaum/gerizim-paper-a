"""Temporary runner for x18 periodicity formal test."""
import sys, numpy as np, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, TIER_APLUS
from lib.results_store import ResultsStore

HARMONIC_STEP_DEG = 3.0
N_PERMS = 10_000
PHASE_SCAN_STEP = 0.01

corpus = load_corpus()
cultural = cultural_sites_with_coords(corpus)
lons = np.array([s.longitude for s in cultural])
N = len(lons)

def beru_dev(lon):
    arc = min(abs(lon - GERIZIM), 360.0 - abs(lon - GERIZIM))
    bv = arc / BERU
    return abs(bv - round(bv * 10) / 10)

is_ap = np.array([beru_dev(lo) <= TIER_APLUS for lo in lons])
ap_lons = lons[is_ap]
N_ap = int(is_ap.sum())

def rayleigh_R(lon_array):
    phases = lon_array % HARMONIC_STEP_DEG
    angles = 2 * np.pi * phases / HARMONIC_STEP_DEG
    return float(abs(np.mean(np.exp(1j * angles))))

phases_scan = np.arange(0.0, HARMONIC_STEP_DEG, PHASE_SCAN_STEP)

def phase_tie_frac_vec(lon_arr):
    diff = np.abs(lon_arr[None,:] - phases_scan[:,None])
    arcs = np.minimum(diff, 360.0 - diff)
    bv   = arcs / BERU
    devs = np.abs(bv - np.round(bv * 10) / 10)
    counts = np.sum(devs <= TIER_APLUS, axis=1)
    return float(np.mean(counts == counts.max()))

obs_full_R = rayleigh_R(lons)
obs_ap_R   = rayleigh_R(ap_lons)
obs_tie    = phase_tie_frac_vec(lons)
n_phases   = int(round(HARMONIC_STEP_DEG / PHASE_SCAN_STEP))
n_tie      = int(round(obs_tie * n_phases))

print("=" * 80)
print("  x.18 PERIODICITY — FORMAL STATISTICAL TESTS")
print(f"  N={N}  N_ap={N_ap}  Harmonic step={HARMONIC_STEP_DEG}  Anchor=Gerizim {GERIZIM}E")
print("=" * 80)
print(f"  Obs A (full-corpus Rayleigh R, mod 3): {obs_full_R:.4f}")
print(f"  Obs B (A+-only Rayleigh R, mod 3):     {obs_ap_R:.4f}")
print(f"  Obs C (phase-tie fraction):             {obs_tie:.4f}  ({n_tie}/{n_phases} phases tie)")
print(f"\n  Running {N_PERMS:,} permutations ...")

rng = np.random.default_rng(42)
perm_full = np.zeros(N_PERMS)
perm_ap   = np.zeros(N_PERMS)
perm_tie  = np.zeros(N_PERMS)

t0 = time.time()
for i in range(N_PERMS):
    pl            = rng.permutation(lons)
    perm_full[i]  = rayleigh_R(pl)
    perm_ap[i]    = rayleigh_R(pl[:N_ap])
    perm_tie[i]   = phase_tie_frac_vec(pl)
    if (i + 1) % 2_000 == 0:
        print(f"    ... {i+1}/{N_PERMS}  ({time.time()-t0:.1f}s)", flush=True)

def sig(p):
    return "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "~" if p<0.10 else "ns"

def _z(obs, arr):
    s = arr.std()
    return float((obs - arr.mean()) / s) if s > 0 else float("nan")

p_full = float(np.mean(perm_full >= obs_full_R))
p_ap   = float(np.mean(perm_ap   >= obs_ap_R))
p_tie  = float(np.mean(perm_tie  >= obs_tie))
z_full = _z(obs_full_R, perm_full)
z_ap   = _z(obs_ap_R,   perm_ap)
z_tie  = _z(obs_tie,    perm_tie)

print(f"\n  PERMUTATION RESULTS ({N_PERMS:,} shuffles, seed=42)")
print(f"  {'Test':<45}  {'Obs':>8}  {'Null_mu':>8}  {'Z':>6}  {'p':>8}  Sig")
for label, obs, mu, z, p in [
    ("A: Rayleigh R full corpus (mod 3)",          obs_full_R, perm_full.mean(), z_full, p_full),
    (f"B: Rayleigh R A+ sites N={N_ap} (mod 3)",  obs_ap_R,   perm_ap.mean(),   z_ap,   p_ap),
    ("C: Phase-tie fraction",                       obs_tie,    perm_tie.mean(),  z_tie,  p_tie),
]:
    print(f"  {label:<45}  {obs:>8.4f}  {mu:>8.4f}  {z:>6.2f}  {p:>8.4f}  {sig(p)}")

print(f"\n  LATEX MACROS (GROUP 11b — x.18° periodicity formal test):")
print("=" * 80)
macro_pairs = [
    ("fullRayleighR",     f"{obs_full_R:.4f}"),
    ("fullRayleighZ",     f"{z_full:.2f}"),
    ("fullRayleighPermP", f"{p_full:.4f}"),
    ("rayleighR",         f"{obs_ap_R:.4f}"),
    ("rayleighZ",         f"{z_ap:.2f}"),
    ("rayleighPermP",     f"{p_ap:.4f}"),
    ("phaseTieFrac",      f"{obs_tie:.4f}"),
    ("phaseTieZ",         f"{z_tie:.2f}"),
    ("phaseTiePermP",     f"{p_tie:.4f}"),
]
for name, val in macro_pairs:
    print(f"\\newcommand{{\\{name}}}{{{val}}}")

ResultsStore().write_many({
    "fullRayleighR":     float(obs_full_R),
    "fullRayleighZ":     float(z_full),
    "fullRayleighPermP": float(p_full),
    "rayleighR":         float(obs_ap_R),
    "rayleighZ":         float(z_ap),
    "rayleighPermP":     float(p_ap),
    "phaseTieFrac":      float(obs_tie),
    "phaseTieZ":         float(z_tie),
    "phaseTiePermP":     float(p_tie),
})
print("Results written to data/store/results.json")
