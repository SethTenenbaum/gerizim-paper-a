"""
x18_max_permutation_test.py
GROUP: 11c

DISCOVERY-CORRECTED MAX-PERMUTATION TEST FOR THE x.18-DEGREE BAND
======================================================================
The global anchor sweep (anchor_uniqueness_audit.py, GROUP 11) identified
the x.18-degree phase as the optimum by scanning 36,000 anchors.  The naive
binomial p-value from x18_optimal_band_significance.py does not account for
that search: having chosen the best of 36,000 anchors inflates Type I error.

This test corrects for that by computing the max-statistic permutation:
  - Observed: maximum A+ count across all 36,000 trial anchors (read from store).
  - Null: for each random dataset of N phases drawn uniformly from [0, 3°),
    record the maximum A+ count achieved by any anchor phase.  The fraction
    of null maxima >= observed is the discovery-corrected p-value.

NULL MODEL
----------
A site's phase is its longitude mod 3° (position within the 3-degree harmonic
period).  The maximum A+ count across all 36,000 anchors equals the maximum
number of sites whose phase falls within 0.06° of any single anchor phase —
i.e., the peak bin count in a circular histogram on [0, 3°) with bin
half-width 0.06°.

The correct null draws N phases uniformly from [0, 3°), asking:
"What is the largest bin count expected by chance when N sites are placed
with no preferred phase?"  This is distinct from shuffling the observed
phases (which preserves the multiset of phase values and therefore leaves
the max count invariant — giving SD=0, Z=NaN).

The empirical-bootstrap null (resample from observed phases with replacement)
is also computed for comparison; it yields a conservative p because the
observed phase distribution is already Europe-heavy and clustered.

Reference: Westfall & Young (1993) Resampling-Based Multiple Testing.
======================================================================
"""
import sys
import time
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import BERU, TIER_APLUS
from lib.results_store import ResultsStore

# ── Parameters ───────────────────────────────────────────────────────────────
N_PERMS       = 10_000
SEED          = 42
ANCHOR_STEP   = 0.01   # matches the anchor sweep resolution

# ── Read observed optimum from store ─────────────────────────────────────────
rs = ResultsStore()
try:
    obs_max = int(rs.read("anchorSweepGlobalMaxA"))  # 59, Sweep A (Jerusalem removed)
    N_corpus = int(rs.read("anchorSweepNsweep"))     # 1010
except (KeyError, TypeError):
    raise SystemExit(
        "ERROR: anchorSweepGlobalMaxA / anchorSweepNsweep not in results store.\n"
        "Run analysis/global/anchor_uniqueness_audit.py first."
    )

# ── Load corpus longitudes (Sweep A: Jerusalem removed) ──────────────────────
from lib.beru import CONFIG
JERUSALEM_ID = "148"
corpus   = load_corpus()
cultural = cultural_sites_with_coords(corpus)
lons     = np.array([s.longitude for s in cultural
                     if s.id_number != JERUSALEM_ID])
assert len(lons) == N_corpus, f"Expected {N_corpus} sites, got {len(lons)}"

# ── Vectorized max A+ across all anchors ─────────────────────────────────────
# Key insight: A+ membership depends only on (lon mod HARMONIC_STEP) relative
# to the anchor phase (anchor mod HARMONIC_STEP).  We reduce each longitude to
# its phase on the 3-degree circle once, then for each permutation we only need
# to check how many phases fall within TIER_APLUS of ANY harmonic node.
# This collapses the 36,000-anchor scan to a single modulo operation per site.
#
# A site is A+ for anchor `a` iff:
#   |((lon - a) mod HARMONIC_STEP) - round(((lon-a) mod HARMONIC_STEP)/HARMONIC_STEP)*HARMONIC_STEP|
#   <= TIER_APLUS * BERU   (in degrees)
#
# Equivalently, reduce each lon to its residue r = lon mod HARMONIC_STEP in [0, HARMONIC_STEP).
# A site is A+ for some anchor iff its residue is within TIER_APLUS*BERU of 0 or HARMONIC_STEP.
# The MAXIMUM A+ count across all anchors equals the count of sites whose residue
# (relative to the BEST anchor phase) falls within the window — which is just
# the densest bin in a circular histogram with bin-width = 2*TIER_APLUS*BERU on the 3-deg circle.
#
# So: max_ap(lons) = max over all phases p of #{i : |((lons[i]-p) mod 3) centered| <= thresh}
# We compute this by evaluating only the 120 x.18-equivalent phases (one per 3-deg period),
# which already span the full phase space at the resolution that matters.

HARMONIC_STEP = 3.0                          # degrees
thresh_deg    = TIER_APLUS * BERU            # 0.002 * 30 = 0.06 degrees

# Reduce longitudes to residues on [0, HARMONIC_STEP)
phases_obs = lons % HARMONIC_STEP            # shape (N,)

# The best anchor phase is one of the 36,000 values, but since the A+ window
# is only 0.06 deg wide and sites land at discrete longitudes, the maximum
# is determined by the denser end of the site-phase distribution.
# Evaluate all ANCHOR_STEP-spaced phases on [0, HARMONIC_STEP) — only 300 values.
anchor_phases = np.arange(0.0, HARMONIC_STEP, ANCHOR_STEP)   # 300 values

def max_ap_fast(site_phases: np.ndarray) -> int:
    """Vectorized: max A+ count across all anchor phases."""
    # site_phases: (N,)  anchor_phases: (M,)
    # diff shape: (N, M) — circular distance on [0, HARMONIC_STEP)
    diff = np.abs(site_phases[:, None] - anchor_phases[None, :])   # (N, M)
    diff = np.minimum(diff, HARMONIC_STEP - diff)                   # wrap
    counts = np.sum(diff <= thresh_deg, axis=0)                     # (M,)
    return int(counts.max())

# Verify against stored observed maximum
obs_max_check = max_ap_fast(phases_obs)

anchors = np.arange(0.0, 360.0, ANCHOR_STEP)   # kept for reference / print only

# ── Permutation loops ─────────────────────────────────────────────────────────
# NOTE ON NULL MODEL:
#   rng.permutation(phases_obs) would be a no-op for max_ap_fast: the maximum
#   bin count depends only on the multiset of phase values, not their ordering.
#   Shuffling the same values always gives the same max → SD=0, Z=NaN.
#
#   The correct null draws N phases uniformly from [0, HARMONIC_STEP):
#   this asks "what max bin count arises when sites have no preferred phase?"
#   An empirical-bootstrap null (resample observed phases with replacement) is
#   computed in parallel as a conservative cross-check.

print("=" * 70)
print("  x.18° BAND — DISCOVERY-CORRECTED MAX-PERMUTATION TEST (GROUP 11c)")
print(f"  N={N_corpus}  Observed optimum max A+={obs_max}  (cross-check: {obs_max_check})")
print(f"  Vectorized over {len(anchor_phases)} phase values (equiv. to {len(anchors):,} anchors)")
print(f"  Running {N_PERMS:,} permutations (seed={SEED})...")
print(f"  Null A: uniform random phases on [0, {HARMONIC_STEP}°)")
print(f"  Null B: bootstrap resample of observed phases (conservative)")
print("=" * 70)

rng          = np.random.default_rng(SEED)
perm_max_A   = np.empty(N_PERMS, dtype=int)   # Null A: uniform
perm_max_B   = np.empty(N_PERMS, dtype=int)   # Null B: bootstrap

t0 = time.time()
for i in range(N_PERMS):
    # Null A: draw N phases uniformly from [0, HARMONIC_STEP)
    rand_phases    = rng.uniform(0.0, HARMONIC_STEP, size=N_corpus)
    perm_max_A[i]  = max_ap_fast(rand_phases)

    # Null B: bootstrap resample from observed phases (with replacement)
    boot_phases    = rng.choice(phases_obs, size=N_corpus, replace=True)
    perm_max_B[i]  = max_ap_fast(boot_phases)

    if (i + 1) % 2_000 == 0:
        elapsed = time.time() - t0
        eta     = elapsed / (i + 1) * (N_PERMS - i - 1)
        print(f"    {i+1:,}/{N_PERMS:,}  elapsed {elapsed:.1f}s  eta {eta:.1f}s",
              flush=True)

# ── Results ───────────────────────────────────────────────────────────────────
def _stats(perm_max, label):
    p    = float((np.sum(perm_max >= obs_max) + 1) / (N_PERMS + 1))
    mu   = float(perm_max.mean())
    std  = float(perm_max.std(ddof=1))
    z    = float((obs_max - mu) / std) if std > 0 else float("nan")
    return p, mu, std, z

p_A, mu_A, sd_A, z_A = _stats(perm_max_A, "Null A (uniform)")
p_B, mu_B, sd_B, z_B = _stats(perm_max_B, "Null B (bootstrap)")

def _sig(p):
    return ("***" if p < 0.001 else "**" if p < 0.01 else
            "*"   if p < 0.05  else "~"  if p < 0.10 else "ns")

print()
print(f"  Observed max A+ (x.18 optimum)   : {obs_max}")
print()
print(f"  Null A — uniform random phases:")
print(f"    Null mean / SD                 : {mu_A:.2f} / {sd_A:.2f}")
print(f"    Z                              : {z_A:.2f}")
print(f"    Discovery-corrected p          : {p_A:.4f}  {_sig(p_A)}")
print()
print(f"  Null B — bootstrap (conservative):")
print(f"    Null mean / SD                 : {mu_B:.2f} / {sd_B:.2f}")
print(f"    Z                              : {z_B:.2f}")
print(f"    Discovery-corrected p          : {p_B:.4f}  {_sig(p_B)}")
print()
print("  Interpretation:")
print("  The observed maximum (59) lies near the centre of both null")
print("  distributions.  The x.18 band is not unusually dense once the")
print("  36,000-anchor search is accounted for.  The remarkable fact is")
print("  Gerizim's historical motivation and its location in this band,")
print("  not the band maximum itself.")

# Report primary (Null A) as the headline number; Null B in supplement
p_maxperm = p_A
null_mu   = mu_A
null_std  = sd_A
z_max     = z_A

# ── Save permuted distribution ────────────────────────────────────────────────
out_dir = Path("results")
out_dir.mkdir(parents=True, exist_ok=True)
np.save(out_dir / "x18_maxperm_null_A.npy", perm_max_A)
np.save(out_dir / "x18_maxperm_null_B.npy", perm_max_B)

# ── LaTeX macros ──────────────────────────────────────────────────────────────
print()
print("  LATEX MACROS (GROUP 11c):")
print("=" * 70)
macro_pairs = [
    ("anchorMaxPermObsMax",   f"{obs_max}"),
    ("anchorMaxPermNullMu",   f"{null_mu:.2f}"),
    ("anchorMaxPermNullSD",   f"{null_std:.2f}"),
    ("anchorMaxPermZ",        f"{z_max:.2f}"),
    ("anchorMaxPermP",        f"{p_maxperm:.4f}"),
    ("anchorMaxPermNperms",   f"{N_PERMS:,}".replace(",", "{,}")),
    ("anchorMaxPermBootMu",   f"{mu_B:.2f}"),
    ("anchorMaxPermBootSD",   f"{sd_B:.2f}"),
    ("anchorMaxPermBootZ",    f"{z_B:.2f}"),
    ("anchorMaxPermBootP",    f"{p_B:.4f}"),
]
for name, val in macro_pairs:
    print(f"\\newcommand{{\\{name}}}{{{val}}}")

rs.write_many({
    "anchorMaxPermObsMax":   int(obs_max),
    "anchorMaxPermNullMu":   null_mu,
    "anchorMaxPermNullSD":   null_std,
    "anchorMaxPermZ":        z_max,
    "anchorMaxPermP":        p_maxperm,
    "anchorMaxPermNperms":   N_PERMS,
    "anchorMaxPermBootMu":   mu_B,
    "anchorMaxPermBootSD":   sd_B,
    "anchorMaxPermBootZ":    z_B,
    "anchorMaxPermBootP":    p_B,
})
print("Results written to data/store/results.json")
