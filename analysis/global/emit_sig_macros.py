"""
emit_sig_macros.py
==================
Reads every p-value entry from the ResultsStore and emits a companion
``\\macroNameSig`` LaTeX macro containing the significance stars or 'ns'.

This script must run AFTER all analysis scripts have populated the store.
It is called by reproduce_all_macros.sh as a final step.

Convention
----------
  * p < 0.001  →  ***
  * p < 0.01   →  **
  * p < 0.05   →  *
  * p < 0.10   →  $^\dagger$  (marginal)
  * otherwise  →  ns

Keys that are p-values are identified by:
  1. Explicit whitelist (P_VALUE_KEYS below), OR
  2. Key starts with a lowercase 'p' followed by an uppercase letter
     (camelCase convention used in this pipeline, e.g. pCircAp, pEvoAp).

Run from repo root:
    python3 analysis/global/emit_sig_macros.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from lib.results_store import ResultsStore

# ── Significance helper ───────────────────────────────────────────────────────

def sig(p: float) -> str:
    """Return LaTeX-safe significance label for p."""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    if p < 0.10:
        return r"$^\dagger$"
    return "ns"


# ── Explicit whitelist of p-value keys that don't match the camelCase heuristic
# (e.g. keys starting with 'full', 'anchor', 'sim', 'geo', 'cluster', 'stupa')
EXPLICIT_P_KEYS = {
    # Bonferroni-adjusted
    "pAdjTestOne", "pAdjTestTwo", "pAdjTestThree", "pAdjTestTwoB",
    "pAdjTestFour",  # in case Test 4 macro is added
    # Cluster asymmetry
    "clusterApBinom", "clusterMWp", "clusterPermP", "clusterHarmonicMWp",
    # Dome / spherical raw sweep
    "pCircAp", "pCircA", "pCircApFisher", "pCircAFisher", "pCircChi",
    # Evolution — Fisher and binomial
    "pEvoAp", "pEvoApFisher",
    "pEvoMoundFisher", "pEvoStupaFisher", "pEvoDomeFisher",
    "pEvoAFisher", "pEvoStupaAFisher", "pEvoDomeAFisher", "pEvoMoundAFisher",
    "pEvoApValidated",
    # Religion
    "pReligUnion", "pReligChiSq", "pBudJoint", "pHinJoint",
    # Anchor comparisons
    "pMeccaAnchor", "pMeruAnchor", "pKailashAnchor",
    # Simulation null models
    "simDomePermP", "simDomeBootP", "simKDEP",
    "simCanonPermP", "simPreTwoKPermP", "simPostTwoKPermP",
    # Geographic concentration
    "geoNullDomeBootP", "geoNullDomeRestrictedP",
    "geoNullDomeRestrictedPtwo", "geoNullDomeRestrictedPten",
    "stupaGeoBootP", "stupaGeoRestrictedP",
    # Rayleigh / periodicity formal tests
    "fullRayleighPermP", "rayleighPermP",
    "anchorShiftPermP", "anchorMaxPermP",
    "targetedPfull", "targetedPAp",
    # Temporal
    "pCanon", "pPreTwoK", "pModern", "pFirstHalf", "pSecondHalf",
    "pHalfFisher", "pFisherCanonModern", "pCochranThree", "pCochranFive",
    "pPreCEbinom", "pPreCEfisher", "pFoundDateSpearman",
    # Corridor
    "pCorridorBinom", "pCorridorFisher", "pCorridorAppFisher", "pCorridorMW",
    # Spatial independence
    "pNeffQuarter", "pNeffHalf", "blockBootZ_p",
    # CMH
    "pCMH",
    # Region
    "pRegionEmpAll",
    # Founding sites
    "pFoundKwFisher",
    # Wikidata
    "wikiJavaTightApP", "stupaIndiaPakCP",
    "wikiStupaATierBinomP", "wikiHeartlandCbandP", "wikiMyanmarCbandP",
    "wikiJavaNodeATierP", "wikiJavaTightATierP",
    "wikiNonHeartlandATierP", "wikiNonHeartlandApP", "wikiNonHeartlandCbandP",
    # Americas
    "AmericasOneSidedP",
    # FDR q-values for the BH test family
    "FDRqTestOne", "FDRqTestTwo", "FDRqTestThree",
    # Unit sweep / slope tests
    "sweepEightDomeP", "sweepEightFullP",
    "permSlopeCanonBestJoint", "permSlopeCanonRank",
    # Cluster harmonic string (has embedded p-value name)
    "clusterHarmonicMWpStr",
}


def is_pvalue_key(key: str) -> bool:
    if key in EXPLICIT_P_KEYS:
        return True
    # camelCase heuristic: starts with lowercase 'p' then uppercase
    if len(key) >= 2 and key[0] == 'p' and key[1].isupper():
        return True
    return False


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    store = ResultsStore()
    data = store.all()

    def _sanitize_macro_name(key: str) -> str:
        """Convert store key to a valid LaTeX macro name (no underscores).
        e.g. pCircA_validated -> pCircAValidated
        """
        parts = key.split("_")
        return parts[0] + "".join(p.capitalize() for p in parts[1:])

    # Build sig-macro companions
    sig_macros = {}
    for key, val in data.items():
        if not is_pvalue_key(key):
            continue
        try:
            p = float(val)
        except (TypeError, ValueError):
            continue
        if not (0.0 <= p <= 1.0):
            continue
        macro_name = _sanitize_macro_name(key) + "Sig"
        sig_macros[macro_name] = sig(p)

    print("=" * 72)
    print("  EMIT SIG MACROS — auto-generated significance labels")
    print(f"  Source: data/store/results.json  ({len(sig_macros)} sig macros)")
    print("=" * 72)
    print()
    for macro_name, label in sorted(sig_macros.items()):
        print(f"  \\newcommand{{\\{macro_name}}}{{{label}}}  % auto sig label")

    # Write back to store so downstream scripts can read them
    store.write_many(sig_macros)

    print()
    print(f"  ✓ {len(sig_macros)} significance-label macros emitted.")


if __name__ == "__main__":
    main()
