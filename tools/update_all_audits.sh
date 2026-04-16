#!/usr/bin/env bash
# update_all_audits.sh
# ====================
# Regenerate every supplementary audit file.
#
# Run from the repo root:
#   bash tools/update_all_audits.sh
#
# The FDR audit (fdr.txt) reads from data/store/results.json, which is
# populated by the main analysis scripts.  If you want the FDR audit to
# reflect the latest p-values, run  bash manuscript/reproduce_all_macros.sh
# first, or pass --full to this script:
#
#   bash tools/update_all_audits.sh --full
#
# With --full, this script runs reproduce_all_macros.sh before regenerating
# the audits so that results.json is up to date.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

FULL=0
for arg in "$@"; do
    [[ "$arg" == "--full" ]] && FULL=1
done

echo "=========================================="
echo "  Gerizim Paper-A — Update All Audits"
echo "  Repo: $REPO_ROOT"
echo "=========================================="

if [[ $FULL -eq 1 ]]; then
    echo ""
    echo "  --full: running reproduce_all_macros.sh first..."
    bash manuscript/reproduce_all_macros.sh
    echo ""
fi

echo ""
echo "  [1/6] Dome keyword audit..."
python3 tools/generate_audit_dome.py

echo "  [2/6] Dome + mound evolution audit..."
python3 tools/generate_audit_dome_mound.py

echo "  [3/6] Founding / sacred-origin audit..."
python3 tools/generate_audit_founding.py

echo "  [4/6] World religion keyword audit..."
python3 tools/generate_audit_religion.py

echo "  [5/6] Corridor precision audit..."
python3 tools/generate_audit_corridor.py

echo "  [6/7] Gerizim vs Jerusalem A-tier site comparison..."
python3 tools/generate_audit_aplus_sites.py

echo "  [7/7] FDR / multiple-comparisons audit..."
python3 tools/generate_audit_fdr.py

echo ""
echo "  Done.  Files written to supplementary/audit/"
echo ""
ls -lh "$REPO_ROOT/supplementary/audit/"*.txt
