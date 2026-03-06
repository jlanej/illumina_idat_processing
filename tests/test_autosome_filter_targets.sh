#!/usr/bin/env bash
#
# test_autosome_filter_targets.sh
#
# Ensure sample QC filtering explicitly targets autosomes and does not rely on
# "exclude sex chromosomes" negation, which can accidentally include alt contigs.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PASS=0
FAIL=0

echo "============================================"
echo "  Autosome Target Filter Tests"
echo "============================================"
echo ""

COLLECT_QC="${REPO_DIR}/scripts/collect_qc_metrics.sh"
COMPUTE_QC="${REPO_DIR}/scripts/compute_sample_qc.py"

if grep -Eq -- "-t \^chrX,chrY,chrM,X,Y,MT" "${COLLECT_QC}" "${COMPUTE_QC}"; then
    echo "  FAIL: legacy negation-based chromosome filter is still present"
    (( FAIL++ )) || true
else
    echo "  PASS: legacy negation-based chromosome filter removed"
    (( PASS++ )) || true
fi

if grep -Eq -- "autosome_targets=.*printf 'chr%s,' \{1\.\.22\}" "${COLLECT_QC}" && \
   grep -Eq -- "printf '%s,' \{1\.\.22\}" "${COLLECT_QC}" && \
   grep -Eq -- "-t \"\\$\\{autosome_targets\\}\"" "${COLLECT_QC}"; then
    echo "  PASS: collect_qc_metrics.sh explicitly targets chromosomes 1-22"
    (( PASS++ )) || true
else
    echo "  FAIL: explicit chromosome 1-22 targeting not found"
    (( FAIL++ )) || true
fi

echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

if [[ "${FAIL}" -gt 0 ]]; then
    exit 1
fi
