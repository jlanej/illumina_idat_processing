#!/usr/bin/env bash
#
# test_compute_sample_qc.sh
#
# Regression tests for scripts/compute_sample_qc.py:
#   - Correct LRR median for odd/even counts
#   - Correct call rate and heterozygosity calculations
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

PASS=0
FAIL=0

echo "============================================"
echo "  compute_sample_qc.py Tests"
echo "============================================"
echo ""

OUT="${TMP_DIR}/qc.tsv"

# Four variants, two samples; all calls present to keep expectations explicit.
cat <<'EOF' | python3 "${REPO_DIR}/scripts/compute_sample_qc.py" --num-samples 2 --output "${OUT}" >/dev/null
	0/1:1.0:0.45	0/0:2.0:0.10
	0/0:3.0:0.20	0/1:4.0:0.60
	0/1:5.0:0.55	1/1:6.0:0.90
	1/1:7.0:0.90	0/1:8.0:0.40
EOF

check_field() {
    local idx="$1" col="$2" expected="$3" label="$4"
    local actual
    actual=$(awk -F'\t' -v i="${idx}" -v c="${col}" 'NR>1 && $1==i {print $c}' "${OUT}")
    if [[ "${actual}" == "${expected}" ]]; then
        echo "  PASS: ${label} = ${expected}"
        (( PASS++ )) || true
    else
        echo "  FAIL: ${label} expected=${expected} actual=${actual}"
        (( FAIL++ )) || true
    fi
}

# idx=0: LRR values 1,3,5,7 -> median 4.000000; 2 het calls out of 4
check_field 0 2 "1.000000" "sample0 call_rate"
check_field 0 5 "4.000000" "sample0 lrr_median"
check_field 0 7 "0.500000" "sample0 het_rate"

# idx=1: LRR values 2,4,6,8 -> median 5.000000; 2 het calls out of 4
check_field 1 2 "1.000000" "sample1 call_rate"
check_field 1 5 "5.000000" "sample1 lrr_median"
check_field 1 7 "0.500000" "sample1 het_rate"

echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

if [[ "${FAIL}" -gt 0 ]]; then
    exit 1
fi
