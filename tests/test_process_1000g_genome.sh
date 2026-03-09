#!/usr/bin/env bash
#
# test_process_1000g_genome.sh
#
# Validate genome argument handling in process_1000g.sh.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PASS=0
FAIL=0

echo "============================================"
echo "  process_1000g Genome Argument Tests"
echo "============================================"
echo ""

PROCESS_1000G="${REPO_DIR}/scripts/process_1000g.sh"

echo "--- Test 1: --help mentions --genome ---"
help_out=$(bash "${PROCESS_1000G}" --help 2>&1) || true
if echo "${help_out}" | grep -q -- "--genome"; then
    echo "  PASS: process_1000g --help mentions --genome"
    (( PASS++ )) || true
else
    echo "  FAIL: process_1000g --help does not mention --genome"
    (( FAIL++ )) || true
fi

echo "--- Test 2: process_1000g passes genome to run_pipeline ---"
if grep -q -- '--genome "\${GENOME}"' "${PROCESS_1000G}" && \
   grep -q -- '--genome)        GENOME="\$2"; shift 2 ;;' "${PROCESS_1000G}"; then
    echo "  PASS: process_1000g parses and passes --genome to run_pipeline"
    (( PASS++ )) || true
else
    echo "  FAIL: process_1000g missing --genome parse/pass-through wiring"
    (( FAIL++ )) || true
fi

echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

if [[ "${FAIL}" -gt 0 ]]; then
    exit 1
fi
