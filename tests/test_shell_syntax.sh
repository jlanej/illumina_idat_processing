#!/usr/bin/env bash
#
# test_shell_syntax.sh
#
# Validate shell script syntax with bash -n and shellcheck.
# Ensures all shell scripts parse correctly and follow best practices.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PASS=0
FAIL=0

echo "============================================"
echo "  Shell Script Syntax Checks"
echo "============================================"
echo ""

# Find all shell scripts
SCRIPTS=()
while IFS= read -r -d '' f; do
    SCRIPTS+=("${f}")
done < <(find "${REPO_DIR}/scripts" "${REPO_DIR}/tests" -name "*.sh" -print0 2>/dev/null | sort -z)

# ---------------------------------------------------------------
# Test 1: bash -n syntax check
# ---------------------------------------------------------------
echo "--- bash -n syntax validation ---"
for script in "${SCRIPTS[@]}"; do
    relpath="${script#"${REPO_DIR}/"}"
    if bash -n "${script}" 2>/dev/null; then
        echo "  PASS: ${relpath}"
        (( PASS++ )) || true
    else
        echo "  FAIL: ${relpath}"
        bash -n "${script}" 2>&1 | sed 's/^/        /'
        (( FAIL++ )) || true
    fi
done

echo ""

# ---------------------------------------------------------------
# Test 2: shellcheck (if available)
# ---------------------------------------------------------------
echo "--- shellcheck analysis ---"
if command -v shellcheck &>/dev/null; then
    for script in "${SCRIPTS[@]}"; do
        relpath="${script#"${REPO_DIR}/"}"
        # SC1091: Not following sourced files (dynamic paths)
        # SC2086: Double quote to prevent globbing (intentional in some cases)
        if shellcheck -e SC1091 -S warning "${script}" 2>/dev/null; then
            echo "  PASS: ${relpath}"
            (( PASS++ )) || true
        else
            echo "  WARN: ${relpath}"
            shellcheck -e SC1091 -S warning "${script}" 2>&1 | head -20 | sed 's/^/        /'
            # Shellcheck warnings are informational, not hard failures
        fi
    done
else
    echo "  SKIP: shellcheck not installed"
fi

echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

if [[ "${FAIL}" -gt 0 ]]; then
    exit 1
fi
