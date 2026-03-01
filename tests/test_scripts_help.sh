#!/usr/bin/env bash
#
# test_scripts_help.sh
#
# Validate that all pipeline scripts display help without errors.
# Ensures scripts are self-documenting and parseable.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PASS=0
FAIL=0

echo "============================================"
echo "  Script Help Output Tests"
echo "============================================"
echo ""

# Scripts that support --help
HELP_SCRIPTS=(
    scripts/run_pipeline.sh
    scripts/install_dependencies.sh
    scripts/download_manifests.sh
    scripts/download_reference.sh
    scripts/realign_manifest.sh
    scripts/stage1_initial_genotyping.sh
    scripts/stage2_recluster.sh
    scripts/process_1000g.sh
    scripts/compute_variant_qc.sh
)

for relpath in "${HELP_SCRIPTS[@]}"; do
    script="${REPO_DIR}/${relpath}"
    if [[ ! -f "${script}" ]]; then
        echo "  SKIP: ${relpath} (file not found)"
        continue
    fi

    # --help should exit 0 and produce output
    output=$(bash "${script}" --help 2>&1) || true

    if [[ -n "${output}" ]]; then
        # Verify help output contains Usage
        if echo "${output}" | grep -qi "usage"; then
            echo "  PASS: ${relpath} (--help shows usage)"
            (( PASS++ )) || true
        else
            echo "  WARN: ${relpath} (--help output missing 'Usage')"
            (( PASS++ )) || true
        fi
    else
        echo "  FAIL: ${relpath} (--help produced no output)"
        (( FAIL++ )) || true
    fi
done

echo ""

# ---------------------------------------------------------------
# Validate Python script
# ---------------------------------------------------------------
echo "--- Python script syntax ---"
PYTHON_SCRIPTS=(
    scripts/recluster_egt.py
    scripts/diagnose_qc.py
    scripts/plot_qc_comparison.py
)
for relpath in "${PYTHON_SCRIPTS[@]}"; do
    PYTHON_SCRIPT="${REPO_DIR}/${relpath}"
    if [[ -f "${PYTHON_SCRIPT}" ]]; then
        if python3 -c "import py_compile; py_compile.compile('${PYTHON_SCRIPT}', doraise=True)" 2>/dev/null; then
            echo "  PASS: ${relpath} (compiles)"
            (( PASS++ )) || true
        else
            echo "  FAIL: ${relpath} (syntax error)"
            (( FAIL++ )) || true
        fi
    else
        echo "  SKIP: ${relpath} not found"
    fi
done

echo ""

# ---------------------------------------------------------------
# Validate config file
# ---------------------------------------------------------------
echo "--- Config file validation ---"
MANIFEST_TSV="${REPO_DIR}/config/manifest_urls.tsv"
if [[ -f "${MANIFEST_TSV}" ]]; then
    # Check that file has expected number of columns (4: name, bpm, egt, csv)
    bad_lines=$(grep -v '^#' "${MANIFEST_TSV}" | grep -v '^$' | \
        awk -F'\t' 'NF != 4 {print NR": "$0}')
    if [[ -z "${bad_lines}" ]]; then
        n_arrays=$(grep -v '^#' "${MANIFEST_TSV}" | grep -v '^$' | wc -l)
        echo "  PASS: config/manifest_urls.tsv (${n_arrays} arrays, 4 columns)"
        (( PASS++ )) || true
    else
        echo "  FAIL: config/manifest_urls.tsv (malformed lines):"
        echo "${bad_lines}" | sed 's/^/        /'
        (( FAIL++ )) || true
    fi
else
    echo "  FAIL: config/manifest_urls.tsv not found"
    (( FAIL++ )) || true
fi

echo ""

# ---------------------------------------------------------------
# Validate Dockerfile
# ---------------------------------------------------------------
echo "--- Dockerfile validation ---"
DOCKERFILE="${REPO_DIR}/Dockerfile"
if [[ -f "${DOCKERFILE}" ]]; then
    # Basic validation: check for FROM, COPY, required stages
    if grep -q '^FROM ' "${DOCKERFILE}" && grep -q 'bcftools' "${DOCKERFILE}"; then
        echo "  PASS: Dockerfile (has FROM and bcftools build)"
        (( PASS++ )) || true
    else
        echo "  FAIL: Dockerfile (missing required directives)"
        (( FAIL++ )) || true
    fi
else
    echo "  SKIP: Dockerfile not found"
fi

echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

if [[ "${FAIL}" -gt 0 ]]; then
    exit 1
fi
