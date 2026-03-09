#!/usr/bin/env bash
#
# test_temp_dir_scoping.sh
#
# Ensure temporary files/directories are created under output-dir scoped tmp paths.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PASS=0
FAIL=0

echo "============================================"
echo "  Temp Directory Scoping Tests"
echo "============================================"
echo ""

STAGE1="${REPO_DIR}/scripts/stage1_initial_genotyping.sh"
STAGE2="${REPO_DIR}/scripts/stage2_recluster.sh"
COLLECT_QC="${REPO_DIR}/scripts/collect_qc_metrics.sh"
VARIANT_QC="${REPO_DIR}/scripts/compute_variant_qc.sh"
ANCESTRY_PCA="${REPO_DIR}/scripts/ancestry_pca.sh"
RUN_PEDDY="${REPO_DIR}/scripts/run_peddy.sh"

if grep -Eq 'TMP_DIR="\$\{OUTPUT_DIR\}/tmp/stage1_initial_genotyping"' "${STAGE1}" && \
   grep -Eq 'mktemp -p "\$\{TMP_DIR\}" stage1_batch\.' "${STAGE1}" && \
   grep -Eq 'mktemp -p "\$\{TMP_DIR\}" stage1_sample\.' "${STAGE1}" && \
   grep -Eq 'bcftools sort .* -T "\$\{TMP_DIR\}/bcftools\."' "${STAGE1}"; then
    echo "  PASS: stage1 temp files are scoped under \${OUTPUT_DIR}/tmp"
    (( PASS++ )) || true
else
    echo "  FAIL: stage1 temp file scoping is missing or incomplete"
    (( FAIL++ )) || true
fi

if grep -Eq 'TMP_DIR="\$\{OUTPUT_DIR\}/tmp/stage2_recluster"' "${STAGE2}" && \
   grep -Eq 'mktemp -p "\$\{TMP_DIR\}" stage2_batch\.' "${STAGE2}" && \
   grep -Eq 'mktemp -p "\$\{TMP_DIR\}" stage2_sample\.' "${STAGE2}" && \
   grep -Eq 'bcftools sort .* -T "\$\{TMP_DIR\}/bcftools\."' "${STAGE2}"; then
    echo "  PASS: stage2 temp files are scoped under \${OUTPUT_DIR}/tmp"
    (( PASS++ )) || true
else
    echo "  FAIL: stage2 temp file scoping is missing or incomplete"
    (( FAIL++ )) || true
fi

if grep -Eq 'tmp_root="\$\(dirname "\$\{output_file\}"\)/tmp/collect_qc_metrics"' "${COLLECT_QC}" && \
   grep -Eq 'mktemp -d -p "\$\{tmp_root\}"' "${COLLECT_QC}" && \
   grep -Eq 'mktemp -p "\$\{tmp_root\}"' "${COLLECT_QC}"; then
    echo "  PASS: collect_qc_metrics temp artifacts are scoped under output file tmp"
    (( PASS++ )) || true
else
    echo "  FAIL: collect_qc_metrics temp artifact scoping is missing or incomplete"
    (( FAIL++ )) || true
fi

if grep -Eq 'TMP_DIR="\$\{OUTPUT_DIR\}/tmp/compute_variant_qc"' "${VARIANT_QC}" && \
   grep -Eq 'mktemp -p "\$\{TMP_DIR\}".*filtered_variant_qc' "${VARIANT_QC}"; then
    echo "  PASS: compute_variant_qc temp files are scoped under \${OUTPUT_DIR}/tmp"
    (( PASS++ )) || true
else
    echo "  FAIL: compute_variant_qc temp file scoping is missing or incomplete"
    (( FAIL++ )) || true
fi

if grep -Eq 'TMP_DIR="\$\{OUTPUT_DIR\}/tmp/ancestry_pca"' "${ANCESTRY_PCA}" && \
   grep -Eq 'mktemp -p "\$\{TMP_DIR\}".*filtered_ancestry_pca' "${ANCESTRY_PCA}"; then
    echo "  PASS: ancestry_pca temp files are scoped under \${OUTPUT_DIR}/tmp"
    (( PASS++ )) || true
else
    echo "  FAIL: ancestry_pca temp file scoping is missing or incomplete"
    (( FAIL++ )) || true
fi

if grep -Eq 'TMP_DIR="\$\{OUTPUT_DIR\}/tmp/run_peddy"' "${RUN_PEDDY}" && \
   grep -Eq 'bcftools sort -T "\$\{TMP_DIR\}/bcftools\."' "${RUN_PEDDY}"; then
    echo "  PASS: run_peddy temp directory is scoped under \${OUTPUT_DIR}/tmp"
    (( PASS++ )) || true
else
    echo "  FAIL: run_peddy temp directory scoping is missing or incomplete"
    (( FAIL++ )) || true
fi

echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

if [[ "${FAIL}" -gt 0 ]]; then
    exit 1
fi
