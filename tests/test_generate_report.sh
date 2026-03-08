#!/usr/bin/env bash
#
# test_generate_report.sh
#
# Tests for scripts/generate_report.py:
#   - GWAS threshold constants are defined
#   - Report generates with all expected interactive elements
#   - Summary statistics include new threshold-based counts
#   - Sample JSON data is correctly serialised for Plotly
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

PASS=0
FAIL=0

echo "============================================"
echo "  generate_report.py Tests"
echo "============================================"
echo ""

# --- Create mock pipeline output ---
mkdir -p "${TMP_DIR}/stage2/qc" "${TMP_DIR}/stage1/qc" "${TMP_DIR}/ancestry_pca" \
         "${TMP_DIR}/realigned_manifests" "${TMP_DIR}/stage2/qc/comparison" \
         "${TMP_DIR}/stage2/qc/variant_qc" "${TMP_DIR}/sex_check"

cat > "${TMP_DIR}/stage2/qc/stage2_sample_qc.tsv" << 'EOF'
sample_id	call_rate	lrr_sd	lrr_mean	lrr_median	baf_sd	het_rate	computed_gender
S001	0.9985	0.1523	-0.0012	-0.0005	0.0321	0.2345	M
S002	0.9654	0.4123	0.0321	0.0215	0.1623	0.2890	F
S003	0.9992	0.1234	-0.0003	-0.0001	0.0287	0.2356	M
EOF

cat > "${TMP_DIR}/stage1/qc/stage1_sample_qc.tsv" << 'EOF'
sample_id	call_rate	lrr_sd	lrr_mean	lrr_median	baf_sd	het_rate	computed_gender
S001	0.9945	0.1723	-0.0018	-0.0009	0.0345	0.2312	M
S002	0.9554	0.4523	0.0389	0.0278	0.1712	0.2934	F
S003	0.9972	0.1434	-0.0009	-0.0004	0.0298	0.2323	M
EOF

cat > "${TMP_DIR}/ancestry_pca/pca_projections.tsv" << 'EOF'
FID	IID	PC1	PC2	PC3	PC4
S001	S001	-0.0200	0.0100	0.0020	0.0010
S002	S002	0.0100	-0.0200	-0.0010	0.0005
S003	S003	0.0300	0.0200	0.0015	-0.0020
EOF

cat > "${TMP_DIR}/ancestry_pca/qc_summary.txt" << 'EOF'
Ancestry PCA QC summary (mock)
EOF

cat > "${TMP_DIR}/realigned_manifests/realign_summary.txt" << 'EOF'
Manifest realignment summary (mock)
EOF

cat > "${TMP_DIR}/stage2/qc/comparison/qc_comparison_report.txt" << 'EOF'
QC comparison report (mock)
EOF

cat > "${TMP_DIR}/stage2/qc/variant_qc/variant_qc_summary.txt" << 'EOF'
Variant QC summary (mock)
EOF

printf 'mock-png' > "${TMP_DIR}/stage2/qc/comparison/qc_comparison_samples.png"
printf 'mock-png' > "${TMP_DIR}/stage2/qc/comparison/qc_comparison_variants.png"
printf 'mock-png' > "${TMP_DIR}/sex_check/sex_check_chrXY_lrr.png"
cat > "${TMP_DIR}/sex_check/sex_check_chrXY_lrr.tsv" << 'EOF'
sample_id	chrx_lrr_median	chry_lrr_median	computed_gender
S001	-0.1200	0.0900	M
S002	0.0400	-0.0800	F
S003	-0.1000	0.0700	M
EOF

cat > "${TMP_DIR}/mock_data_notice.txt" << 'EOF'
This report is an example generated from deterministic mock data for demonstration purposes only.
EOF

# --- Test 1: Generate report ---
echo "--- Test 1: Report generation ---"
python3 "${REPO_DIR}/scripts/generate_report.py" --output-dir "${TMP_DIR}" >/dev/null 2>&1
REPORT="${TMP_DIR}/pipeline_report.html"

if [[ -f "${REPORT}" ]]; then
    echo "  PASS: Report file generated"
    (( PASS++ )) || true
else
    echo "  FAIL: Report file not generated"
    (( FAIL++ )) || true
fi

# --- Test 2: Report contains Plotly CDN reference ---
echo "--- Test 2: Plotly.js CDN reference ---"
if grep -q 'plotly' "${REPORT}"; then
    echo "  PASS: Plotly CDN reference found"
    (( PASS++ )) || true
else
    echo "  FAIL: Plotly CDN reference missing"
    (( FAIL++ )) || true
fi

# --- Test 3: Report contains interactive plot containers ---
echo "--- Test 3: Interactive plot containers ---"
for div_id in plot-scatter plot-cr plot-lrr plot-pca12 plot-pca34 plot-sex-check; do
    if grep -q "id=\"${div_id}\"" "${REPORT}"; then
        echo "  PASS: Plot container '${div_id}' found"
        (( PASS++ )) || true
    else
        echo "  FAIL: Plot container '${div_id}' missing"
        (( FAIL++ )) || true
    fi
done

# --- Test 4: Mock data banner ---
echo "--- Test 4: Mock data banner ---"
if grep -q 'Example Output (Mock Data)' "${REPORT}"; then
    echo "  PASS: Mock/example banner present"
    (( PASS++ )) || true
else
    echo "  FAIL: Mock/example banner missing"
    (( FAIL++ )) || true
fi

# --- Test 4b: Scatter/density controls ---
echo "--- Test 4b: Scatter/density controls ---"
if grep -q "label:'Scatter'" "${REPORT}" && grep -q "label:'Density'" "${REPORT}"; then
    echo "  PASS: Scatter/density mode toggle controls present"
    (( PASS++ )) || true
else
    echo "  FAIL: Scatter/density mode toggle controls missing"
    (( FAIL++ )) || true
fi

if grep -q "size:5" "${REPORT}"; then
    echo "  PASS: PCA default point size updated"
    (( PASS++ )) || true
else
    echo "  FAIL: PCA default point size not updated"
    (( FAIL++ )) || true
fi

# --- Test 5: GWAS QC thresholds reference table ---
echo "--- Test 5: GWAS QC thresholds section ---"
if grep -q 'GWAS QC Best Practice Thresholds' "${REPORT}"; then
    echo "  PASS: GWAS QC thresholds section present"
    (( PASS++ )) || true
else
    echo "  FAIL: GWAS QC thresholds section missing"
    (( FAIL++ )) || true
fi

if grep -q 'threshold-table' "${REPORT}"; then
    echo "  PASS: Threshold table CSS class present"
    (( PASS++ )) || true
else
    echo "  FAIL: Threshold table CSS class missing"
    (( FAIL++ )) || true
fi

if grep -q 'doi:10.1038/nprot.2010.116' "${REPORT}" && \
   grep -q 'doi:10.1002/mpr.1608' "${REPORT}" && \
   grep -q 'doi:10.1002/0471142905.hg0119s68' "${REPORT}"; then
    echo "  PASS: Threshold evidence citations present"
    (( PASS++ )) || true
else
    echo "  FAIL: Threshold evidence citations missing"
    (( FAIL++ )) || true
fi

if grep -q 'PCA step uses MAF ≥ 0.05 for stable loadings' "${REPORT}"; then
    echo "  PASS: MAF context note present"
    (( PASS++ )) || true
else
    echo "  FAIL: MAF context note missing"
    (( FAIL++ )) || true
fi

# --- Test 6: Per-sample QC table ---
echo "--- Test 6: Per-sample QC table ---"
if grep -q 'sample-tbody' "${REPORT}"; then
    echo "  PASS: Sample table body element found"
    (( PASS++ )) || true
else
    echo "  FAIL: Sample table body element missing"
    (( FAIL++ )) || true
fi

if grep -q 'sample-search' "${REPORT}"; then
    echo "  PASS: Sample search input found"
    (( PASS++ )) || true
else
    echo "  FAIL: Sample search input missing"
    (( FAIL++ )) || true
fi

# --- Test 7: QC data JSON embedded ---
echo "--- Test 7: QC data JSON ---"
if grep -q 'id="qc-data"' "${REPORT}"; then
    echo "  PASS: QC data JSON block found"
    (( PASS++ )) || true
else
    echo "  FAIL: QC data JSON block missing"
    (( FAIL++ )) || true
fi

# Verify JSON contains sample IDs
if grep -q '"S001"' "${REPORT}" && grep -q '"S002"' "${REPORT}"; then
    echo "  PASS: Sample IDs present in embedded JSON"
    (( PASS++ )) || true
else
    echo "  FAIL: Sample IDs missing from embedded JSON"
    (( FAIL++ )) || true
fi

# --- Test 7b: PCA data JSON embedded ---
echo "--- Test 7b: PCA data JSON ---"
if grep -q 'id="pca-data"' "${REPORT}" && grep -q '"pc1"' "${REPORT}" && grep -q '"S001"' "${REPORT}"; then
    echo "  PASS: PCA data JSON block found with sample data"
    (( PASS++ )) || true
else
    echo "  FAIL: PCA data JSON block missing or incomplete"
    (( FAIL++ )) || true
fi

# --- Test 7c: Sex check data JSON embedded ---
echo "--- Test 7c: Sex check data JSON ---"
if grep -q 'id="sex-data"' "${REPORT}" && grep -q '"chrx_lrr_median"' "${REPORT}" && grep -q '"S001"' "${REPORT}"; then
    echo "  PASS: Sex check data JSON block found with sample data"
    (( PASS++ )) || true
else
    echo "  FAIL: Sex check data JSON block missing or incomplete"
    (( FAIL++ )) || true
fi

# --- Test 8: Pass/fail classification in JSON ---
echo "--- Test 8: Pass/fail classification ---"
# S001 and S003 should pass (CR >= 0.97, LRR SD <= 0.35)
# S002 should fail (CR 0.9654 < 0.97, LRR SD 0.4123 > 0.35)
if grep -q '"qc_pass": true' "${REPORT}" && grep -q '"qc_pass": false' "${REPORT}"; then
    echo "  PASS: Pass/fail classification present in JSON"
    (( PASS++ )) || true
else
    echo "  FAIL: Pass/fail classification missing"
    (( FAIL++ )) || true
fi

# --- Test 9: Stage comparison table ---
echo "--- Test 9: Stage comparison ---"
if grep -q 'Stage 1' "${REPORT}" && grep -q 'Stage 2' "${REPORT}"; then
    echo "  PASS: Stage comparison present"
    (( PASS++ )) || true
else
    echo "  FAIL: Stage comparison missing"
    (( FAIL++ )) || true
fi

# --- Test 10: Summary statistics TSV ---
echo "--- Test 10: Summary statistics TSV ---"
STATS="${TMP_DIR}/summary_statistics.tsv"
if [[ -f "${STATS}" ]]; then
    echo "  PASS: Summary statistics TSV generated"
    (( PASS++ )) || true
else
    echo "  FAIL: Summary statistics TSV not generated"
    (( FAIL++ )) || true
fi

# --- Test 11: Methods text ---
echo "--- Test 11: Methods text ---"
METHODS="${TMP_DIR}/methods_text.txt"
if [[ -f "${METHODS}" ]]; then
    echo "  PASS: Methods text file generated"
    (( PASS++ )) || true
else
    echo "  FAIL: Methods text file not generated"
    (( FAIL++ )) || true
fi

# --- Test 11b: Compiled citations summary ---
echo "--- Test 11b: Compiled citations summary ---"
CITATIONS="${TMP_DIR}/citations_summary.tsv"
if [[ -f "${CITATIONS}" ]]; then
    echo "  PASS: Citations summary TSV generated"
    (( PASS++ )) || true
else
    echo "  FAIL: Citations summary TSV not generated"
    (( FAIL++ )) || true
fi

if grep -q 'GWAS Methods and Best-Practice Citations' "${REPORT}" && \
   grep -q 'doi:10.1093/gigascience/giab008' "${REPORT}"; then
    echo "  PASS: Citations section and DOI links present in report"
    (( PASS++ )) || true
else
    echo "  FAIL: Citations section or DOI links missing in report"
    (( FAIL++ )) || true
fi

# --- Test 12: Modern CSS variables ---
echo "--- Test 12: Modern CSS ---"
if grep -q '\-\-primary' "${REPORT}" && grep -q '\-\-success' "${REPORT}"; then
    echo "  PASS: CSS custom properties (variables) present"
    (( PASS++ )) || true
else
    echo "  FAIL: CSS custom properties missing"
    (( FAIL++ )) || true
fi

# --- Test 13: Consolidated summary directory ---
echo "--- Test 13: Consolidated summary directory ---"
SUMMARY_DIR="${TMP_DIR}/summary"
if [[ -d "${SUMMARY_DIR}" ]]; then
    echo "  PASS: Summary directory generated"
    (( PASS++ )) || true
else
    echo "  FAIL: Summary directory missing"
    (( FAIL++ )) || true
fi

for required_file in \
    pipeline_report.html \
    summary_statistics.tsv \
    methods_text.txt \
    citations_summary.tsv \
    pca_projections.tsv \
    ancestry_pca_qc_summary.txt \
    realign_summary.txt \
    qc_comparison_report.txt \
    variant_qc_summary_stage2.txt \
    sex_check_chrXY_lrr.png \
    sex_check_chrXY_lrr.tsv; do
    if [[ -f "${SUMMARY_DIR}/${required_file}" ]]; then
        echo "  PASS: Summary artifact '${required_file}' present"
        (( PASS++ )) || true
    else
        echo "  FAIL: Summary artifact '${required_file}' missing"
        (( FAIL++ )) || true
    fi
done

# --- Test 14: BAF flag count in JSON ---
echo "--- Test 14: BAF SD flag ---"
# S002 has BAF SD 0.1623 > 0.15, should be flagged
if grep -q '"baf_flag": true' "${REPORT}"; then
    echo "  PASS: BAF SD flag present in JSON"
    (( PASS++ )) || true
else
    echo "  FAIL: BAF SD flag missing"
    (( FAIL++ )) || true
fi

# --- Test 15: GWAS threshold constants ---
echo "--- Test 15: GWAS threshold constants in Python ---"
if python3 -c "
import sys
sys.path.insert(0, '${REPO_DIR}/scripts')
from generate_report import GWAS_THRESHOLDS
assert GWAS_THRESHOLDS['call_rate_min'] == 0.97
assert GWAS_THRESHOLDS['lrr_sd_max'] == 0.35
assert GWAS_THRESHOLDS['baf_sd_max'] == 0.15
assert GWAS_THRESHOLDS['het_rate_sd_multiplier'] == 3
assert GWAS_THRESHOLDS['hwe_pvalue_min'] == 1e-6
assert GWAS_THRESHOLDS['maf_min'] == 0.01
print('OK')
" 2>/dev/null; then
    echo "  PASS: GWAS threshold constants valid"
    (( PASS++ )) || true
else
    echo "  FAIL: GWAS threshold constants error"
    (( FAIL++ )) || true
fi

# --- Test 16: compute_summary_stats enhancements ---
echo "--- Test 16: compute_summary_stats enhancements ---"
if python3 -c "
import sys
sys.path.insert(0, '${REPO_DIR}/scripts')
from generate_report import compute_summary_stats, read_tsv
_, rows = read_tsv('${TMP_DIR}/stage2/qc/stage2_sample_qc.tsv')
stats = compute_summary_stats(rows)
assert 'n_baf_flag' in stats, 'n_baf_flag missing'
assert 'n_het_outlier' in stats, 'n_het_outlier missing'
assert 'pass_rate' in stats, 'pass_rate missing'
assert stats['n_pass'] == 2, f'n_pass expected 2, got {stats[\"n_pass\"]}'
assert stats['n_fail'] == 1, f'n_fail expected 1, got {stats[\"n_fail\"]}'
assert stats['n_baf_flag'] == 1, f'n_baf_flag expected 1, got {stats[\"n_baf_flag\"]}'
print('OK')
" 2>/dev/null; then
    echo "  PASS: compute_summary_stats has threshold-based counts"
    (( PASS++ )) || true
else
    echo "  FAIL: compute_summary_stats missing enhancements"
    (( FAIL++ )) || true
fi

# --- Test 17: Boundary value pass/fail classification ---
echo "--- Test 17: Boundary value pass/fail classification ---"
BOUNDARY_DIR="${TMP_DIR}/boundary"
mkdir -p "${BOUNDARY_DIR}"
cat > "${BOUNDARY_DIR}/boundary_qc.tsv" << 'BEOF'
sample_id	call_rate	lrr_sd	lrr_mean	lrr_median	baf_sd	het_rate	computed_gender
EXACT_PASS	0.9700	0.3500	0.0000	0.0000	0.0100	0.2500	M
JUST_BELOW_CR	0.9699	0.3500	0.0000	0.0000	0.0100	0.2500	F
JUST_ABOVE_LRR	0.9700	0.3501	0.0000	0.0000	0.0100	0.2500	M
BOTH_FAIL	0.9699	0.3501	0.0000	0.0000	0.0100	0.2500	F
BEOF

if python3 -c "
import sys, json
sys.path.insert(0, '${REPO_DIR}/scripts')
from generate_report import compute_summary_stats, read_tsv, _prepare_sample_json

_, rows = read_tsv('${BOUNDARY_DIR}/boundary_qc.tsv')
stats = compute_summary_stats(rows)
data = json.loads(_prepare_sample_json(rows, stats))

results = {s['id']: s['qc_pass'] for s in data}

# CR=0.97 and LRR_SD=0.35 should pass (>= and <=)
assert results['EXACT_PASS'] == True, f'EXACT_PASS should pass, got {results[\"EXACT_PASS\"]}'
# CR=0.9699 is below 0.97 threshold
assert results['JUST_BELOW_CR'] == False, f'JUST_BELOW_CR should fail, got {results[\"JUST_BELOW_CR\"]}'
# LRR_SD=0.3501 is above 0.35 threshold
assert results['JUST_ABOVE_LRR'] == False, f'JUST_ABOVE_LRR should fail, got {results[\"JUST_ABOVE_LRR\"]}'
# Both below threshold
assert results['BOTH_FAIL'] == False, f'BOTH_FAIL should fail, got {results[\"BOTH_FAIL\"]}'

assert stats['n_pass'] == 1, f'expected 1 pass, got {stats[\"n_pass\"]}'
assert stats['n_fail'] == 3, f'expected 3 fail, got {stats[\"n_fail\"]}'
print('OK')
" 2>/dev/null; then
    echo "  PASS: Boundary values correctly classified"
    (( PASS++ )) || true
else
    echo "  FAIL: Boundary value classification incorrect"
    (( FAIL++ )) || true
fi

echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

if [[ "${FAIL}" -gt 0 ]]; then
    exit 1
fi
