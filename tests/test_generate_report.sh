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

# --- Create mock peddy output ---
mkdir -p "${TMP_DIR}/peddy"

cat > "${TMP_DIR}/peddy/peddy.het_check.csv" << 'EOF'
sample_id,sampled_sites,mean_depth,median_depth,depth_outlier,het_count,het_ratio,ratio_outlier,idr_baf,p10,p90,PC1,PC2,PC3,PC4,ancestry-prediction,ancestry-prob
S001,18000,32.5,32.5,False,4200,0.2333,False,0.02,0.01,0.88,-0.015,0.008,0.001,0.002,EUR,0.95
S002,17500,28.3,28.3,False,5100,0.2914,False,0.03,0.02,0.90,0.008,-0.015,-0.001,0.001,AFR,0.88
S003,19000,35.1,35.1,False,4400,0.2316,False,0.02,0.01,0.89,0.025,0.018,0.002,-0.003,EAS,0.92
EOF

cat > "${TMP_DIR}/peddy/peddy.sex_check.csv" << 'EOF'
sample_id,error,het_count,hom_alt_count,hom_ref_count,het_ratio,ped_sex,predicted_sex
S001,False,15,300,5000,0.0028,male,male
S002,False,2200,1500,3000,0.3284,female,female
S003,False,10,350,5200,0.0018,male,male
EOF

cat > "${TMP_DIR}/peddy/peddy.ped_check.csv" << 'EOF'
sample_a,sample_b,n,rel,pedigree_relatedness,rel_difference,ibs0,ibs2,shared_hets,hets_a,hets_b,pedigree_parents
S001,S002,18000,0.0123,0.0,0.0123,2000,10000,3000,4000,4500,False
S001,S003,18000,-0.0034,0.0,-0.0034,2500,9500,2800,4000,4200,False
EOF

cat > "${TMP_DIR}/peddy/peddy.peddy.ped" << 'EOF'
S001	S001	0	0	1	-9
S002	S002	0	0	2	-9
S003	S003	0	0	1	-9
EOF

cat > "${TMP_DIR}/peddy/peddy_final.ped" << 'EOF'
S001	S001	0	0	1	-9
S002	S002	0	0	2	-9
S003	S003	0	0	1	-9
EOF

# --- Create mock ancestry-stratified QC output ---
mkdir -p "${TMP_DIR}/ancestry_stratified_qc/EUR/variant_qc" \
         "${TMP_DIR}/ancestry_stratified_qc/AFR/variant_qc"

cat > "${TMP_DIR}/ancestry_stratified_qc/EUR/variant_qc/variant_qc_summary.txt" << 'EOF'
Variant QC summary (EUR subset, mock)
EOF

cat > "${TMP_DIR}/ancestry_stratified_qc/AFR/variant_qc/variant_qc_summary.txt" << 'EOF'
Variant QC summary (AFR subset, mock)
EOF

cat > "${TMP_DIR}/ancestry_stratified_qc/ancestry_stratified_summary.txt" << 'EOF'
Ancestry-Stratified QC Summary (mock)
EUR: 1 sample
AFR: 1 sample
EOF

cat > "${TMP_DIR}/ancestry_stratified_qc/collated_variant_qc.tsv" << 'EOF'
variant_id	all_call_rate	all_hwe_p	all_maf	AFR_call_rate	AFR_hwe_p	AFR_maf	EUR_call_rate	EUR_hwe_p	EUR_maf	all_ancestries_call_rate_pass	all_ancestries_hwe_pass	all_ancestries_maf_pass	all_ancestries_qc_pass
rs100000	0.985	1.2e-03	0.15	0.98	2.5e-04	0.12	0.99	5.1e-02	0.18	1	1	1	1
rs100001	0.992	4.5e-01	0.32	0.99	3.2e-01	0.28	0.99	6.7e-01	0.35	1	1	1	1
rs100002	0.978	8.9e-08	0.05	0.96	1.1e-06	0.08	0.99	2.3e-01	0.03	0	1	1	0
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
for div_id in plot-scatter plot-cr plot-lrr plot-pca12 plot-pca34 plot-pca3d plot-sex-check; do
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
if grep -q 'id="pca-scatter-btn"' "${REPORT}" && grep -q 'id="pca-density-btn"' "${REPORT}" && \
   grep -q 'id="sex-scatter-btn"' "${REPORT}" && grep -q 'id="sex-density-btn"' "${REPORT}"; then
    echo "  PASS: Scatter/density external button controls present for both PCA and sex-check"
    (( PASS++ )) || true
else
    echo "  FAIL: Scatter/density external button controls missing"
    (( FAIL++ )) || true
fi

if grep -q 'id="pca-density-bins"' "${REPORT}" && \
   grep -q 'id="sex-density-bins"' "${REPORT}"; then
    echo "  PASS: Density bin slider controls present"
    (( PASS++ )) || true
else
    echo "  FAIL: Density bin slider controls missing"
    (( FAIL++ )) || true
fi

if grep -q "size:5" "${REPORT}"; then
    echo "  PASS: PCA default point size updated"
    (( PASS++ )) || true
else
    echo "  FAIL: PCA default point size not updated"
    (( FAIL++ )) || true
fi

if grep -q 'id="pca-view-toggle"' "${REPORT}" && \
   grep -q 'id="pca-cluster-run"' "${REPORT}" && \
   grep -q 'id="pca-cluster-export"' "${REPORT}" && \
   grep -q 'ancestry_cluster_assignments.tsv' "${REPORT}"; then
    echo "  PASS: PCA 3D toggle and clustering/export controls present"
    (( PASS++ )) || true
else
    echo "  FAIL: PCA 3D toggle and clustering/export controls missing"
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

# --- Test 18: Peddy data JSON ---
echo "--- Test 18: Peddy data JSON ---"
if grep -q 'id="peddy-data"' "${TMP_DIR}/pipeline_report.html" && \
   grep -q '"ancestry"' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Peddy data JSON block found with ancestry data"
    (( PASS++ )) || true
else
    echo "  FAIL: Peddy data JSON block missing"
    (( FAIL++ )) || true
fi

# --- Test 18b: Peddy section ---
echo "--- Test 18b: Peddy QC section ---"
if grep -q 'id="peddy"' "${TMP_DIR}/pipeline_report.html" && \
   grep -q 'Peddy QC' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Peddy QC section present in report"
    (( PASS++ )) || true
else
    echo "  FAIL: Peddy QC section missing from report"
    (( FAIL++ )) || true
fi

# --- Test 18c: Peddy PCA plot containers ---
echo "--- Test 18c: Peddy PCA plot containers ---"
if grep -q 'plot-peddy-pca12' "${TMP_DIR}/pipeline_report.html" && \
   grep -q 'plot-peddy-pca34' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Peddy PCA plot containers found"
    (( PASS++ )) || true
else
    echo "  FAIL: Peddy PCA plot containers missing"
    (( FAIL++ )) || true
fi

# --- Test 18d: Peddy ancestry color toggle button ---
echo "--- Test 18d: Peddy ancestry color toggle ---"
if grep -q 'peddy-overlay-toggle' "${TMP_DIR}/pipeline_report.html" && \
   grep -q 'Color by peddy ancestry' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Peddy ancestry color toggle present"
    (( PASS++ )) || true
else
    echo "  FAIL: Peddy ancestry color toggle missing"
    (( FAIL++ )) || true
fi

# --- Test 18e: Peddy nav link ---
echo "--- Test 18e: Peddy nav link ---"
if grep -q 'href="#peddy"' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Peddy nav link present"
    (( PASS++ )) || true
else
    echo "  FAIL: Peddy nav link missing"
    (( FAIL++ )) || true
fi

# --- Test 18ea: Peddy before Ancestry PCA in nav/section order ---
echo "--- Test 18ea: Peddy before Ancestry PCA order ---"
peddy_nav_line="$(grep -n 'href="#peddy"' "${TMP_DIR}/pipeline_report.html" | head -n1 | cut -d: -f1 || true)"
pca_nav_line="$(grep -n 'href="#pca"' "${TMP_DIR}/pipeline_report.html" | head -n1 | cut -d: -f1 || true)"
peddy_section_line="$(grep -n '<section id="peddy">' "${TMP_DIR}/pipeline_report.html" | head -n1 | cut -d: -f1 || true)"
pca_section_line="$(grep -n '<section id="pca">' "${TMP_DIR}/pipeline_report.html" | head -n1 | cut -d: -f1 || true)"
if [[ -n "${peddy_nav_line}" && -n "${pca_nav_line}" && -n "${peddy_section_line}" && -n "${pca_section_line}" && \
      "${peddy_nav_line}" -lt "${pca_nav_line}" && "${peddy_section_line}" -lt "${pca_section_line}" ]]; then
    echo "  PASS: Report orders peddy before ancestry PCA in nav and sections"
    (( PASS++ )) || true
else
    echo "  FAIL: Report order should be peddy before ancestry PCA"
    (( FAIL++ )) || true
fi

# --- Test 18f: Peddy citation ---
echo "--- Test 18f: Peddy citation ---"
if grep -q 'Pedersen and Quinlan' "${TMP_DIR}/pipeline_report.html" && \
   grep -q '10.1016/j.ajhg.2017.01.017' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Peddy citation (Pedersen and Quinlan 2017) present"
    (( PASS++ )) || true
else
    echo "  FAIL: Peddy citation missing"
    (( FAIL++ )) || true
fi

# --- Test 18g: Peddy summary artifacts ---
echo "--- Test 18g: Peddy summary artifacts ---"
for artifact in peddy.het_check.csv peddy.sex_check.csv peddy.ped_check.csv peddy.peddy.ped peddy_final.ped; do
    if [[ -f "${TMP_DIR}/summary/${artifact}" ]]; then
        echo "  PASS: Summary artifact '${artifact}' present"
        (( PASS++ )) || true
    else
        echo "  FAIL: Summary artifact '${artifact}' missing"
        (( FAIL++ )) || true
    fi
done

# --- Test 18h: Peddy PCs section toggle ---
echo "--- Test 18h: Peddy PCs section toggle ---"
if grep -q 'peddy-pcs-toggle' "${TMP_DIR}/pipeline_report.html" && \
   grep -q 'Show peddy PCs' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Peddy PCs toggle button present"
    (( PASS++ )) || true
else
    echo "  FAIL: Peddy PCs toggle button missing"
    (( FAIL++ )) || true
fi

# --- Test 19: Ancestry-stratified variant QC toggles ---
echo "--- Test 19: Ancestry-stratified variant QC toggles ---"
if grep -q 'vqc-toggle-bar' "${TMP_DIR}/pipeline_report.html" && \
   grep -q 'class="vqc-toggle"' "${TMP_DIR}/pipeline_report.html" && \
   grep -q 'value="all"' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Variant QC toggle bar with 'All' option present"
    (( PASS++ )) || true
else
    echo "  FAIL: Variant QC toggle bar missing"
    (( FAIL++ )) || true
fi

# --- Test 19b: Per-ancestry variant QC toggles ---
echo "--- Test 19b: Per-ancestry variant QC toggles ---"
for anc in EUR AFR; do
    if grep -q "value=\"${anc}\"" "${TMP_DIR}/pipeline_report.html"; then
        echo "  PASS: Variant QC toggle for ${anc} present"
        (( PASS++ )) || true
    else
        echo "  FAIL: Variant QC toggle for ${anc} missing"
        (( FAIL++ )) || true
    fi
done

# --- Test 19c: Collated variant QC JSON data ---
echo "--- Test 19c: Collated variant QC JSON data ---"
if grep -q 'collated-vqc-data' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Collated variant QC JSON block present"
    (( PASS++ )) || true
else
    echo "  FAIL: Collated variant QC JSON block missing"
    (( FAIL++ )) || true
fi

# --- Test 19d: Variant QC multi-ancestry plot containers ---
echo "--- Test 19d: Variant QC multi-ancestry plot containers ---"
for plot_id in plot-vqc-cr plot-vqc-maf plot-vqc-hwe; do
    if grep -q "${plot_id}" "${TMP_DIR}/pipeline_report.html"; then
        echo "  PASS: Variant QC plot container ${plot_id} present"
        (( PASS++ )) || true
    else
        echo "  FAIL: Variant QC plot container ${plot_id} missing"
        (( FAIL++ )) || true
    fi
done

# --- Test 19e: Ancestry-stratified summary in report ---
echo "--- Test 19e: Ancestry-stratified summary ---"
if grep -q 'Ancestry-Stratified QC Summary' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Ancestry-stratified summary present in report"
    (( PASS++ )) || true
else
    echo "  FAIL: Ancestry-stratified summary missing from report"
    (( FAIL++ )) || true
fi

# --- Test 19f: Peterson et al. 2019 citation ---
echo "--- Test 19f: Peterson citation ---"
if grep -q 'Peterson' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Peterson et al. 2019 citation present"
    (( PASS++ )) || true
else
    echo "  FAIL: Peterson et al. 2019 citation missing"
    (( FAIL++ )) || true
fi

# --- Test 19g: Methods text mentions ancestry-stratified QC ---
echo "--- Test 19g: Methods text ancestry-stratified QC ---"
if grep -q 'ancestry-stratified' "${TMP_DIR}/methods_text.txt"; then
    echo "  PASS: Methods text mentions ancestry-stratified QC"
    (( PASS++ )) || true
else
    echo "  FAIL: Methods text does not mention ancestry-stratified QC"
    (( FAIL++ )) || true
fi

# --- Test 19h: Summary bundle includes collated variant QC ---
echo "--- Test 19h: Summary artifacts for ancestry-stratified QC ---"
for art in collated_variant_qc.tsv ancestry_stratified_summary.txt; do
    if [[ -f "${TMP_DIR}/summary/${art}" ]]; then
        echo "  PASS: Summary artifact '${art}' present"
        (( PASS++ )) || true
    else
        echo "  FAIL: Summary artifact '${art}' missing"
        (( FAIL++ )) || true
    fi
done

# --- Test 19i: Collated VQC uses full variant set for histograms ---
echo "--- Test 19i: Complete collated VQC histogram values ---"
if REPO_DIR="${REPO_DIR}" TMP_DIR="${TMP_DIR}" python3 - <<'PY'
import json
import os
import sys

repo_dir = os.environ['REPO_DIR']
tmp_dir = os.environ['TMP_DIR']
sys.path.insert(0, os.path.join(repo_dir, 'scripts'))
from generate_report import _prepare_collated_vqc_json  # noqa: E402

test_file = os.path.join(tmp_dir, 'collated_bias_test.tsv')
with open(test_file, 'w') as f:
    f.write('variant_id\tall_call_rate\tall_hwe_p\tall_maf\n')
    # Bimodal distribution verifies the payload retains every row (no subsampling).
    for i in range(10000):
        maf = 0.0 if i < 5000 else 0.5
        f.write(f'rs{i:06d}\t0.99\t0.5\t{maf}\n')

payload = json.loads(_prepare_collated_vqc_json(test_file))
hist = payload['all']['maf_hist']
total = sum(hist['counts'])
assert total == 10000, total
# Bin at edge 0.0 should hold 5000 (maf=0.0), bin at edge 0.5 should hold 5000 (maf=0.5)
assert hist['counts'][0] == 5000, hist['counts'][0]
assert hist['counts'][-1] == 5000, hist['counts'][-1]
PY
then
    echo "  PASS: Histogram payload includes all variants without subsampling"
    (( PASS++ )) || true
else
    echo "  FAIL: Histogram payload appears subsampled"
    (( FAIL++ )) || true
fi

# --- Test 19j: Cross-ancestry pass summary rendered ---
echo "--- Test 19j: Cross-ancestry pass summary in report ---"
if grep -q 'Cross-Ancestry Variant QC Pass Summary' "${TMP_DIR}/pipeline_report.html" && \
   grep -q 'vqc-cross-ancestry-summary' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Cross-ancestry pass summary card present"
    (( PASS++ )) || true
else
    echo "  FAIL: Cross-ancestry pass summary card missing"
    (( FAIL++ )) || true
fi

# --- Test 19ja: Cross-ancestry missing-data explanation text ---
echo "--- Test 19ja: Cross-ancestry missing-data explanations present ---"
if grep -q 'cross-ancestry pass flags were not found in collated variant QC' "${TMP_DIR}/pipeline_report.html" && \
   grep -q 'No values available for the selected ancestry group' "${TMP_DIR}/pipeline_report.html" && \
   grep -q 'need at least two groups in ancestry_stratified_qc/collated_variant_qc.tsv' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Report includes explicit missing-data explanations for variant QC plots"
    (( PASS++ )) || true
else
    echo "  FAIL: Missing-data explanations for variant QC plots not found"
    (( FAIL++ )) || true
fi

# --- Test 19jb: Collated VQC JSON excludes NaN/Inf tokens ---
echo "--- Test 19jb: Collated VQC JSON sanitizes NaN/Inf values ---"
if REPO_DIR="${REPO_DIR}" TMP_DIR="${TMP_DIR}" python3 - <<'PY'
import json
import os
import sys

repo_dir = os.environ['REPO_DIR']
tmp_dir = os.environ['TMP_DIR']
sys.path.insert(0, os.path.join(repo_dir, 'scripts'))
from generate_report import _prepare_collated_vqc_json  # noqa: E402

test_file = os.path.join(tmp_dir, 'collated_nan_test.tsv')
with open(test_file, 'w') as f:
    f.write('variant_id\tall_call_rate\tall_hwe_p\tall_maf\tEUR_call_rate\tEUR_hwe_p\tEUR_maf\n')
    f.write('rs1\tnan\t0.5\t0.1\t0.99\t0.4\t0.1\n')
    f.write('rs2\t0.98\tinf\t0.2\t0.97\t-inf\t0.2\n')
    f.write('rs3\t0.97\t0.3\tinf\t0.96\t0.3\t-inf\n')

payload_text = _prepare_collated_vqc_json(test_file)
assert 'NaN' not in payload_text and 'Infinity' not in payload_text and '-Infinity' not in payload_text, payload_text
payload = json.loads(payload_text)
assert payload['all']['n_cr'] == 2, payload
assert payload['all']['n_hwe'] == 2, payload
assert payload['all']['n_maf'] == 2, payload
assert payload['EUR']['n_cr'] == 3, payload
assert payload['EUR']['n_hwe'] == 2, payload
assert payload['EUR']['n_maf'] == 2, payload
PY
then
    echo "  PASS: Collated VQC JSON payload is valid JSON when TSV contains NaN/Inf tokens"
    (( PASS++ )) || true
else
    echo "  FAIL: Collated VQC JSON payload still emits non-finite numeric tokens"
    (( FAIL++ )) || true
fi

# --- Test 19m: Variant pass-rate comparison plot container ---
echo "--- Test 19m: Variant pass-rate comparison plot container ---"
if grep -q 'plot-vqc-pass-compare' "${TMP_DIR}/pipeline_report.html"; then
    echo "  PASS: Variant pass-rate comparison plot container present"
    (( PASS++ )) || true
else
    echo "  FAIL: Variant pass-rate comparison plot container missing"
    (( FAIL++ )) || true
fi

# --- Test 19n: Variant QC toggles can be derived from collated data ---
echo "--- Test 19n: Variant QC toggle fallback to collated groups ---"
if REPO_DIR="${REPO_DIR}" python3 - <<'PY'
import json
import os
import sys

sys.path.insert(0, os.path.join(os.environ['REPO_DIR'], 'scripts'))
from generate_report import _build_html  # noqa: E402

stats = {
    'n_samples': 1,
    'n_pass': 1,
    'n_fail': 0,
    'n_male': 1,
    'n_female': 0,
    'pass_rate': 100.0,
}
tool_versions = {'python3': '3.x'}
collated = json.dumps({
    'all': {'n_variants': 1, 'cr_hist': {'edges': [0.99], 'counts': [1]}, 'hwe_hist': {'edges': [0.3], 'counts': [1]}, 'maf_hist': {'edges': [0.12], 'counts': [1]}},
    'EAS': {'n_variants': 1, 'cr_hist': {'edges': [0.98], 'counts': [1]}, 'hwe_hist': {'edges': [0.1], 'counts': [1]}, 'maf_hist': {'edges': [0.2], 'counts': [1]}},
})
html = _build_html(
    stats, None, {}, '', '', '', '', '',
    '', tool_versions, 'stage2',
    ancestry_vqc_texts={}, collated_vqc_json=collated
)
assert 'value="EAS"' in html
assert 'plot-vqc-cr' in html
assert 'class="vqc-toggle"' in html
PY
then
    echo "  PASS: Variant QC toggles created from collated groups even without text summary"
    (( PASS++ )) || true
else
    echo "  FAIL: Variant QC toggle fallback from collated groups did not work"
    (( FAIL++ )) || true
fi

# --- Test 19k: Cross-ancestry summary JSON computation ---
echo "--- Test 19k: Cross-ancestry summary JSON computation ---"
if REPO_DIR="${REPO_DIR}" TMP_DIR="${TMP_DIR}" python3 - <<'PY'
import json
import os
import sys

repo_dir = os.environ['REPO_DIR']
tmp_dir = os.environ['TMP_DIR']
sys.path.insert(0, os.path.join(repo_dir, 'scripts'))
from generate_report import _prepare_collated_vqc_json  # noqa: E402

test_file = os.path.join(tmp_dir, 'ancestry_stratified_qc', 'collated_variant_qc.tsv')
payload = json.loads(_prepare_collated_vqc_json(test_file))
cross = payload.get('_meta', {}).get('cross_ancestry')
assert cross is not None
assert cross['n_variants'] == 3, cross
assert cross['n_call_rate_pass'] == 2, cross
assert cross['n_hwe_pass'] == 3, cross
assert cross['n_maf_pass'] == 3, cross
assert cross['n_qc_pass'] == 2, cross
PY
then
    echo "  PASS: Cross-ancestry pass counts are computed correctly"
    (( PASS++ )) || true
else
    echo "  FAIL: Cross-ancestry pass counts not computed as expected"
    (( FAIL++ )) || true
fi

# --- Test 19l: Collate script emits cross-ancestry pass flags ---
echo "--- Test 19l: Collate script cross-ancestry flags ---"
if REPO_DIR="${REPO_DIR}" TMP_DIR="${TMP_DIR}" python3 - <<'PY'
import csv
import os
import subprocess
import sys

repo_dir = os.environ['REPO_DIR']
tmp_dir = os.environ['TMP_DIR']
base = os.path.join(tmp_dir, 'collate_flags_test')
os.makedirs(base, exist_ok=True)

def write_qc(prefix, rows):
    qc_dir = os.path.join(base, prefix)
    os.makedirs(qc_dir, exist_ok=True)
    with open(os.path.join(qc_dir, 'variant_qc.vmiss'), 'w') as f:
        f.write('#CHROM\tID\tMISSING_CT\tOBS_CT\tF_MISS\n')
        for vid, vals in rows.items():
            f.write(f"1\t{vid}\t0\t100\t{vals['fmiss']}\n")
    with open(os.path.join(qc_dir, 'variant_qc.hardy'), 'w') as f:
        f.write('#CHROM\tID\tA1\tA2\tCT\tP\n')
        for vid, vals in rows.items():
            f.write(f"1\t{vid}\tA\tG\t100\t{vals['hwe']}\n")
    with open(os.path.join(qc_dir, 'variant_qc.afreq'), 'w') as f:
        f.write('#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n')
        for vid, vals in rows.items():
            f.write(f"1\t{vid}\tA\tG\t{vals['af']}\t200\n")
    return qc_dir

eur = write_qc('EUR', {
    'rs_pass_all': {'fmiss': '0.01', 'hwe': '1e-4', 'af': '0.20'},
    'rs_missing_in_afr': {'fmiss': '0.01', 'hwe': '0.2', 'af': '0.10'},
})
afr = write_qc('AFR', {
    'rs_pass_all': {'fmiss': '0.01', 'hwe': '0.2', 'af': '0.30'},
    'rs_fail_some': {'fmiss': '0.03', 'hwe': '1e-7', 'af': '0.005'},
})

out_tsv = os.path.join(base, 'collated.tsv')
cmd = [
    sys.executable,
    os.path.join(repo_dir, 'scripts', 'collate_variant_qc.py'),
    '--ancestry-variant-qc', f'EUR:{eur}',
    '--ancestry-variant-qc', f'AFR:{afr}',
    '--output', out_tsv,
]
subprocess.check_call(cmd)

with open(out_tsv) as f:
    rows = list(csv.DictReader(f, delimiter='\t'))

flags = {row['variant_id']: row for row in rows}
assert flags['rs_pass_all']['all_ancestries_call_rate_pass'] == '1'
assert flags['rs_pass_all']['all_ancestries_hwe_pass'] == '1'
assert flags['rs_pass_all']['all_ancestries_maf_pass'] == '1'
assert flags['rs_pass_all']['all_ancestries_qc_pass'] == '1'

assert flags['rs_fail_some']['all_ancestries_call_rate_pass'] == '0'
assert flags['rs_fail_some']['all_ancestries_hwe_pass'] == '0'
assert flags['rs_fail_some']['all_ancestries_maf_pass'] == '0'
assert flags['rs_fail_some']['all_ancestries_qc_pass'] == '0'

assert flags['rs_missing_in_afr']['all_ancestries_qc_pass'] == '0'
PY
then
    echo "  PASS: Collated TSV includes correct cross-ancestry pass flags"
    (( PASS++ )) || true
else
    echo "  FAIL: Collated TSV cross-ancestry pass flags incorrect"
    (( FAIL++ )) || true
fi

echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

if [[ "${FAIL}" -gt 0 ]]; then
    exit 1
fi
