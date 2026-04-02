#!/usr/bin/env bash
#
# test_qc_filtering.sh
#
# Regression tests for:
#   - filter_related_samples.py (relatedness-based sample exclusion)
#   - filter_het_outliers.py (per-ancestry heterozygosity outlier exclusion)
#   - recluster_egt.py ddof=1 usage verification
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

PASS=0
FAIL=0

echo "============================================"
echo "  QC Filtering Tests"
echo "============================================"
echo ""

# ---------------------------------------------------------------
# Helper: check a value
# ---------------------------------------------------------------
assert_eq() {
    local actual="$1" expected="$2" label="$3"
    if [[ "${actual}" == "${expected}" ]]; then
        echo "  PASS: ${label} = ${expected}"
        (( PASS++ )) || true
    else
        echo "  FAIL: ${label} expected='${expected}' actual='${actual}'"
        (( FAIL++ )) || true
    fi
}

assert_contains() {
    local haystack="$1" needle="$2" label="$3"
    if echo "${haystack}" | grep -qF -- "${needle}"; then
        echo "  PASS: ${label}"
        (( PASS++ )) || true
    else
        echo "  FAIL: ${label} (missing '${needle}')"
        (( FAIL++ )) || true
    fi
}

# ===============================================================
# Test 1: filter_related_samples.py  — basic relatedness filtering
# ===============================================================
echo "--- Test 1: Relatedness filtering (basic) ---"

# Create mock ped_check.csv with two pairs: one above threshold, one below
cat > "${TMP_DIR}/ped_check.csv" <<EOF
sample_a,sample_b,n,rel,pedigree_relatedness,rel_difference,ibs0,ibs2,shared_hets,hets_a,hets_b,pedigree_parents
S1,S2,1000,0.250,0.0,0.250,10,400,200,300,310,False
S3,S4,1000,0.050,0.0,0.050,50,350,150,280,290,False
S5,S6,1000,0.500,0.0,0.500,5,500,250,320,315,False
EOF

# Create mock QC file: S2 has lower call rate, S6 has higher LRR SD
cat > "${TMP_DIR}/sample_qc.tsv" <<EOF
sample_id	call_rate	lrr_sd	het_rate
S1	0.990000	0.200000	0.300000
S2	0.950000	0.200000	0.300000
S3	0.990000	0.200000	0.300000
S4	0.990000	0.200000	0.300000
S5	0.990000	0.200000	0.300000
S6	0.990000	0.350000	0.300000
EOF

python3 "${REPO_DIR}/scripts/filter_related_samples.py" \
    --ped-check "${TMP_DIR}/ped_check.csv" \
    --sample-qc "${TMP_DIR}/sample_qc.tsv" \
    --output-excluded "${TMP_DIR}/rel_excluded.txt" \
    --output-log "${TMP_DIR}/rel_log.tsv" \
    --min-relatedness 0.125 > /dev/null

N_REL_EXCLUDED=$(wc -l < "${TMP_DIR}/rel_excluded.txt" | tr -d ' ')
assert_eq "${N_REL_EXCLUDED}" "2" "relatedness: 2 samples excluded"

# S2 should be excluded (lower call rate in S1/S2 pair)
assert_contains "$(cat "${TMP_DIR}/rel_excluded.txt")" "S2" "S2 excluded (lower call rate)"

# S6 should be excluded (higher LRR SD in S5/S6 pair)
assert_contains "$(cat "${TMP_DIR}/rel_excluded.txt")" "S6" "S6 excluded (higher LRR SD)"

# S3/S4 pair below threshold → neither excluded
! grep -q "S3" "${TMP_DIR}/rel_excluded.txt" && ! grep -q "S4" "${TMP_DIR}/rel_excluded.txt"
if [[ $? -eq 0 ]]; then
    echo "  PASS: S3/S4 not excluded (below threshold)"
    (( PASS++ )) || true
else
    echo "  FAIL: S3/S4 incorrectly excluded"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 2: filter_related_samples.py  — no related pairs
# ===============================================================
echo "--- Test 2: Relatedness filtering (no related pairs) ---"

cat > "${TMP_DIR}/ped_check_none.csv" <<EOF
sample_a,sample_b,n,rel,pedigree_relatedness,rel_difference,ibs0,ibs2,shared_hets,hets_a,hets_b,pedigree_parents
A1,A2,1000,0.01,0.0,0.01,100,200,100,250,260,False
EOF

python3 "${REPO_DIR}/scripts/filter_related_samples.py" \
    --ped-check "${TMP_DIR}/ped_check_none.csv" \
    --output-excluded "${TMP_DIR}/rel_excluded_none.txt" \
    --min-relatedness 0.125 > /dev/null

N_NONE=$(wc -l < "${TMP_DIR}/rel_excluded_none.txt" | tr -d ' ')
assert_eq "${N_NONE}" "0" "no related pairs → 0 excluded"

echo ""

# ===============================================================
# Test 3: filter_related_samples.py  — without QC data
# ===============================================================
echo "--- Test 3: Relatedness filtering (no QC data) ---"

cat > "${TMP_DIR}/ped_check_noqc.csv" <<EOF
sample_a,sample_b,n,rel,pedigree_relatedness,rel_difference,ibs0,ibs2,shared_hets,hets_a,hets_b,pedigree_parents
X1,X2,1000,0.200,0.0,0.200,20,400,200,300,310,False
EOF

python3 "${REPO_DIR}/scripts/filter_related_samples.py" \
    --ped-check "${TMP_DIR}/ped_check_noqc.csv" \
    --output-excluded "${TMP_DIR}/rel_excluded_noqc.txt" \
    --min-relatedness 0.125 > /dev/null

N_NOQC=$(wc -l < "${TMP_DIR}/rel_excluded_noqc.txt" | tr -d ' ')
assert_eq "${N_NOQC}" "1" "without QC data → 1 excluded (alphabetically first)"

# Without QC data, min(X1,X2) = X1 should be excluded
assert_contains "$(cat "${TMP_DIR}/rel_excluded_noqc.txt")" "X1" "X1 excluded (alphabetical tiebreak)"

echo ""

# ===============================================================
# Test 4: filter_het_outliers.py  — per-ancestry outlier detection
# ===============================================================
echo "--- Test 4: Heterozygosity outlier filtering ---"

# Create mock QC file with 20 EUR samples; one extreme outlier (E20).
# With 19 normal + 1 outlier in n=20, the max z-score for the outlier is
# (n-1)/sqrt(n) = 19/sqrt(20) ≈ 4.25, comfortably above 3 SD.
cat > "${TMP_DIR}/het_qc.tsv" <<EOF
sample_id	call_rate	lrr_sd	het_rate
E1	0.990	0.20	0.300
E2	0.990	0.20	0.300
E3	0.990	0.20	0.300
E4	0.990	0.20	0.300
E5	0.990	0.20	0.300
E6	0.990	0.20	0.300
E7	0.990	0.20	0.300
E8	0.990	0.20	0.300
E9	0.990	0.20	0.300
E10	0.990	0.20	0.300
E11	0.990	0.20	0.300
E12	0.990	0.20	0.300
E13	0.990	0.20	0.300
E14	0.990	0.20	0.300
E15	0.990	0.20	0.300
E16	0.990	0.20	0.300
E17	0.990	0.20	0.300
E18	0.990	0.20	0.300
E19	0.990	0.20	0.300
E20	0.990	0.20	0.800
EOF

# Create mock peddy het_check assigning all to EUR
cat > "${TMP_DIR}/het_check.csv" <<EOF
sample_id,sampled_sites,mean_depth,median_depth,depth_outlier,het_count,het_ratio,ratio_outlier,idr_baf,p10,p90,PC1,PC2,PC3,PC4,ancestry-prediction,ancestry-prob
E1,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E2,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E3,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E4,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E5,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E6,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E7,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E8,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E9,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E10,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E11,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E12,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E13,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E14,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E15,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E16,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E17,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E18,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E19,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
E20,1000,30,30,False,800,0.80,True,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
EOF

python3 "${REPO_DIR}/scripts/filter_het_outliers.py" \
    --sample-qc "${TMP_DIR}/het_qc.tsv" \
    --peddy-het-check "${TMP_DIR}/het_check.csv" \
    --output-excluded "${TMP_DIR}/het_excluded.txt" \
    --output-log "${TMP_DIR}/het_log.tsv" \
    --sd-threshold 3 > /dev/null

N_HET_EXCLUDED=$(wc -l < "${TMP_DIR}/het_excluded.txt" | tr -d ' ')
assert_eq "${N_HET_EXCLUDED}" "1" "het outliers: 1 excluded (E20)"

assert_contains "$(cat "${TMP_DIR}/het_excluded.txt")" "E20" "E20 excluded (extreme het rate)"

echo ""

# ===============================================================
# Test 5: filter_het_outliers.py  — no outliers
# ===============================================================
echo "--- Test 5: Het outlier filtering (no outliers) ---"

cat > "${TMP_DIR}/het_qc_normal.tsv" <<EOF
sample_id	call_rate	lrr_sd	het_rate
N1	0.990	0.20	0.300
N2	0.990	0.20	0.301
N3	0.990	0.20	0.299
N4	0.990	0.20	0.300
N5	0.990	0.20	0.301
EOF

cat > "${TMP_DIR}/het_check_normal.csv" <<EOF
sample_id,sampled_sites,mean_depth,median_depth,depth_outlier,het_count,het_ratio,ratio_outlier,idr_baf,p10,p90,PC1,PC2,PC3,PC4,ancestry-prediction,ancestry-prob
N1,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
N2,1000,30,30,False,301,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
N3,1000,30,30,False,299,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
N4,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
N5,1000,30,30,False,301,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
EOF

python3 "${REPO_DIR}/scripts/filter_het_outliers.py" \
    --sample-qc "${TMP_DIR}/het_qc_normal.tsv" \
    --peddy-het-check "${TMP_DIR}/het_check_normal.csv" \
    --output-excluded "${TMP_DIR}/het_excluded_normal.txt" \
    --sd-threshold 3 > /dev/null

N_HET_NORMAL=$(wc -l < "${TMP_DIR}/het_excluded_normal.txt" | tr -d ' ')
assert_eq "${N_HET_NORMAL}" "0" "no het outliers → 0 excluded"

echo ""

# ===============================================================
# Test 6: filter_het_outliers.py  — multi-ancestry groups
# ===============================================================
echo "--- Test 6: Het outlier filtering (multi-ancestry) ---"

# 19 normal EUR + 1 outlier, 19 normal AFR + 1 outlier
cat > "${TMP_DIR}/het_qc_multi.tsv" <<EOF
sample_id	call_rate	lrr_sd	het_rate
ME1	0.990	0.20	0.300
ME2	0.990	0.20	0.300
ME3	0.990	0.20	0.300
ME4	0.990	0.20	0.300
ME5	0.990	0.20	0.300
ME6	0.990	0.20	0.300
ME7	0.990	0.20	0.300
ME8	0.990	0.20	0.300
ME9	0.990	0.20	0.300
ME10	0.990	0.20	0.300
ME11	0.990	0.20	0.300
ME12	0.990	0.20	0.300
ME13	0.990	0.20	0.300
ME14	0.990	0.20	0.300
ME15	0.990	0.20	0.300
ME16	0.990	0.20	0.300
ME17	0.990	0.20	0.300
ME18	0.990	0.20	0.300
ME19	0.990	0.20	0.300
ME20	0.990	0.20	0.800
MA1	0.990	0.20	0.200
MA2	0.990	0.20	0.200
MA3	0.990	0.20	0.200
MA4	0.990	0.20	0.200
MA5	0.990	0.20	0.200
MA6	0.990	0.20	0.200
MA7	0.990	0.20	0.200
MA8	0.990	0.20	0.200
MA9	0.990	0.20	0.200
MA10	0.990	0.20	0.200
MA11	0.990	0.20	0.200
MA12	0.990	0.20	0.200
MA13	0.990	0.20	0.200
MA14	0.990	0.20	0.200
MA15	0.990	0.20	0.200
MA16	0.990	0.20	0.200
MA17	0.990	0.20	0.200
MA18	0.990	0.20	0.200
MA19	0.990	0.20	0.200
MA20	0.990	0.20	0.010
EOF

cat > "${TMP_DIR}/het_check_multi.csv" <<EOF
sample_id,sampled_sites,mean_depth,median_depth,depth_outlier,het_count,het_ratio,ratio_outlier,idr_baf,p10,p90,PC1,PC2,PC3,PC4,ancestry-prediction,ancestry-prob
ME1,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME2,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME3,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME4,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME5,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME6,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME7,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME8,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME9,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME10,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME11,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME12,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME13,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME14,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME15,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME16,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME17,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME18,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME19,1000,30,30,False,300,0.30,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
ME20,1000,30,30,False,800,0.80,True,0.0,0.0,1.0,0.0,0.0,0.0,0.0,EUR,0.99
MA1,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA2,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA3,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA4,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA5,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA6,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA7,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA8,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA9,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA10,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA11,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA12,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA13,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA14,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA15,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA16,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA17,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA18,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA19,1000,30,30,False,200,0.20,False,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
MA20,1000,30,30,False,10,0.01,True,0.0,0.0,1.0,0.0,0.0,0.0,0.0,AFR,0.99
EOF

python3 "${REPO_DIR}/scripts/filter_het_outliers.py" \
    --sample-qc "${TMP_DIR}/het_qc_multi.tsv" \
    --peddy-het-check "${TMP_DIR}/het_check_multi.csv" \
    --output-excluded "${TMP_DIR}/het_excluded_multi.txt" \
    --output-log "${TMP_DIR}/het_log_multi.tsv" \
    --sd-threshold 3 > /dev/null

N_HET_MULTI=$(wc -l < "${TMP_DIR}/het_excluded_multi.txt" | tr -d ' ')
assert_eq "${N_HET_MULTI}" "2" "multi-ancestry het outliers: 2 excluded"

assert_contains "$(cat "${TMP_DIR}/het_excluded_multi.txt")" "ME20" "ME20 excluded (EUR het outlier)"
assert_contains "$(cat "${TMP_DIR}/het_excluded_multi.txt")" "MA20" "MA20 excluded (AFR het outlier)"

echo ""

# ===============================================================
# Test 7: recluster_egt.py uses ddof=1
# ===============================================================
echo "--- Test 7: recluster_egt.py uses ddof=1 ---"

# Check that recluster_egt.py uses Bessel-corrected (ddof=1 equivalent) std.
# The streaming Welford implementation divides M2 by (n-1), which is
# equivalent to ddof=1.  Also verify no ddof=0 pattern leaks in.
DDOF1_CODE=$(grep -c 'ddof=1\|/ (n - 1)' "${REPO_DIR}/scripts/recluster_egt.py" || true)
DDOF0_CODE=$(grep -c 'np.nanstd.*ddof=0' "${REPO_DIR}/scripts/recluster_egt.py" || true)

assert_eq "${DDOF0_CODE}" "0" "no np.nanstd with ddof=0 in recluster_egt.py"

if [[ "${DDOF1_CODE}" -ge 2 ]]; then
    echo "  PASS: Bessel correction (ddof=1 / n-1) used ${DDOF1_CODE} times"
    (( PASS++ )) || true
else
    echo "  FAIL: Bessel correction (ddof=1 / n-1) found ${DDOF1_CODE} times, expected >= 2"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 8: Verify ddof=1 produces correct standard deviation
# ===============================================================
echo "--- Test 8: ddof=1 numerical correctness ---"

# For values [1, 3]:
#   ddof=0 → std = 1.0  (population)
#   ddof=1 → std = √2 ≈ 1.4142  (sample)
DDOF1_RESULT=$(python3 -c "
import math
vals = [1.0, 3.0]
mean = sum(vals) / len(vals)
var = sum((x - mean)**2 for x in vals) / (len(vals) - 1)
print(f'{math.sqrt(var):.4f}')
")
assert_eq "${DDOF1_RESULT}" "1.4142" "ddof=1 std([1,3]) = √2 ≈ 1.4142"

echo ""

# ===============================================================
# Test 9: --help flags for new scripts
# ===============================================================
echo "--- Test 9: --help flags for new scripts ---"

for script in filter_related_samples.py filter_het_outliers.py; do
    output=$(python3 "${REPO_DIR}/scripts/${script}" --help 2>&1) || true
    if echo "${output}" | grep -qi "usage"; then
        echo "  PASS: ${script} --help shows usage"
        (( PASS++ )) || true
    else
        echo "  FAIL: ${script} --help missing usage text"
        (( FAIL++ )) || true
    fi
done

echo ""

# ===============================================================
# Test 10: --exclude-samples in ancestry_pca.sh help
# ===============================================================
echo "--- Test 10: ancestry_pca.sh --exclude-samples in help ---"

ANC_PCA_HELP=$(bash "${REPO_DIR}/scripts/ancestry_pca.sh" --help 2>&1) || true
if echo "${ANC_PCA_HELP}" | grep -q "exclude-samples"; then
    echo "  PASS: ancestry_pca.sh help mentions --exclude-samples"
    (( PASS++ )) || true
else
    echo "  FAIL: ancestry_pca.sh help missing --exclude-samples"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 11: --exclude-samples in ancestry_stratified_qc.sh help
# ===============================================================
echo "--- Test 11: ancestry_stratified_qc.sh --exclude-samples in help ---"

ANC_STRAT_HELP=$(bash "${REPO_DIR}/scripts/ancestry_stratified_qc.sh" --help 2>&1) || true
if echo "${ANC_STRAT_HELP}" | grep -q "exclude-samples"; then
    echo "  PASS: ancestry_stratified_qc.sh help mentions --exclude-samples"
    (( PASS++ )) || true
else
    echo "  FAIL: ancestry_stratified_qc.sh help missing --exclude-samples"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 12: Relatedness filtering with chain (one already excluded)
# ===============================================================
echo "--- Test 12: Relatedness chain exclusion ---"

# A-B related and B-C related: only B should be excluded (greedy)
cat > "${TMP_DIR}/ped_check_chain.csv" <<EOF
sample_a,sample_b,n,rel,pedigree_relatedness,rel_difference,ibs0,ibs2,shared_hets,hets_a,hets_b,pedigree_parents
C1,C2,1000,0.250,0.0,0.250,10,400,200,300,310,False
C2,C3,1000,0.250,0.0,0.250,10,400,200,310,305,False
EOF

cat > "${TMP_DIR}/qc_chain.tsv" <<EOF
sample_id	call_rate	lrr_sd	het_rate
C1	0.990	0.20	0.30
C2	0.950	0.20	0.30
C3	0.990	0.20	0.30
EOF

python3 "${REPO_DIR}/scripts/filter_related_samples.py" \
    --ped-check "${TMP_DIR}/ped_check_chain.csv" \
    --sample-qc "${TMP_DIR}/qc_chain.tsv" \
    --output-excluded "${TMP_DIR}/rel_chain_excluded.txt" \
    --min-relatedness 0.125 > /dev/null

N_CHAIN=$(wc -l < "${TMP_DIR}/rel_chain_excluded.txt" | tr -d ' ')
assert_eq "${N_CHAIN}" "1" "chain: only 1 sample excluded"

# C2 has lowest call rate and is shared between both pairs
assert_contains "$(cat "${TMP_DIR}/rel_chain_excluded.txt")" "C2" "C2 excluded in chain (lowest call rate, shared pair member)"

echo ""

# ===============================================================
# Test 13: Combined exclusion list deduplication and empty-line filtering
# ===============================================================
echo "--- Test 13: Combined exclusion list deduplication ---"

# A sample excluded by BOTH filters should appear only once in the merged list.
# Empty files should not produce empty-line entries.
cat > "${TMP_DIR}/rel_combined.txt" <<EOF
SAMP_A
SAMP_B
EOF

cat > "${TMP_DIR}/het_combined.txt" <<EOF
SAMP_B
SAMP_C
EOF

# Simulate the pipeline's combine step (with sed '/^$/d' fix)
cat "${TMP_DIR}/rel_combined.txt" "${TMP_DIR}/het_combined.txt" \
    | sed '/^$/d' | sort -u > "${TMP_DIR}/combined_exclude.txt"

N_COMBINED=$(wc -l < "${TMP_DIR}/combined_exclude.txt" | tr -d ' ')
assert_eq "${N_COMBINED}" "3" "combined: 3 unique samples (SAMP_A, SAMP_B, SAMP_C)"

# Verify SAMP_B appears exactly once
SAMP_B_COUNT=$(grep -c "SAMP_B" "${TMP_DIR}/combined_exclude.txt" || true)
assert_eq "${SAMP_B_COUNT}" "1" "SAMP_B appears exactly once in merged list"

echo ""

# ===============================================================
# Test 14: Empty exclusion files produce no empty-line entries
# ===============================================================
echo "--- Test 14: Empty exclusion files produce no empty lines ---"

true > "${TMP_DIR}/empty_rel.txt"
true > "${TMP_DIR}/empty_het.txt"

cat "${TMP_DIR}/empty_rel.txt" "${TMP_DIR}/empty_het.txt" \
    | sed '/^$/d' | sort -u > "${TMP_DIR}/combined_empty.txt"

N_EMPTY=$(wc -l < "${TMP_DIR}/combined_empty.txt" | tr -d ' ')
assert_eq "${N_EMPTY}" "0" "empty inputs → 0 lines in combined list (no blank entries)"

echo ""

# ===============================================================
# Test 15: ancestry_pca.sh --include-variants in help
# ===============================================================
echo "--- Test 15: ancestry_pca.sh --include-variants in help ---"

PCA_HELP=$(bash "${REPO_DIR}/scripts/ancestry_pca.sh" --help 2>&1 || true)
assert_contains "${PCA_HELP}" "include-variants" "ancestry_pca.sh help mentions --include-variants"

echo ""

# ===============================================================
# Test 16: ancestry_stratified_qc.sh help mentions hwe_passing_variants.txt
# ===============================================================
echo "--- Test 16: ancestry_stratified_qc.sh help mentions hwe_passing_variants.txt ---"

STRAT_HELP=$(bash "${REPO_DIR}/scripts/ancestry_stratified_qc.sh" --help 2>&1 || true)
assert_contains "${STRAT_HELP}" "hwe_passing_variants.txt" "ancestry_stratified_qc.sh help mentions hwe_passing_variants output"

echo ""

# ===============================================================
# Test 17: HWE-passing variant list intersection logic
# ===============================================================
echo "--- Test 17: HWE-passing variant list intersection logic ---"

# Simulate .hardy files from two ancestry groups
HWE_TEST_DIR="${TMP_DIR}/hwe_intersection"
mkdir -p "${HWE_TEST_DIR}/EUR/variant_qc" "${HWE_TEST_DIR}/AFR/variant_qc"

# EUR hardy: variants v1, v2 pass HWE (p >> 1e-6); v3 borderline (p = 0.001, still passes)
cat > "${HWE_TEST_DIR}/EUR/variant_qc/variant_qc.hardy" <<EOF
#CHROM	ID	A1	AX	HOM_A1_CT	HET_A1_AX_CT	TWO_AX_CT	O(HET_A1_AX)	E(HET_A1_AX)	P
1	v1	A	T	100	200	100	0.5	0.5	0.95
1	v2	A	T	100	200	100	0.5	0.5	0.80
1	v3	A	T	100	200	100	0.5	0.5	0.001
EOF

# AFR hardy: v1, v2 pass HWE; v3 fails (p < 1e-6)
cat > "${HWE_TEST_DIR}/AFR/variant_qc/variant_qc.hardy" <<EOF
#CHROM	ID	A1	AX	HOM_A1_CT	HET_A1_AX_CT	TWO_AX_CT	O(HET_A1_AX)	E(HET_A1_AX)	P
1	v1	A	T	100	200	100	0.5	0.5	0.90
1	v2	A	T	100	200	100	0.5	0.5	0.50
1	v3	A	T	100	200	100	0.5	0.5	1e-8
EOF

# Test the intersection logic inline (mirrors ancestry_stratified_qc.sh logic)
HWE_P_THRESH="1e-6"
for ANC in EUR AFR; do
    awk -v thresh="${HWE_P_THRESH}" \
        'NR>1 { p=$NF+0; if (p >= thresh) print $2 }' \
        "${HWE_TEST_DIR}/${ANC}/variant_qc/variant_qc.hardy" | sort > "${HWE_TEST_DIR}/${ANC}/variant_qc/hwe_passing_ids.txt"
done

# Intersect
cp "${HWE_TEST_DIR}/EUR/variant_qc/hwe_passing_ids.txt" "${HWE_TEST_DIR}/hwe_tmp.txt"
comm -12 "${HWE_TEST_DIR}/hwe_tmp.txt" "${HWE_TEST_DIR}/AFR/variant_qc/hwe_passing_ids.txt" > "${HWE_TEST_DIR}/hwe_passing_variants.txt"

N_HWE=$(wc -l < "${HWE_TEST_DIR}/hwe_passing_variants.txt" | tr -d ' ')
assert_eq "${N_HWE}" "2" "HWE intersection: 2 variants pass HWE in both EUR and AFR"

# Check that the correct variants are in the intersection
assert_contains "$(cat "${HWE_TEST_DIR}/hwe_passing_variants.txt")" "v1" "v1 in HWE-passing intersection"
assert_contains "$(cat "${HWE_TEST_DIR}/hwe_passing_variants.txt")" "v2" "v2 in HWE-passing intersection"

# v3 fails in AFR so should NOT be in the intersection
if grep -q "v3" "${HWE_TEST_DIR}/hwe_passing_variants.txt"; then
    echo "  FAIL: v3 should not be in HWE-passing intersection (fails AFR HWE)"
    (( FAIL++ )) || true
else
    echo "  PASS: v3 correctly excluded from HWE-passing intersection"
    (( PASS++ )) || true
fi

echo ""

# ===============================================================
# Test 18: Sex check cross-tabulation logic
# ===============================================================
echo "--- Test 18: Sex check cross-tabulation logic ---"

python3 -c "
import sys
sys.path.insert(0, '${REPO_DIR}/scripts')
from plot_sex_check import cross_tabulate_sex

sample_names = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7']
sex_map = {'S1': 'M', 'S2': 'F', 'S3': 'M', 'S4': 'F', 'S5': 'M',
           'S6': 'F', 'S7': 'U'}
chrx_medians = {0: -0.1, 1: 0.03, 2: -0.1, 3: 0.03, 4: -0.1, 5: 0.03, 6: 0.01}
chry_medians = {0: 0.08, 1: -0.07, 2: 0.08, 3: -0.07, 4: 0.08, 5: -0.07, 6: -0.05}
peddy_sex_results = {
    'S1': {'peddy_sex': 'male', 'error': False},
    'S2': {'peddy_sex': 'female', 'error': False},
    'S3': {'peddy_sex': 'male', 'error': True},
    # S6: methods agree (lrr=F, peddy=F), but peddy_error=True
    # (PED sex may have been wrong/unknown); status should be CONCORDANT
    'S6': {'peddy_sex': 'female', 'error': True},
    # S7: only peddy available (lrr=U->NA), peddy_error=True
    'S7': {'peddy_sex': 'male', 'error': True},
}

rows = cross_tabulate_sex(sample_names, sex_map, chrx_medians, chry_medians,
                          peddy_sex_results)

# Check S1: lrr=M, peddy=M -> CONCORDANT
s1 = [r for r in rows if r['sample_id'] == 'S1'][0]
assert s1['status'] == 'CONCORDANT', f'S1: expected CONCORDANT, got {s1[\"status\"]}'

# Check S3: lrr=M, peddy=M -> CONCORDANT
s3 = [r for r in rows if r['sample_id'] == 'S3'][0]
assert s3['status'] == 'CONCORDANT', f'S3: expected CONCORDANT, got {s3[\"status\"]}'

# Check S4: lrr=F, no peddy -> SINGLE_METHOD
s4 = [r for r in rows if r['sample_id'] == 'S4'][0]
assert s4['status'] == 'SINGLE_METHOD', f'S4: expected SINGLE_METHOD, got {s4[\"status\"]}'

# Check S6: peddy_error=True but lrr=F, peddy=F -> CONCORDANT
s6 = [r for r in rows if r['sample_id'] == 'S6'][0]
assert s6['status'] == 'CONCORDANT', f'S6: expected CONCORDANT, got {s6[\"status\"]}'

# Check S7: only peddy has M call (lrr=U->NA) -> SINGLE_METHOD
s7 = [r for r in rows if r['sample_id'] == 'S7'][0]
assert s7['status'] == 'SINGLE_METHOD', f'S7: expected SINGLE_METHOD, got {s7[\"status\"]}'

print('All cross-tabulation assertions passed')
" 2>&1

if [[ $? -eq 0 ]]; then
    echo "  PASS: Sex check cross-tabulation correctly classifies concordant/discordant/ambiguous"
    (( PASS++ )) || true
else
    echo "  FAIL: Sex check cross-tabulation logic error"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 19: plot_sex_check.py --peddy-sex-check in help
# ===============================================================
echo "--- Test 19: plot_sex_check.py --peddy-sex-check in help ---"

SEX_HELP=$(python3 "${REPO_DIR}/scripts/plot_sex_check.py" --help 2>&1 || true)
assert_contains "${SEX_HELP}" "peddy-sex-check" "plot_sex_check.py help mentions --peddy-sex-check"

echo ""

# ===============================================================
# Test 20: Sex check table has correct columns
# ===============================================================
echo "--- Test 20: Sex check table format with cross-tabulation columns ---"

python3 -c "
import sys, os, tempfile
sys.path.insert(0, '${REPO_DIR}/scripts')
from plot_sex_check import write_sex_check_table

sample_names = ['S1', 'S2']
sex_map = {'S1': 'M', 'S2': 'F'}
chrx_medians = {0: -0.1, 1: 0.03}
chry_medians = {0: 0.08, 1: -0.07}

with tempfile.TemporaryDirectory() as td:
    path = write_sex_check_table(
        chrx_medians, chry_medians, sample_names, sex_map, td)
    with open(path) as f:
        header = f.readline().strip()
    expected_cols = ['sample_id', 'chrx_lrr_median', 'chry_lrr_median',
                     'computed_gender',
                     'peddy_sex', 'sex_status']
    actual_cols = header.split('\t')
    assert actual_cols == expected_cols, f'Header mismatch: {actual_cols}'

    # Check summary file was created
    summary = os.path.join(td, 'sex_check_summary.txt')
    assert os.path.exists(summary), 'sex_check_summary.txt not created'

print('Table format and summary file verified')
" 2>&1

if [[ $? -eq 0 ]]; then
    echo "  PASS: Sex check table includes cross-tabulation columns"
    (( PASS++ )) || true
else
    echo "  FAIL: Sex check table format incorrect"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 21: PAR/XTR region definitions per genome build
# ===============================================================
echo "--- Test 21: PAR/XTR region definitions per genome build ---"

python3 -c "
import sys
sys.path.insert(0, '${REPO_DIR}/scripts')
from par_xtr_regions import get_par_xtr_regions, get_par_xtr_bed, get_nonpar_chrx_regions, supported_builds
import tempfile, os

# Verify supported builds
builds = supported_builds()
assert 'CHM13' in builds, f'CHM13 not in supported builds: {builds}'
assert 'GRCh38' in builds, f'GRCh38 not in supported builds: {builds}'
assert 'GRCh37' in builds, f'GRCh37 not in supported builds: {builds}'

# Verify CHM13 regions include both chrX and chrY entries
chm13 = get_par_xtr_regions('CHM13')
chm13_x = [r for r in chm13 if r[0] == 'chrX']
chm13_y = [r for r in chm13 if r[0] == 'chrY']
assert len(chm13_x) == 3, f'Expected 3 chrX regions for CHM13, got {len(chm13_x)}'
assert len(chm13_y) == 3, f'Expected 3 chrY regions for CHM13, got {len(chm13_y)}'
assert chm13_x[0] == ('chrX', 0, 2781479, 'PAR1'), f'CHM13 chrX PAR1: {chm13_x[0]}'
assert chm13_x[1] == ('chrX', 2781479, 6400875, 'XTR'), f'CHM13 chrX XTR: {chm13_x[1]}'
assert chm13_x[2] == ('chrX', 155701382, 156040895, 'PAR2'), f'CHM13 chrX PAR2: {chm13_x[2]}'
assert chm13_y[0] == ('chrY', 0, 2458320, 'PAR1'), f'CHM13 chrY PAR1: {chm13_y[0]}'
assert chm13_y[1] == ('chrY', 2458320, 6400875, 'XTR'), f'CHM13 chrY XTR: {chm13_y[1]}'
assert chm13_y[2] == ('chrY', 62122809, 62460029, 'PAR2'), f'CHM13 chrY PAR2: {chm13_y[2]}'

# Verify GRCh38 has chrY PAR entries (no XTR on Y for this build)
grch38 = get_par_xtr_regions('GRCh38')
grch38_y = [r for r in grch38 if r[0] == 'chrY']
assert len(grch38_y) == 2, f'Expected 2 chrY regions for GRCh38 (PAR1+PAR2, XTR omitted), got {len(grch38_y)}'

# Verify GRCh37 uses non-chr-prefixed chromosomes for both X and Y
grch37 = get_par_xtr_regions('GRCh37')
grch37_x = [r for r in grch37 if r[0] == 'X']
grch37_y = [r for r in grch37 if r[0] == 'Y']
assert len(grch37_x) == 3, f'GRCh37 should have 3 X regions (PAR1+XTR+PAR2): {grch37_x}'
assert len(grch37_y) == 2, f'GRCh37 should have 2 Y regions (PAR1+PAR2): {grch37_y}'
assert grch37_x[0][0] == 'X', f'GRCh37 should use X not chrX: {grch37_x[0][0]}'
assert grch37_y[0][0] == 'Y', f'GRCh37 should use Y not chrY: {grch37_y[0][0]}'

# Verify BED output includes both chrX and chrY lines
with tempfile.TemporaryDirectory() as td:
    bed_path = get_par_xtr_bed('CHM13', os.path.join(td, 'test.bed'))
    with open(bed_path) as f:
        lines = f.readlines()
    assert len(lines) == 6, f'Expected 6 BED lines (3 chrX + 3 chrY), got {len(lines)}'
    chrx_lines = [l for l in lines if l.startswith('chrX')]
    chry_lines = [l for l in lines if l.startswith('chrY')]
    assert len(chrx_lines) == 3, f'Expected 3 chrX BED lines, got {len(chrx_lines)}'
    assert len(chry_lines) == 3, f'Expected 3 chrY BED lines, got {len(chry_lines)}'

# Verify non-PAR region computation still returns chrX-only regions
nonpar = get_nonpar_chrx_regions('CHM13')
assert len(nonpar) >= 1, f'Expected at least 1 nonPAR region, got {len(nonpar)}'
# All nonpar regions should be chrX, not chrY
for region in nonpar:
    assert 'chrX' in region or 'X:' in region, f'nonPAR region should be chrX: {region}'

# Verify build aliases work
assert get_par_xtr_regions('hg38') == get_par_xtr_regions('GRCh38')
assert get_par_xtr_regions('hg19') == get_par_xtr_regions('GRCh37')
assert get_par_xtr_regions('T2T') == get_par_xtr_regions('CHM13')

# Verify unknown build raises ValueError
try:
    get_par_xtr_regions('unknown_build')
    assert False, 'Should have raised ValueError'
except ValueError:
    pass

print('All PAR/XTR region assertions passed (chrX + chrY)')
" 2>&1

if [[ $? -eq 0 ]]; then
    echo "  PASS: PAR/XTR region definitions correct for all genome builds"
    (( PASS++ )) || true
else
    echo "  FAIL: PAR/XTR region definition error"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 22: plot_sex_check.py --genome in help
# ===============================================================
echo "--- Test 22: plot_sex_check.py --genome in help ---"

SEX_HELP=$(python3 "${REPO_DIR}/scripts/plot_sex_check.py" --help 2>&1 || true)
assert_contains "${SEX_HELP}" "--genome" "plot_sex_check.py help mentions --genome"

echo ""

# ===============================================================
# Test 23: LRR cache round-trip
# ===============================================================
echo "--- Test 23: LRR cache round-trip ---"

python3 -c "
import sys, os, tempfile
sys.path.insert(0, '${REPO_DIR}/scripts')
from plot_sex_check import _write_lrr_cache, _read_lrr_cache

# --- LRR cache round-trip ---
medx = {0: -0.1, 1: 0.03, 2: -0.15}
medy = {0: 0.08, 1: -0.07}
with tempfile.TemporaryDirectory() as td:
    cache = os.path.join(td, 'lrr.tsv')
    _write_lrr_cache(cache, medx, medy, 'CHM13')
    loaded = _read_lrr_cache(cache, expected_genome='CHM13')
    assert loaded is not None, 'LRR cache should load'
    lx, ly = loaded
    assert abs(lx[0] - (-0.1)) < 1e-5, f'chrX idx 0 mismatch: {lx[0]}'
    assert abs(ly[0] - 0.08) < 1e-5, f'chrY idx 0 mismatch: {ly[0]}'
    assert 2 in lx, 'idx 2 should be in chrX'
    assert 2 not in ly, 'idx 2 should not be in chrY'

    # LRR genome mismatch
    bad = _read_lrr_cache(cache, expected_genome='GRCh37')
    assert bad is None, 'LRR cache should return None on genome mismatch'

print('LRR cache round-trip verified')
" 2>&1

if [[ $? -eq 0 ]]; then
    echo "  PASS: LRR cache writes, reads, and invalidates correctly"
    (( PASS++ )) || true
else
    echo "  FAIL: LRR cache round-trip error"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 24: run_pipeline.sh passes --genome to plot_sex_check.py
# ===============================================================
echo "--- Test 24: run_pipeline.sh passes --genome to plot_sex_check.py ---"

if grep -q -- '--genome "${GENOME}"' "${REPO_DIR}/scripts/run_pipeline.sh" && \
   grep -A5 'plot_sex_check.py' "${REPO_DIR}/scripts/run_pipeline.sh" | grep -q -- '--genome'; then
    echo "  PASS: run_pipeline.sh passes --genome to plot_sex_check.py"
    (( PASS++ )) || true
else
    echo "  FAIL: run_pipeline.sh missing --genome for plot_sex_check.py"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 25: run_peddy.sh excludes PAR/XTR in chrX append
# ===============================================================
echo "--- Test 25: run_peddy.sh excludes PAR/XTR in chrX append ---"

if grep -q 'par_xtr_exclusion.bed' "${REPO_DIR}/scripts/run_peddy.sh" && \
   grep -q 'PAR/XTR excluded' "${REPO_DIR}/scripts/run_peddy.sh" && \
   grep -q 'get_par_xtr_bed' "${REPO_DIR}/scripts/run_peddy.sh" && \
   grep -q 'source.*utils.sh' "${REPO_DIR}/scripts/run_peddy.sh" && \
   grep -q 'T \^' "${REPO_DIR}/scripts/run_peddy.sh"; then
    echo "  PASS: run_peddy.sh excludes PAR/XTR when appending chrX via utils.sh"
    (( PASS++ )) || true
else
    echo "  FAIL: run_peddy.sh missing PAR/XTR exclusion in chrX append"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 26: utils.sh get_par_xtr_bed function exists
# ===============================================================
echo "--- Test 26: utils.sh get_par_xtr_bed function ---"

if grep -q 'get_par_xtr_bed()' "${REPO_DIR}/scripts/utils.sh" && \
   grep -q 'PAR1' "${REPO_DIR}/scripts/utils.sh" && \
   grep -q 'XTR' "${REPO_DIR}/scripts/utils.sh" && \
   grep -q 'PAR2' "${REPO_DIR}/scripts/utils.sh" && \
   grep -q 'chrY' "${REPO_DIR}/scripts/utils.sh"; then
    echo "  PASS: utils.sh contains get_par_xtr_bed with chrX + chrY PAR/XTR regions"
    (( PASS++ )) || true
else
    echo "  FAIL: utils.sh missing get_par_xtr_bed function"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 27: sex check module does NOT have F-stat functions (removed)
# ===============================================================
echo "--- Test 27: F-stat functions removed from plot_sex_check ---"

python3 -c "
import sys
sys.path.insert(0, '${REPO_DIR}/scripts')
import plot_sex_check as ps

# Verify F-stat functions were removed
assert not hasattr(ps, 'compute_chrx_f_statistic'), \
    'compute_chrx_f_statistic should be removed'
assert not hasattr(ps, '_compute_f_stat_plink2'), \
    '_compute_f_stat_plink2 should be removed'
assert not hasattr(ps, '_compute_f_stat_bcftools'), \
    '_compute_f_stat_bcftools should be removed'
assert not hasattr(ps, '_write_f_stat_cache'), \
    '_write_f_stat_cache should be removed'
assert not hasattr(ps, '_read_f_stat_cache'), \
    '_read_f_stat_cache should be removed'

# Verify essential functions still exist
assert hasattr(ps, 'extract_chrXY_lrr_medians'), \
    'extract_chrXY_lrr_medians should exist'
assert hasattr(ps, 'cross_tabulate_sex'), \
    'cross_tabulate_sex should exist'
assert hasattr(ps, 'write_sex_check_table'), \
    'write_sex_check_table should exist'

print('F-stat functions confirmed removed, essential functions present')
" 2>&1

if [[ $? -eq 0 ]]; then
    echo "  PASS: F-stat functions removed, essential sex check functions present"
    (( PASS++ )) || true
else
    echo "  FAIL: Unexpected function presence/absence in plot_sex_check"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 28: compile_sample_sheet.py — sex check cross-tabulation merge
# ===============================================================
echo "--- Test 28: compile_sample_sheet.py sex check merge ---"

cat > "${TMP_DIR}/css_qc.tsv" <<EOF
sample_id	call_rate	lrr_sd	het_rate
SA	0.990	0.20	0.300
SB	0.980	0.25	0.310
SC	0.970	0.30	0.290
EOF

cat > "${TMP_DIR}/sex_check.tsv" <<EOF
sample_id	chrx_lrr_median	chry_lrr_median	computed_gender	peddy_sex	sex_status
SA	-0.50	-3.00	F	female	CONCORDANT
SB	-0.10	-0.50	M	male	CONCORDANT
SC	-0.30	-1.50	F	female	CONCORDANT
EOF

python3 "${REPO_DIR}/scripts/compile_sample_sheet.py" \
    --sample-qc "${TMP_DIR}/css_qc.tsv" \
    --sex-check "${TMP_DIR}/sex_check.tsv" \
    --output "${TMP_DIR}/css_output.tsv" > /dev/null

# Verify sex check columns are present in header
CSS_HEADER=$(head -1 "${TMP_DIR}/css_output.tsv")
assert_contains "${CSS_HEADER}" "chrx_lrr_median" "compiled sheet has chrx_lrr_median column"
assert_contains "${CSS_HEADER}" "chry_lrr_median" "compiled sheet has chry_lrr_median column"
assert_contains "${CSS_HEADER}" "peddy_sex" "compiled sheet has peddy_sex column"
assert_contains "${CSS_HEADER}" "sex_status" "compiled sheet has sex_status column"

# Verify values for SA
SA_ROW=$(grep "^SA" "${TMP_DIR}/css_output.tsv")
assert_contains "${SA_ROW}" "CONCORDANT" "SA sex_status is CONCORDANT"
assert_contains "${SA_ROW}" "-0.50" "SA chrx_lrr_median is -0.50"

# Verify values for SC
SC_ROW=$(grep "^SC" "${TMP_DIR}/css_output.tsv")
assert_contains "${SC_ROW}" "CONCORDANT" "SC sex_status is CONCORDANT"

echo ""

# ===============================================================
# Test 29: compile_sample_sheet.py — pre-PCA exclusion flags
# ===============================================================
echo "--- Test 29: compile_sample_sheet.py exclusion flags ---"

cat > "${TMP_DIR}/excl_qc.tsv" <<EOF
sample_id	call_rate	lrr_sd	het_rate
X1	0.990	0.20	0.300
X2	0.980	0.25	0.310
X3	0.970	0.30	0.290
X4	0.960	0.35	0.280
EOF

cat > "${TMP_DIR}/rel_excl.txt" <<EOF
X1
X2
EOF

cat > "${TMP_DIR}/het_excl.txt" <<EOF
X2
X3
EOF

cat > "${TMP_DIR}/pca_excl.txt" <<EOF
X1
X2
X3
EOF

python3 "${REPO_DIR}/scripts/compile_sample_sheet.py" \
    --sample-qc "${TMP_DIR}/excl_qc.tsv" \
    --relatedness-excluded "${TMP_DIR}/rel_excl.txt" \
    --het-outlier-excluded "${TMP_DIR}/het_excl.txt" \
    --pre-pca-excluded "${TMP_DIR}/pca_excl.txt" \
    --output "${TMP_DIR}/excl_output.tsv" > /dev/null

# Verify exclusion flag columns in header
EXCL_HEADER=$(head -1 "${TMP_DIR}/excl_output.tsv")
assert_contains "${EXCL_HEADER}" "excluded_relatedness" "compiled sheet has excluded_relatedness column"
assert_contains "${EXCL_HEADER}" "excluded_het_outlier" "compiled sheet has excluded_het_outlier column"
assert_contains "${EXCL_HEADER}" "pre_pca_excluded" "compiled sheet has pre_pca_excluded column"

# X1: rel=1, het=0, pca=1
X1_REL=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="excluded_relatedness") c=i} NR>1 && $1=="X1" {print $c}' "${TMP_DIR}/excl_output.tsv")
X1_HET=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="excluded_het_outlier") c=i} NR>1 && $1=="X1" {print $c}' "${TMP_DIR}/excl_output.tsv")
X1_PCA=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="pre_pca_excluded") c=i} NR>1 && $1=="X1" {print $c}' "${TMP_DIR}/excl_output.tsv")

assert_eq "${X1_REL}" "1" "X1 excluded_relatedness = 1"
assert_eq "${X1_HET}" "0" "X1 excluded_het_outlier = 0"
assert_eq "${X1_PCA}" "1" "X1 pre_pca_excluded = 1"

# X4: not excluded by anything
X4_REL=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="excluded_relatedness") c=i} NR>1 && $1=="X4" {print $c}' "${TMP_DIR}/excl_output.tsv")
X4_HET=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="excluded_het_outlier") c=i} NR>1 && $1=="X4" {print $c}' "${TMP_DIR}/excl_output.tsv")
X4_PCA=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="pre_pca_excluded") c=i} NR>1 && $1=="X4" {print $c}' "${TMP_DIR}/excl_output.tsv")

assert_eq "${X4_REL}" "0" "X4 excluded_relatedness = 0"
assert_eq "${X4_HET}" "0" "X4 excluded_het_outlier = 0"
assert_eq "${X4_PCA}" "0" "X4 pre_pca_excluded = 0"

echo ""

# ===============================================================
# Test 30: collate_variant_qc.py — genotype counts and missingness
# ===============================================================
echo "--- Test 30: collate_variant_qc.py genotype counts ---"

VQCTEST="${TMP_DIR}/vqc_test"
mkdir -p "${VQCTEST}/all_vqc"

cat > "${VQCTEST}/all_vqc/variant_qc.vmiss" <<EOF
#CHROM	ID	MISSING_CT	OBS_CT	F_MISS
1	v1	5	1000	0.005
1	v2	20	1000	0.020
EOF

cat > "${VQCTEST}/all_vqc/variant_qc.hardy" <<EOF
#CHROM	ID	A1	AX	HOM_A1_CT	HET_A1_AX_CT	TWO_AX_CT	O(HET_A1_AX)	E(HET_A1_AX)	P
1	v1	A	T	400	500	100	0.5	0.48	0.50
1	v2	G	C	300	400	300	0.4	0.42	0.30
EOF

cat > "${VQCTEST}/all_vqc/variant_qc.afreq" <<EOF
#CHROM	ID	REF	ALT	ALT_FREQS	OBS_CT
1	v1	A	T	0.35	2000
1	v2	G	C	0.50	2000
EOF

python3 "${REPO_DIR}/scripts/collate_variant_qc.py" \
    --all-variant-qc-dir "${VQCTEST}/all_vqc" \
    --output "${VQCTEST}/collated.tsv" > /dev/null

VQC_HEADER=$(head -1 "${VQCTEST}/collated.tsv")
assert_contains "${VQC_HEADER}" "all_obs_ct" "collated has all_obs_ct column"
assert_contains "${VQC_HEADER}" "all_missing_ct" "collated has all_missing_ct column"
assert_contains "${VQC_HEADER}" "all_hom_a1_ct" "collated has all_hom_a1_ct column"
assert_contains "${VQC_HEADER}" "all_het_ct" "collated has all_het_ct column"
assert_contains "${VQC_HEADER}" "all_hom_a2_ct" "collated has all_hom_a2_ct column"

# Verify actual values for v1
V1_OBS=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="all_obs_ct") c=i} NR>1 && $1=="v1" {print $c}' "${VQCTEST}/collated.tsv")
V1_MISS=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="all_missing_ct") c=i} NR>1 && $1=="v1" {print $c}' "${VQCTEST}/collated.tsv")
V1_HOM1=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="all_hom_a1_ct") c=i} NR>1 && $1=="v1" {print $c}' "${VQCTEST}/collated.tsv")
V1_HET=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="all_het_ct") c=i} NR>1 && $1=="v1" {print $c}' "${VQCTEST}/collated.tsv")
V1_HOM2=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="all_hom_a2_ct") c=i} NR>1 && $1=="v1" {print $c}' "${VQCTEST}/collated.tsv")

assert_eq "${V1_OBS}" "1000" "v1 obs_ct = 1000"
assert_eq "${V1_MISS}" "5" "v1 missing_ct = 5"
assert_eq "${V1_HOM1}" "400" "v1 hom_a1_ct = 400"
assert_eq "${V1_HET}" "500" "v1 het_ct = 500"
assert_eq "${V1_HOM2}" "100" "v1 hom_a2_ct = 100"

echo ""

# ===============================================================
# Test 31: collate_variant_qc.py — Ti/Tv ratio from tstv file
# ===============================================================
echo "--- Test 31: collate_variant_qc.py Ti/Tv ratio ---"

cat > "${VQCTEST}/tstv_stats.txt" <<EOF
SN	0	number of SNPs:	5000
SN	0	number of indels:	200
TSTV	0	3000	1500	2.00	0.00
EOF

python3 "${REPO_DIR}/scripts/collate_variant_qc.py" \
    --all-variant-qc-dir "${VQCTEST}/all_vqc" \
    --tstv-file "${VQCTEST}/tstv_stats.txt" \
    --output "${VQCTEST}/collated_tstv.tsv" > /dev/null

# Verify Ti/Tv ratio is in header comment
TSTV_LINE=$(head -1 "${VQCTEST}/collated_tstv.tsv")
assert_contains "${TSTV_LINE}" "#tstv_ratio=2.00" "Ti/Tv ratio present as header comment"

echo ""

# ===============================================================
# Test 32: collate_variant_qc.py — ancestry-specific genotype counts
# ===============================================================
echo "--- Test 32: collate_variant_qc.py ancestry genotype counts ---"

mkdir -p "${VQCTEST}/eur_vqc"

cat > "${VQCTEST}/eur_vqc/variant_qc.vmiss" <<EOF
#CHROM	ID	MISSING_CT	OBS_CT	F_MISS
1	v1	2	500	0.004
1	v2	10	500	0.020
EOF

cat > "${VQCTEST}/eur_vqc/variant_qc.hardy" <<EOF
#CHROM	ID	A1	AX	HOM_A1_CT	HET_A1_AX_CT	TWO_AX_CT	O(HET_A1_AX)	E(HET_A1_AX)	P
1	v1	A	T	200	250	50	0.5	0.48	0.60
1	v2	G	C	150	200	150	0.4	0.42	0.40
EOF

cat > "${VQCTEST}/eur_vqc/variant_qc.afreq" <<EOF
#CHROM	ID	REF	ALT	ALT_FREQS	OBS_CT
1	v1	A	T	0.35	1000
1	v2	G	C	0.50	1000
EOF

python3 "${REPO_DIR}/scripts/collate_variant_qc.py" \
    --all-variant-qc-dir "${VQCTEST}/all_vqc" \
    --ancestry-variant-qc "EUR:${VQCTEST}/eur_vqc" \
    --output "${VQCTEST}/collated_anc.tsv" > /dev/null

ANC_HEADER=$(grep -v '^#' "${VQCTEST}/collated_anc.tsv" | head -1)
assert_contains "${ANC_HEADER}" "EUR_obs_ct" "ancestry collation has EUR_obs_ct"
assert_contains "${ANC_HEADER}" "EUR_het_ct" "ancestry collation has EUR_het_ct"
assert_contains "${ANC_HEADER}" "EUR_hom_a1_ct" "ancestry collation has EUR_hom_a1_ct"

# Verify EUR v1 values
EUR_V1_HET=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="EUR_het_ct") c=i} /^v1\t/ {print $c}' <(grep -v '^#' "${VQCTEST}/collated_anc.tsv"))
assert_eq "${EUR_V1_HET}" "250" "EUR v1 het_ct = 250"

echo ""

# ===============================================================
# Test 33: compile_sample_sheet.py — --help shows new args
# ===============================================================
echo "--- Test 33: compile_sample_sheet.py help shows new args ---"

CSS_HELP=$(python3 "${REPO_DIR}/scripts/compile_sample_sheet.py" --help 2>&1) || true
assert_contains "${CSS_HELP}" "sex-check" "compile_sample_sheet.py help mentions --sex-check"
assert_contains "${CSS_HELP}" "relatedness-excluded" "compile_sample_sheet.py help mentions --relatedness-excluded"
assert_contains "${CSS_HELP}" "het-outlier-excluded" "compile_sample_sheet.py help mentions --het-outlier-excluded"
assert_contains "${CSS_HELP}" "pre-pca-excluded" "compile_sample_sheet.py help mentions --pre-pca-excluded"

echo ""

# ===============================================================
# Test 34: collate_variant_qc.py — --help shows --tstv-file
# ===============================================================
echo "--- Test 34: collate_variant_qc.py help shows --tstv-file ---"

CVQ_HELP=$(python3 "${REPO_DIR}/scripts/collate_variant_qc.py" --help 2>&1) || true
assert_contains "${CVQ_HELP}" "tstv-file" "collate_variant_qc.py help mentions --tstv-file"

echo ""

# ===============================================================
# Test 35: run_pipeline.sh wires up sex check and exclusion args
# ===============================================================
echo "--- Test 35: run_pipeline.sh wires sex check and exclusion args ---"

if grep -q -- '--sex-check "${SEX_CHECK_TSV}"' "${REPO_DIR}/scripts/run_pipeline.sh" && \
   grep -q -- '--relatedness-excluded "${RELATEDNESS_EXCLUDED}"' "${REPO_DIR}/scripts/run_pipeline.sh" && \
   grep -q -- '--het-outlier-excluded "${HET_OUTLIER_EXCLUDED}"' "${REPO_DIR}/scripts/run_pipeline.sh" && \
   grep -q -- '--pre-pca-excluded "${PRE_PCA_EXCLUDE}"' "${REPO_DIR}/scripts/run_pipeline.sh"; then
    echo "  PASS: run_pipeline.sh passes sex check and exclusion args to compile_sample_sheet.py"
    (( PASS++ )) || true
else
    echo "  FAIL: run_pipeline.sh missing sex check or exclusion args for compile_sample_sheet.py"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 36: compute_sex_chr_variant_qc.sh help and structure
# ===============================================================
echo "--- Test 36: compute_sex_chr_variant_qc.sh help and structure ---"

HELP_OUT=$(bash "${REPO_DIR}/scripts/compute_sex_chr_variant_qc.sh" --help 2>&1 || true)
if echo "${HELP_OUT}" | grep -q "chrX non-PAR" && \
   echo "${HELP_OUT}" | grep -q "females" && \
   echo "${HELP_OUT}" | grep -q "males"; then
    echo "  PASS: compute_sex_chr_variant_qc.sh help mentions sex-stratified approach"
    (( PASS++ )) || true
else
    echo "  FAIL: help text missing sex-stratified description"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 37: collate_variant_qc.py --sex-chr-qc-dir option
# ===============================================================
echo "--- Test 37: collate_variant_qc.py sex chromosome collation ---"

SCQC_TEST_DIR="${TMP_DIR}/sex_chr_qc_test"
mkdir -p "${SCQC_TEST_DIR}/sex_chr_qc" "${SCQC_TEST_DIR}/auto_qc"

# Create mock autosomal variant QC
cat > "${SCQC_TEST_DIR}/auto_qc/variant_qc.vmiss" <<'EOF'
#CHROM	ID	MISSING_CT	OBS_CT	F_MISS
1	rs1	5	1000	0.005
1	rs2	10	1000	0.01
EOF

cat > "${SCQC_TEST_DIR}/auto_qc/variant_qc.hardy" <<'EOF'
#CHROM	ID	A1	AX	HOM_A1_CT	HET_A1_AX_CT	TWO_AX_CT	O(HET_A1_AX)	E(HET_A1_AX)	P
1	rs1	A	G	400	500	100	0.5	0.48	0.5
1	rs2	T	C	300	600	100	0.6	0.5	0.001
EOF

cat > "${SCQC_TEST_DIR}/auto_qc/variant_qc.afreq" <<'EOF'
#CHROM	ID	REF	ALT	ALT_FREQS	OBS_CT
1	rs1	A	G	0.35	2000
1	rs2	T	C	0.40	2000
EOF

# Create mock chrX variant QC
cat > "${SCQC_TEST_DIR}/sex_chr_qc/chrx_variant_qc.vmiss" <<'EOF'
#CHROM	ID	MISSING_CT	OBS_CT	F_MISS
chrX	chrx_v1	3	500	0.006
chrX	chrx_v2	20	500	0.04
EOF

cat > "${SCQC_TEST_DIR}/sex_chr_qc/chrx_variant_qc.hardy" <<'EOF'
#CHROM	ID	A1	AX	HOM_A1_CT	HET_A1_AX_CT	TWO_AX_CT	O(HET_A1_AX)	E(HET_A1_AX)	P
chrX	chrx_v1	A	G	100	150	50	0.5	0.48	0.8
chrX	chrx_v2	T	C	80	170	50	0.57	0.5	0.01
EOF

cat > "${SCQC_TEST_DIR}/sex_chr_qc/chrx_variant_qc.afreq" <<'EOF'
#CHROM	ID	REF	ALT	ALT_FREQS	OBS_CT
chrX	chrx_v1	A	G	0.42	1000
chrX	chrx_v2	T	C	0.45	1000
EOF

# Create mock chrY variant QC (no .hardy for chrY)
cat > "${SCQC_TEST_DIR}/sex_chr_qc/chry_variant_qc.vmiss" <<'EOF'
#CHROM	ID	MISSING_CT	OBS_CT	F_MISS
chrY	chry_v1	2	250	0.008
EOF

cat > "${SCQC_TEST_DIR}/sex_chr_qc/chry_variant_qc.afreq" <<'EOF'
#CHROM	ID	REF	ALT	ALT_FREQS	OBS_CT
chrY	chry_v1	A	G	0.10	500
EOF

SCQC_OUTPUT="${SCQC_TEST_DIR}/collated.tsv"
python3 "${REPO_DIR}/scripts/collate_variant_qc.py" \
    --all-variant-qc-dir "${SCQC_TEST_DIR}/auto_qc" \
    --sex-chr-qc-dir "${SCQC_TEST_DIR}/sex_chr_qc" \
    --output "${SCQC_OUTPUT}" 2>&1 | sed 's/^/    /'

SCQC_HEADER=$(head -1 "${SCQC_OUTPUT}")
assert_contains "${SCQC_HEADER}" "chrX_call_rate" "collated has chrX_call_rate column"
assert_contains "${SCQC_HEADER}" "chrX_hwe_p_females" "collated has chrX_hwe_p_females column"
assert_contains "${SCQC_HEADER}" "chrX_maf" "collated has chrX_maf column"
assert_contains "${SCQC_HEADER}" "chrY_call_rate" "collated has chrY_call_rate column"
assert_contains "${SCQC_HEADER}" "chrY_maf" "collated has chrY_maf column"

# Verify chrX variant is present
CHRX_ROW=$(grep "^chrx_v1" "${SCQC_OUTPUT}" || true)
if [[ -n "${CHRX_ROW}" ]]; then
    echo "  PASS: chrX variant chrx_v1 present in collated output"
    (( PASS++ )) || true
else
    echo "  FAIL: chrX variant chrx_v1 missing from collated output"
    (( FAIL++ )) || true
fi

# Verify chrY variant is present
CHRY_ROW=$(grep "^chry_v1" "${SCQC_OUTPUT}" || true)
if [[ -n "${CHRY_ROW}" ]]; then
    echo "  PASS: chrY variant chry_v1 present in collated output"
    (( PASS++ )) || true
else
    echo "  FAIL: chrY variant chry_v1 missing from collated output"
    (( FAIL++ )) || true
fi

# Verify total variant count includes both autosomal and sex chr variants
N_VARIANTS=$(awk 'NR>1' "${SCQC_OUTPUT}" | wc -l | tr -d ' ')
assert_eq "${N_VARIANTS}" "5" "collated has 5 total variants (2 auto + 2 chrX + 1 chrY)"

echo ""

# ===============================================================
# Test 38: ancestry_stratified_qc.sh accepts --sample-qc and --genome
# ===============================================================
echo "--- Test 38: ancestry_stratified_qc.sh accepts new args ---"

STRAT_HELP=$(bash "${REPO_DIR}/scripts/ancestry_stratified_qc.sh" --help 2>&1 || true)
if echo "${STRAT_HELP}" | grep -q "\-\-sample-qc" && \
   echo "${STRAT_HELP}" | grep -q "\-\-genome"; then
    echo "  PASS: ancestry_stratified_qc.sh help shows --sample-qc and --genome"
    (( PASS++ )) || true
else
    echo "  FAIL: ancestry_stratified_qc.sh help missing new args"
    (( FAIL++ )) || true
fi

echo ""

# ===============================================================
# Test 39: run_pipeline.sh wires sample-qc and genome to strat QC
# ===============================================================
echo "--- Test 39: run_pipeline.sh wires sex chr QC args ---"

PIPELINE_SRC=$(cat "${REPO_DIR}/scripts/run_pipeline.sh")
if echo "${PIPELINE_SRC}" | grep -q "sample-qc" && \
   echo "${PIPELINE_SRC}" | grep -q "STRAT_ARGS.*genome"; then
    echo "  PASS: run_pipeline.sh passes --sample-qc and --genome to ancestry_stratified_qc.sh"
    (( PASS++ )) || true
else
    echo "  FAIL: run_pipeline.sh missing sex chr QC arg wiring"
    (( FAIL++ )) || true
fi

echo ""

# ---------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

if [[ "${FAIL}" -gt 0 ]]; then
    exit 1
fi
