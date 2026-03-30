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

sample_names = ['S1', 'S2', 'S3', 'S4', 'S5']
sex_map = {'S1': 'M', 'S2': 'F', 'S3': 'M', 'S4': 'F', 'S5': 'M'}
chrx_medians = {0: -0.1, 1: 0.03, 2: -0.1, 3: 0.03, 4: -0.1}
chry_medians = {0: 0.08, 1: -0.07, 2: 0.08, 3: -0.07, 4: 0.08}
f_stat_results = {
    'S1': {'f_stat': 0.95, 'n_het': 5, 'n_called': 1000, 'f_sex': 'M'},
    'S2': {'f_stat': 0.05, 'n_het': 450, 'n_called': 1000, 'f_sex': 'F'},
    'S3': {'f_stat': 0.10, 'n_het': 400, 'n_called': 1000, 'f_sex': 'F'},  # discordant
    'S4': {'f_stat': 0.50, 'n_het': 200, 'n_called': 1000, 'f_sex': 'ambiguous'},  # ambiguous
    'S5': {'f_stat': 0.90, 'n_het': 10, 'n_called': 1000, 'f_sex': 'M'},
}
peddy_sex_results = {
    'S1': {'peddy_sex': 'male', 'error': False},
    'S2': {'peddy_sex': 'female', 'error': False},
    'S3': {'peddy_sex': 'male', 'error': True},
}

rows = cross_tabulate_sex(sample_names, sex_map, chrx_medians, chry_medians,
                          f_stat_results, peddy_sex_results)

# Check S1: all agree male → CONCORDANT
s1 = [r for r in rows if r['sample_id'] == 'S1'][0]
assert s1['status'] == 'CONCORDANT', f'S1: expected CONCORDANT, got {s1[\"status\"]}'

# Check S3: LRR says M, F says F → DISCORDANT
s3 = [r for r in rows if r['sample_id'] == 'S3'][0]
assert s3['status'] == 'DISCORDANT', f'S3: expected DISCORDANT, got {s3[\"status\"]}'

# Check S4: F-stat ambiguous, LRR says F → AMBIGUOUS
s4 = [r for r in rows if r['sample_id'] == 'S4'][0]
assert s4['status'] == 'AMBIGUOUS', f'S4: expected AMBIGUOUS, got {s4[\"status\"]}'

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
# Test 20: Sex check table has new columns
# ===============================================================
echo "--- Test 20: Sex check table format with F-statistic columns ---"

python3 -c "
import sys, os, tempfile
sys.path.insert(0, '${REPO_DIR}/scripts')
from plot_sex_check import write_sex_check_table

sample_names = ['S1', 'S2']
sex_map = {'S1': 'M', 'S2': 'F'}
chrx_medians = {0: -0.1, 1: 0.03}
chry_medians = {0: 0.08, 1: -0.07}
f_stat_results = {
    'S1': {'f_stat': 0.95, 'n_het': 5, 'n_called': 1000, 'f_sex': 'M'},
    'S2': {'f_stat': 0.05, 'n_het': 450, 'n_called': 1000, 'f_sex': 'F'},
}

with tempfile.TemporaryDirectory() as td:
    path = write_sex_check_table(
        chrx_medians, chry_medians, sample_names, sex_map, td,
        f_stat_results=f_stat_results)
    with open(path) as f:
        header = f.readline().strip()
    expected_cols = ['sample_id', 'chrx_lrr_median', 'chry_lrr_median',
                     'computed_gender', 'chrx_f_stat', 'f_sex',
                     'peddy_sex', 'sex_status']
    actual_cols = header.split('\t')
    assert actual_cols == expected_cols, f'Header mismatch: {actual_cols}'

    # Check summary file was created
    summary = os.path.join(td, 'sex_check_summary.txt')
    assert os.path.exists(summary), 'sex_check_summary.txt not created'

print('Table format and summary file verified')
" 2>&1

if [[ $? -eq 0 ]]; then
    echo "  PASS: Sex check table includes F-statistic and cross-tabulation columns"
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

# Verify CHM13 regions match the reference values from the issue
chm13 = get_par_xtr_regions('CHM13')
assert len(chm13) == 3, f'Expected 3 regions for CHM13, got {len(chm13)}'
assert chm13[0] == ('chrX', 0, 2781479, 'PAR1'), f'CHM13 PAR1: {chm13[0]}'
assert chm13[1] == ('chrX', 2781479, 6400875, 'XTR'), f'CHM13 XTR: {chm13[1]}'
assert chm13[2] == ('chrX', 155701382, 156040895, 'PAR2'), f'CHM13 PAR2: {chm13[2]}'

# Verify GRCh37 uses non-chr-prefixed chromosomes
grch37 = get_par_xtr_regions('GRCh37')
assert grch37[0][0] == 'X', f'GRCh37 should use X not chrX: {grch37[0][0]}'

# Verify BED output
with tempfile.TemporaryDirectory() as td:
    bed_path = get_par_xtr_bed('CHM13', os.path.join(td, 'test.bed'))
    with open(bed_path) as f:
        lines = f.readlines()
    assert len(lines) == 3, f'Expected 3 BED lines, got {len(lines)}'
    assert lines[0].startswith('chrX\t0\t2781479'), f'Line 0: {lines[0]}'

# Verify non-PAR region computation returns chrX regions
nonpar = get_nonpar_chrx_regions('CHM13')
assert len(nonpar) >= 1, f'Expected at least 1 nonPAR region, got {len(nonpar)}'

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

print('All PAR/XTR region assertions passed')
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
# Test 23: F-stat cache round-trip
# ===============================================================
echo "--- Test 23: F-stat cache round-trip ---"

python3 -c "
import sys, os, tempfile
sys.path.insert(0, '${REPO_DIR}/scripts')
from plot_sex_check import _write_f_stat_cache, _read_f_stat_cache

original = {
    'S1': {'f_stat': 0.95, 'n_het': 5, 'n_called': 1000, 'f_sex': 'M'},
    'S2': {'f_stat': 0.05, 'n_het': 450, 'n_called': 1000, 'f_sex': 'F'},
    'S3': {'f_stat': 0.50, 'n_het': 200, 'n_called': 1000, 'f_sex': 'ambiguous'},
}

with tempfile.TemporaryDirectory() as td:
    cache = os.path.join(td, 'cache.tsv')
    _write_f_stat_cache(cache, original)
    loaded = _read_f_stat_cache(cache)

    assert len(loaded) == 3, f'Expected 3 samples, got {len(loaded)}'
    for name in original:
        assert name in loaded, f'Missing sample {name}'
        assert loaded[name]['f_stat'] == original[name]['f_stat'], \
            f'{name} f_stat mismatch: {loaded[name][\"f_stat\"]} != {original[name][\"f_stat\"]}'
        assert loaded[name]['f_sex'] == original[name]['f_sex'], \
            f'{name} f_sex mismatch'

print('Cache round-trip verified')
" 2>&1

if [[ $? -eq 0 ]]; then
    echo "  PASS: F-stat cache writes and reads correctly"
    (( PASS++ )) || true
else
    echo "  FAIL: F-stat cache round-trip error"
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
   grep -q 'T \^' "${REPO_DIR}/scripts/run_peddy.sh"; then
    echo "  PASS: run_peddy.sh excludes PAR/XTR when appending chrX"
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
   grep -q 'PAR2' "${REPO_DIR}/scripts/utils.sh"; then
    echo "  PASS: utils.sh contains get_par_xtr_bed with PAR/XTR regions"
    (( PASS++ )) || true
else
    echo "  FAIL: utils.sh missing get_par_xtr_bed function"
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
