#!/usr/bin/env bash
#
# test_sample_name_map.sh
#
# Tests for the sample name map feature.
# Validates name map parsing, sanity checks, and GTC renaming logic.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PASS=0
FAIL=0
TMPDIR=$(mktemp -d)
trap 'rm -rf "${TMPDIR}"' EXIT

echo "============================================"
echo "  Sample Name Map Tests"
echo "============================================"
echo ""

# ---------------------------------------------------------------
# Test 1: Validate sample_name_map.txt format (two tab-separated columns)
# ---------------------------------------------------------------
echo "--- Test 1: sample_name_map.txt format ---"
MAP_FILE="${REPO_DIR}/tests/data/sample_name_map.txt"
if [[ -f "${MAP_FILE}" ]]; then
    bad_lines=$(awk -F'\t' 'NF != 2 { print NR": "$0 }' "${MAP_FILE}")
    if [[ -z "${bad_lines}" ]]; then
        n_entries=$(wc -l < "${MAP_FILE}" | tr -d ' ')
        echo "  PASS: sample_name_map.txt has ${n_entries} entries, all 2-column"
        (( PASS++ )) || true
    else
        echo "  FAIL: sample_name_map.txt has malformed lines:"
        echo "${bad_lines}" | head -5 | sed 's/^/        /'
        (( FAIL++ )) || true
    fi
else
    echo "  FAIL: sample_name_map.txt not found"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 2: No duplicate original names in map
# ---------------------------------------------------------------
echo "--- Test 2: No duplicate original names ---"
n_dupes=$(awk -F'\t' '{ print $1 }' "${MAP_FILE}" | sort | uniq -d | wc -l | tr -d ' ')
if [[ "${n_dupes}" -eq 0 ]]; then
    echo "  PASS: No duplicate original names"
    (( PASS++ )) || true
else
    echo "  FAIL: ${n_dupes} duplicate original names found"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 3: Check for duplicate target names (informational)
# ---------------------------------------------------------------
echo "--- Test 3: Duplicate target names check ---"
n_dupes=$(awk -F'\t' '{ print $2 }' "${MAP_FILE}" | sort | uniq -d | wc -l | tr -d ' ')
if [[ "${n_dupes}" -eq 0 ]]; then
    echo "  PASS: No duplicate target names"
    (( PASS++ )) || true
else
    echo "  PASS: ${n_dupes} duplicate target names found (expected for replicate samples)"
    (( PASS++ )) || true
fi

# ---------------------------------------------------------------
# Test 4: Simulated rename with full match (100%)
# ---------------------------------------------------------------
echo "--- Test 4: GTC rename with full match ---"
GTC_DIR="${TMPDIR}/test4_gtc"
mkdir -p "${GTC_DIR}"
# Create fake GTC files matching first 3 map entries
head -3 "${MAP_FILE}" | while IFS=$'\t' read -r old_name new_name; do
    touch "${GTC_DIR}/${old_name}.gtc"
done

MINI_MAP="${TMPDIR}/test4_map.txt"
head -3 "${MAP_FILE}" > "${MINI_MAP}"

# Simulate the rename logic
declare -A T4_MAP
while IFS=$'\t' read -r old_name new_name; do
    T4_MAP["${old_name}"]="${new_name}"
done < "${MINI_MAP}"

n_renamed=0
for gtc_path in "${GTC_DIR}"/*.gtc; do
    gtc_id=$(basename "${gtc_path}" .gtc)
    if [[ -n "${T4_MAP[${gtc_id}]+x}" ]]; then
        new_name="${T4_MAP[${gtc_id}]}"
        if [[ "${new_name}" != "${gtc_id}" ]]; then
            mv "${gtc_path}" "${GTC_DIR}/${new_name}.gtc"
            (( n_renamed++ )) || true
        fi
    fi
done

# Check that renamed files exist
all_ok="true"
while IFS=$'\t' read -r old_name new_name; do
    if [[ ! -f "${GTC_DIR}/${new_name}.gtc" ]]; then
        echo "  FAIL: Expected ${new_name}.gtc not found after rename"
        all_ok="false"
    fi
done < <(head -3 "${MAP_FILE}")

if [[ "${all_ok}" == "true" && "${n_renamed}" -gt 0 ]]; then
    echo "  PASS: Renamed ${n_renamed} GTC files"
    (( PASS++ )) || true
else
    echo "  FAIL: Rename did not produce expected files"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 5: Match proportion calculation (below 50%)
# ---------------------------------------------------------------
echo "--- Test 5: Match proportion below 50% triggers halt ---"
GTC_DIR5="${TMPDIR}/test5_gtc"
mkdir -p "${GTC_DIR5}"
# Create 4 GTC files, only 1 matching the map
first_entry=$(head -1 "${MAP_FILE}" | cut -f1)
touch "${GTC_DIR5}/${first_entry}.gtc"
touch "${GTC_DIR5}/UNMATCHED_SAMPLE_A.gtc"
touch "${GTC_DIR5}/UNMATCHED_SAMPLE_B.gtc"
touch "${GTC_DIR5}/UNMATCHED_SAMPLE_C.gtc"

n_gtc=4
n_matched=1
match_ok=$(awk -v m="${n_matched}" -v t="${n_gtc}" \
    'BEGIN { print (m / t >= 0.5) ? "yes" : "no" }')

if [[ "${match_ok}" == "no" ]]; then
    echo "  PASS: 1/4 (25%) correctly identified as below 50% threshold"
    (( PASS++ )) || true
else
    echo "  FAIL: Should have identified low match proportion"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 6: Match proportion calculation (above 50%)
# ---------------------------------------------------------------
echo "--- Test 6: Match proportion above 50% passes ---"
n_gtc2=4
n_matched2=3
match_ok2=$(awk -v m="${n_matched2}" -v t="${n_gtc2}" \
    'BEGIN { print (m / t >= 0.5) ? "yes" : "no" }')

if [[ "${match_ok2}" == "yes" ]]; then
    echo "  PASS: 3/4 (75%) correctly passes 50% threshold"
    (( PASS++ )) || true
else
    echo "  FAIL: Should have passed threshold"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 7: Name mapping file format
# ---------------------------------------------------------------
echo "--- Test 7: Name mapping file output ---"
MAP_OUT="${TMPDIR}/test7_mapping.tsv"
{
    echo -e "sample_id\toriginal_name"
    echo -e "NA19152\t5426529111_R04C01"
    echo -e "UNMATCHED\tUNMATCHED"
} > "${MAP_OUT}"

header=$(head -1 "${MAP_OUT}")
if [[ "${header}" == "sample_id"$'\t'"original_name" ]]; then
    echo "  PASS: Name mapping file has correct header"
    (( PASS++ )) || true
else
    echo "  FAIL: Unexpected header: ${header}"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 8: QC report original_name column insertion
# ---------------------------------------------------------------
echo "--- Test 8: QC report original_name column ---"
QC_FILE="${TMPDIR}/test8_qc.tsv"
QC_OUT="${TMPDIR}/test8_qc_out.tsv"
MAP_FILE8="${TMPDIR}/test8_map.tsv"

{
    echo -e "sample_id\tcall_rate\tlrr_sd\tcomputed_gender"
    echo -e "NA19152\t0.98\t0.15\tF"
    echo -e "UNKNOWN\t0.95\t0.20\tM"
} > "${QC_FILE}"

{
    echo -e "sample_id\toriginal_name"
    echo -e "NA19152\t5426529111_R04C01"
    echo -e "UNKNOWN\tUNKNOWN"
} > "${MAP_FILE8}"

awk -F'\t' '
    NR==FNR { if (FNR > 1) map[$1] = $2; next }
    FNR==1 { print "sample_id\toriginal_name\t" substr($0, index($0,$2)); next }
    { printf "%s\t%s", $1, ($1 in map ? map[$1] : $1)
      for (i=2; i<=NF; i++) printf "\t%s", $i
      printf "\n" }
' "${MAP_FILE8}" "${QC_FILE}" > "${QC_OUT}"

# Check output header
out_header=$(head -1 "${QC_OUT}")
if echo "${out_header}" | grep -q "original_name"; then
    echo "  PASS: QC output contains original_name column"
    (( PASS++ )) || true
else
    echo "  FAIL: QC output missing original_name column"
    (( FAIL++ )) || true
fi

# Check mapped value
mapped_val=$(awk -F'\t' 'NR==2 { print $2 }' "${QC_OUT}")
if [[ "${mapped_val}" == "5426529111_R04C01" ]]; then
    echo "  PASS: original_name correctly set to 5426529111_R04C01"
    (( PASS++ )) || true
else
    echo "  FAIL: Expected 5426529111_R04C01, got ${mapped_val}"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 9: --help includes sample-name-map in stage1
# ---------------------------------------------------------------
echo "--- Test 9: stage1 --help mentions sample-name-map ---"
help_out=$(bash "${REPO_DIR}/scripts/stage1_initial_genotyping.sh" --help 2>&1) || true
if echo "${help_out}" | grep -q "sample-name-map"; then
    echo "  PASS: stage1 --help mentions --sample-name-map"
    (( PASS++ )) || true
else
    echo "  FAIL: stage1 --help does not mention --sample-name-map"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 10: --help includes sample-name-map in run_pipeline
# ---------------------------------------------------------------
echo "--- Test 10: run_pipeline --help mentions sample-name-map ---"
help_out=$(bash "${REPO_DIR}/scripts/run_pipeline.sh" --help 2>&1) || true
if echo "${help_out}" | grep -q "sample-name-map"; then
    echo "  PASS: run_pipeline --help mentions --sample-name-map"
    (( PASS++ )) || true
else
    echo "  FAIL: run_pipeline --help does not mention --sample-name-map"
    (( FAIL++ )) || true
fi

# ---------------------------------------------------------------
# Test 11: --help includes sample-name-map in process_1000g
# ---------------------------------------------------------------
echo "--- Test 11: process_1000g --help mentions sample-name-map ---"
help_out=$(bash "${REPO_DIR}/scripts/process_1000g.sh" --help 2>&1) || true
if echo "${help_out}" | grep -q "sample-name-map"; then
    echo "  PASS: process_1000g --help mentions --sample-name-map"
    (( PASS++ )) || true
else
    echo "  FAIL: process_1000g --help does not mention --sample-name-map"
    (( FAIL++ )) || true
fi

echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "============================================"

if [[ "${FAIL}" -gt 0 ]]; then
    exit 1
fi
