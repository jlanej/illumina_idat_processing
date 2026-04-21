#!/usr/bin/env bash
#
# test_plink2_sex_chr.sh
#
# Integration test for compute_sex_chr_variant_qc.sh using actual
# plink2 and bcftools binaries.  Validates chrX and chrY analysis
# pipeline with mock BCF data including PAR/XTR variants.
#
# Tests exercise the exact code paths that failed in production:
#   - chrX: plink2 PAR detection requiring --split-par
#   - chrY: sample ID matching via #IID header for BCF input
#
# Requires: plink2, bcftools, python3
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

PASS=0
FAIL=0

echo "============================================"
echo "  plink2 Sex Chromosome Integration Tests"
echo "============================================"
echo ""

# ---------------------------------------------------------------
# Helper
# ---------------------------------------------------------------
assert_file_exists() {
    local path="$1" label="$2"
    if [[ -f "${path}" ]]; then
        echo "  PASS: ${label}"
        (( PASS++ )) || true
    else
        echo "  FAIL: ${label} — file not found: ${path}"
        (( FAIL++ )) || true
    fi
}

assert_file_nonempty() {
    local path="$1" label="$2"
    if [[ -s "${path}" ]]; then
        echo "  PASS: ${label}"
        (( PASS++ )) || true
    else
        echo "  FAIL: ${label} — file missing or empty: ${path}"
        (( FAIL++ )) || true
    fi
}

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

assert_ge() {
    local actual="$1" expected="$2" label="$3"
    if [[ "${actual}" -ge "${expected}" ]]; then
        echo "  PASS: ${label} (${actual} >= ${expected})"
        (( PASS++ )) || true
    else
        echo "  FAIL: ${label} expected >= ${expected}, got ${actual}"
        (( FAIL++ )) || true
    fi
}

# ---------------------------------------------------------------
# Pre-flight: check required tools
# ---------------------------------------------------------------
REQUIRE_PLINK2="${REQUIRE_PLINK2:-0}"
REQUIRE_BCFTOOLS="${REQUIRE_BCFTOOLS:-0}"

if ! command -v plink2 &>/dev/null; then
    if [[ "${REQUIRE_PLINK2}" == "1" ]]; then
        echo "FAIL: plink2 not available but REQUIRE_PLINK2=1"
        exit 1
    fi
    echo "SKIP: plink2 not available — install plink2 to run this test"
    exit 0
fi
if ! command -v bcftools &>/dev/null; then
    if [[ "${REQUIRE_BCFTOOLS}" == "1" ]]; then
        echo "FAIL: bcftools not available but REQUIRE_BCFTOOLS=1"
        exit 1
    fi
    echo "SKIP: bcftools not available — install bcftools to run this test"
    exit 0
fi

echo "Tools:"
echo "  plink2:   $(plink2 --version 2>/dev/null | head -1)"
echo "  bcftools: $(bcftools --version 2>/dev/null | head -1)"
echo ""

# =================================================================
# Test 1: Full pipeline — CHM13 genome with chrX PAR and chrY data
# =================================================================
echo "--- Test 1: Full pipeline with chrX + chrY (CHM13) ---"

# Create mock VCF with chrX and chrY variants.
# CHM13 PAR/XTR boundaries:
#   PAR1: chrX 0-2781479      chrY 0-2458320
#   XTR:  chrX 2781479-6400875  chrY 2458320-6400875
#   PAR2: chrX 155701382-156040895  chrY 62122809-62460029
#
# Samples: S1 (male), S2 (female), S3 (male), S4 (female)
cat > "${TMP_DIR}/mock.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=62460029>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chrX	1500000	xpar1_v1	A	G	100	PASS	.	GT	0/1	0/1	0/0	0/1
chrX	2000000	xpar1_v2	C	T	100	PASS	.	GT	0/0	0/1	0/1	0/0
chrX	3500000	xxtr_v1	G	A	100	PASS	.	GT	0/1	0/0	0/1	0/1
chrX	5000000	xxtr_v2	T	C	100	PASS	.	GT	0/0	0/1	0/0	0/1
chrX	10000000	xnonpar_v1	A	T	100	PASS	.	GT	0/0	0/1	0/0	0/0
chrX	50000000	xnonpar_v2	G	C	100	PASS	.	GT	0/0	0/0	0/0	0/1
chrX	100000000	xnonpar_v3	T	A	100	PASS	.	GT	0/0	0/1	0/0	0/1
chrX	155800000	xpar2_v1	C	G	100	PASS	.	GT	0/1	0/1	0/0	0/1
chrX	155900000	xpar2_v2	A	C	100	PASS	.	GT	0/0	0/1	0/1	0/0
chrY	1000000	ypar1_v1	A	C	100	PASS	.	GT	0/0	./.	0/0	./.
chrY	3000000	yxtr_v1	T	G	100	PASS	.	GT	0/0	./.	0/0	./.
chrY	30000000	ynonpar_v1	C	A	100	PASS	.	GT	0/0	./.	0/0	./.
chrY	40000000	ynonpar_v2	G	T	100	PASS	.	GT	0/0	./.	0/0	./.
chrY	50000000	ynonpar_v3	A	G	100	PASS	.	GT	0/0	./.	0/0	./.
EOF

bcftools view -Ob -o "${TMP_DIR}/mock.bcf" "${TMP_DIR}/mock.vcf" 2>/dev/null
bcftools index "${TMP_DIR}/mock.bcf" 2>/dev/null

cat > "${TMP_DIR}/sample_qc.tsv" <<'EOF'
sample_id	call_rate	lrr_sd	computed_gender
S1	0.990	0.15	M
S2	0.985	0.18	F
S3	0.992	0.12	M
S4	0.988	0.20	F
EOF

OUTPUT_DIR="${TMP_DIR}/sex_chr_qc"

# Run the pipeline — should succeed without error
echo "  Running compute_sex_chr_variant_qc.sh..."
if bash "${REPO_DIR}/scripts/compute_sex_chr_variant_qc.sh" \
    --vcf "${TMP_DIR}/mock.bcf" \
    --sample-qc "${TMP_DIR}/sample_qc.tsv" \
    --output-dir "${OUTPUT_DIR}" \
    --genome CHM13 \
    --threads 1 > "${TMP_DIR}/pipeline.log" 2>&1; then
    echo "  PASS: pipeline completed without error"
    (( PASS++ )) || true
else
    echo "  FAIL: pipeline exited with error (rc=$?)"
    cat "${TMP_DIR}/pipeline.log"
    (( FAIL++ )) || true
fi

# --- chrX output files ---
echo ""
echo "  Checking chrX output files..."
assert_file_exists "${OUTPUT_DIR}/chrX/chrX_all.vmiss" "chrX_all.vmiss exists"
assert_file_exists "${OUTPUT_DIR}/chrX/chrX_all.afreq" "chrX_all.afreq exists"
assert_file_exists "${OUTPUT_DIR}/chrX/chrX_female.vmiss" "chrX_female.vmiss exists"
assert_file_exists "${OUTPUT_DIR}/chrX/chrX_female.hardy" "chrX_female.hardy exists"
assert_file_exists "${OUTPUT_DIR}/chrX/chrX_male.vmiss" "chrX_male.vmiss exists"

# --- chrY output files ---
echo ""
echo "  Checking chrY output files..."
assert_file_exists "${OUTPUT_DIR}/chrY/chrY_male.vmiss" "chrY_male.vmiss exists"
assert_file_exists "${OUTPUT_DIR}/chrY/chrY_male.afreq" "chrY_male.afreq exists"

# --- Verify no plink2 errors in logs ---
echo ""
echo "  Checking plink2 logs for errors..."
CHRX_LOGS=("${OUTPUT_DIR}/chrX/chrX_all.log"
           "${OUTPUT_DIR}/chrX/chrX_female.log"
           "${OUTPUT_DIR}/chrX/chrX_male.log")
CHRY_LOGS=("${OUTPUT_DIR}/chrY/chrY_male.log")

for logf in "${CHRX_LOGS[@]}" "${CHRY_LOGS[@]}"; do
    if [[ -f "${logf}" ]]; then
        if grep -q '^Error:' "${logf}"; then
            echo "  FAIL: plink2 error in $(basename "${logf}"): $(grep '^Error:' "${logf}" | head -1)"
            (( FAIL++ )) || true
        else
            echo "  PASS: no plink2 errors in $(basename "${logf}")"
            (( PASS++ )) || true
        fi
    fi
done

# --- Verify --split-par was used (check chrX logs for PAR reclassification) ---
echo ""
echo "  Verifying --split-par reclassification in chrX logs..."
for logf in "${CHRX_LOGS[@]}"; do
    if [[ -f "${logf}" ]]; then
        if grep -q '\-\-split-par' "${logf}"; then
            echo "  PASS: --split-par present in $(basename "${logf}")"
            (( PASS++ )) || true
        else
            echo "  FAIL: --split-par not found in $(basename "${logf}")"
            (( FAIL++ )) || true
        fi
    fi
done

# --- Verify --update-sex actually matched samples ---
echo ""
echo "  Verifying --update-sex matched samples..."
for logf in "${CHRX_LOGS[@]}" "${CHRY_LOGS[@]}"; do
    if [[ -f "${logf}" ]]; then
        # Check for "N samples updated" where N > 0
        UPDATED=$(grep -oP '\d+ sample.? updated' "${logf}" | head -1 || true)
        if [[ -n "${UPDATED}" ]]; then
            N=$(echo "${UPDATED}" | grep -oP '^\d+')
            if [[ "${N}" -gt 0 ]]; then
                echo "  PASS: --update-sex matched ${N} samples in $(basename "${logf}")"
                (( PASS++ )) || true
            else
                echo "  FAIL: --update-sex matched 0 samples in $(basename "${logf}")"
                (( FAIL++ )) || true
            fi
        fi
    fi
done

# --- Verify --keep retained correct number of samples ---
echo ""
echo "  Verifying --keep retained samples..."
# chrX female log should show 2 females
if [[ -f "${OUTPUT_DIR}/chrX/chrX_female.log" ]]; then
    FEMALE_REMAINING=$(grep -oP '\d+ sample.? remaining' "${OUTPUT_DIR}/chrX/chrX_female.log" | head -1 || true)
    if [[ -n "${FEMALE_REMAINING}" ]]; then
        N=$(echo "${FEMALE_REMAINING}" | grep -oP '^\d+')
        assert_ge "${N}" "1" "chrX_female --keep retained samples"
    fi
fi
# chrY male log should show 2 males
if [[ -f "${OUTPUT_DIR}/chrY/chrY_male.log" ]]; then
    MALE_REMAINING=$(grep -oP '\d+ sample.? remaining' "${OUTPUT_DIR}/chrY/chrY_male.log" | head -1 || true)
    if [[ -n "${MALE_REMAINING}" ]]; then
        N=$(echo "${MALE_REMAINING}" | grep -oP '^\d+')
        assert_ge "${N}" "1" "chrY_male --keep retained samples"
    fi
fi

echo ""

# =================================================================
# Test 2: Verify #IID header format in generated files
# =================================================================
echo "--- Test 2: Verify #IID header format in generated files ---"

# Re-run pipeline with a fresh output dir to inspect tmp files before cleanup
OUTPUT_DIR2="${TMP_DIR}/sex_chr_qc2"

# We patch the script to skip cleanup so we can inspect tmp files.
# Use a broad pattern to match the rm -rf cleanup line robustly.
PATCHED_SCRIPT="${TMP_DIR}/patched_sex_qc.sh"
sed '/rm -rf.*TMP_DIR/s/^/# /' \
    "${REPO_DIR}/scripts/compute_sex_chr_variant_qc.sh" > "${PATCHED_SCRIPT}"
# copy utils.sh so that SCRIPT_DIR-relative source works from the temp dir
cp "${REPO_DIR}/scripts/utils.sh" "${TMP_DIR}/utils.sh"
chmod +x "${PATCHED_SCRIPT}"

bash "${PATCHED_SCRIPT}" \
    --vcf "${TMP_DIR}/mock.bcf" \
    --sample-qc "${TMP_DIR}/sample_qc.tsv" \
    --output-dir "${OUTPUT_DIR2}" \
    --genome CHM13 \
    --threads 1 > /dev/null 2>&1 || true

TMP2="${OUTPUT_DIR2}/tmp"

# sex_update.txt must have #IID header
if [[ -f "${TMP2}/sex_update.txt" ]]; then
    HEADER=$(head -1 "${TMP2}/sex_update.txt")
    assert_eq "${HEADER}" "#IID	SEX" "sex_update.txt has #IID<tab>SEX header"

    # Data lines should have 2 columns (IID, SEX)
    DATA_COLS=$(awk -F'\t' 'NR==2{print NF}' "${TMP2}/sex_update.txt")
    assert_eq "${DATA_COLS}" "2" "sex_update.txt data has 2 columns (IID, SEX)"

    # Males coded as 1, females as 2
    MALE_COUNT=$(awk -F'\t' 'NR>1 && $2==1' "${TMP2}/sex_update.txt" | wc -l | tr -d ' ')
    assert_eq "${MALE_COUNT}" "2" "sex_update.txt has 2 males"
    FEMALE_COUNT=$(awk -F'\t' 'NR>1 && $2==2' "${TMP2}/sex_update.txt" | wc -l | tr -d ' ')
    assert_eq "${FEMALE_COUNT}" "2" "sex_update.txt has 2 females"
else
    echo "  FAIL: sex_update.txt not found in tmp dir"
    (( FAIL++ )) || true
fi

# keep files must have #IID header
for keep_file in "${TMP2}/male_keep.txt" "${TMP2}/female_keep.txt"; do
    if [[ -f "${keep_file}" ]]; then
        HEADER=$(head -1 "${keep_file}")
        assert_eq "${HEADER}" "#IID" "$(basename "${keep_file}") has #IID header"

        # Data lines should have 1 column (IID only)
        DATA_COLS=$(awk 'NR==2{print NF}' "${keep_file}")
        assert_eq "${DATA_COLS}" "1" "$(basename "${keep_file}") data has 1 column"
    else
        echo "  FAIL: $(basename "${keep_file}") not found"
        (( FAIL++ )) || true
    fi
done

echo ""

# =================================================================
# Test 3: GRCh38 genome build
# =================================================================
echo "--- Test 3: Pipeline with GRCh38 genome build ---"

# Create a minimal VCF with GRCh38-style chrX PAR
cat > "${TMP_DIR}/mock_grch38.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chrX,length=156030895>
##contig=<ID=chrY,length=57217415>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SA	SB	SC
chrX	100000	grch38_xpar1	A	G	100	PASS	.	GT	0/1	0/1	0/0
chrX	50000000	grch38_xnonpar	T	C	100	PASS	.	GT	0/0	0/1	0/0
chrX	155800000	grch38_xpar2	G	A	100	PASS	.	GT	0/1	0/0	0/1
chrY	30000000	grch38_ynonpar	C	T	100	PASS	.	GT	0/0	./.	0/0
EOF

bcftools view -Ob -o "${TMP_DIR}/mock_grch38.bcf" "${TMP_DIR}/mock_grch38.vcf" 2>/dev/null
bcftools index "${TMP_DIR}/mock_grch38.bcf" 2>/dev/null

cat > "${TMP_DIR}/sample_qc_grch38.tsv" <<'EOF'
sample_id	call_rate	lrr_sd	computed_gender
SA	0.990	0.15	M
SB	0.985	0.18	F
SC	0.992	0.12	M
EOF

OUTPUT_DIR3="${TMP_DIR}/sex_chr_qc_grch38"
if bash "${REPO_DIR}/scripts/compute_sex_chr_variant_qc.sh" \
    --vcf "${TMP_DIR}/mock_grch38.bcf" \
    --sample-qc "${TMP_DIR}/sample_qc_grch38.tsv" \
    --output-dir "${OUTPUT_DIR3}" \
    --genome GRCh38 \
    --threads 1 > "${TMP_DIR}/grch38.log" 2>&1; then
    echo "  PASS: GRCh38 pipeline completed without error"
    (( PASS++ )) || true
else
    echo "  FAIL: GRCh38 pipeline exited with error"
    cat "${TMP_DIR}/grch38.log"
    (( FAIL++ )) || true
fi

# Check no plink2 errors
for logf in "${OUTPUT_DIR3}/chrX/chrX_all.log" \
            "${OUTPUT_DIR3}/chrX/chrX_female.log" \
            "${OUTPUT_DIR3}/chrX/chrX_male.log" \
            "${OUTPUT_DIR3}/chrY/chrY_male.log"; do
    if [[ -f "${logf}" ]]; then
        if grep -q '^Error:' "${logf}"; then
            echo "  FAIL: plink2 error in $(basename "${logf}"): $(grep '^Error:' "${logf}" | head -1)"
            (( FAIL++ )) || true
        else
            echo "  PASS: no plink2 errors in $(basename "${logf}")"
            (( PASS++ )) || true
        fi
    fi
done

echo ""

# =================================================================
# Test 4: Verify non-PAR chrX variants survive after --split-par
# =================================================================
echo "--- Test 4: Non-PAR chrX variants present in output ---"

# Check that chrX_all.afreq has non-PAR variants
if [[ -f "${OUTPUT_DIR}/chrX/chrX_all.afreq" ]]; then
    # Count data lines (skip header)
    N_VARIANTS=$(awk 'NR>1' "${OUTPUT_DIR}/chrX/chrX_all.afreq" | wc -l | tr -d ' ')
    # We expect at least the 3 non-PAR variants (xnonpar_v1, v2, v3)
    # PAR and XTR variants should be reclassified/excluded
    assert_ge "${N_VARIANTS}" "1" "chrX_all.afreq has non-PAR variant(s) (got ${N_VARIANTS})"
else
    echo "  FAIL: chrX_all.afreq not found"
    (( FAIL++ )) || true
fi

# Check chrY_male.afreq has non-PAR variants
if [[ -f "${OUTPUT_DIR}/chrY/chrY_male.afreq" ]]; then
    N_VARIANTS=$(awk 'NR>1' "${OUTPUT_DIR}/chrY/chrY_male.afreq" | wc -l | tr -d ' ')
    assert_ge "${N_VARIANTS}" "1" "chrY_male.afreq has non-PAR variant(s) (got ${N_VARIANTS})"
else
    echo "  FAIL: chrY_male.afreq not found"
    (( FAIL++ )) || true
fi

echo ""

# =================================================================
# Test 5: .hardy.x (non-PAR chrX) present, and collation parses HWE
# =================================================================
#
# Regression test for the bug where load_sex_chr_qc() only read
# chrX_female.hardy (PAR/XTR variants reclassified by --split-par)
# and never read chrX_female.hardy.x (the non-PAR chrX variants).
# The result was that chrX_female_hwe_p was NA for the vast majority
# of chrX variants in the collated output, which in turn made
# all_ancestries_chrX_female_hwe_pass silently fail for nearly
# every chrX variant.
#
# The fix: read both <prefix>.hardy and <prefix>.hardy.x and merge.
# The fix also normalises column-name handling so HET_A1_CT
# (biallelic autosomal), HET_A1_AX_CT (multiallelic autosomal), and
# FEMALE_HET_A1_CT/FEMALE_HET_A1_AX_CT (chrX-specific) are all
# recognised as the het-count column.
echo "--- Test 5: .hardy.x presence and HWE collation ---"

# The OUTPUT_DIR from Test 1 had xpar1, xpar2, xxtr and xnonpar_v1..3
# chrX variants.  After --split-par, PAR/XTR variants move to PAR1/PAR2
# (landing in chrX_female.hardy); only non-PAR chrX variants land in
# chrX_female.hardy.x.
HARDY_X="${OUTPUT_DIR}/chrX/chrX_female.hardy.x"

assert_file_nonempty "${HARDY_X}" "chrX_female.hardy.x exists and is non-empty"

if [[ -s "${HARDY_X}" ]]; then
    # Every data row (skip header) must have a numeric p-value in the
    # final column.  "nan" is acceptable only when flagged as such by
    # plink2 (e.g. monomorphic sites).
    N_X=$(awk 'NR>1' "${HARDY_X}" | wc -l | tr -d ' ')
    N_X_NUMERIC=$(awk 'NR>1 { if ($NF+0 == $NF) n++ } END { print n+0 }' "${HARDY_X}")
    assert_ge "${N_X}" "1" "chrX_female.hardy.x has non-PAR chrX variants"
    # At least one non-PAR variant must yield a parseable numeric p-value
    assert_ge "${N_X_NUMERIC}" "1" "chrX_female.hardy.x has numeric p-values"
fi

# Now run collate_variant_qc.py on the sex-chr QC dir and verify the
# resulting TSV reports non-NA HWE p-values for non-PAR chrX variants.
COLLATED="${TMP_DIR}/collated_variant_qc.tsv"
if python3 "${REPO_DIR}/scripts/collate_variant_qc.py" \
        --sex-chr-qc-dir "${OUTPUT_DIR}" \
        --output "${COLLATED}" > "${TMP_DIR}/collate.log" 2>&1; then
    echo "  PASS: collate_variant_qc.py succeeded on sex-chr QC dir"
    (( PASS++ )) || true
else
    echo "  FAIL: collate_variant_qc.py failed"
    cat "${TMP_DIR}/collate.log"
    (( FAIL++ )) || true
fi

if [[ -s "${COLLATED}" ]]; then
    # Locate the hwe / het columns in the data header (skip any
    # leading '#'-prefixed comment lines).
    HDR=$(awk '!/^#/ { print; exit }' "${COLLATED}")
    HWE_IDX=$(echo "${HDR}" | awk -F'\t' '{ for(i=1;i<=NF;i++) if($i=="chrX_female_hwe_p") { print i; exit } }')
    HET_IDX=$(echo "${HDR}" | awk -F'\t' '{ for(i=1;i<=NF;i++) if($i=="chrX_female_het_ct") { print i; exit } }')
    VID_IDX=1

    if [[ -z "${HWE_IDX}" || -z "${HET_IDX}" ]]; then
        echo "  FAIL: chrX_female_hwe_p / chrX_female_het_ct columns missing from collated output"
        (( FAIL++ )) || true
    else
        # Count rows for the non-PAR chrX variants (xnonpar_v1..v3)
        # where hwe_p is populated (not 'NA') and numeric.
        N_NONPAR_HWE=$(awk -F'\t' -v h="${HWE_IDX}" -v vid="${VID_IDX}" '
            !/^#/ && NR>1 && $vid ~ /^xnonpar/ {
                if ($h != "NA" && $h+0 == $h) n++
            } END { print n+0 }' "${COLLATED}")
        assert_ge "${N_NONPAR_HWE}" "1" \
            "collated TSV reports numeric chrX_female_hwe_p for non-PAR chrX variants"

        # Het count must also be populated (was always NA before the
        # HET_A1_CT column-name fix — the reader was looking for
        # HET_A1_AX_CT which plink2 biallelic output never emits).
        N_NONPAR_HET=$(awk -F'\t' -v h="${HET_IDX}" -v vid="${VID_IDX}" '
            !/^#/ && NR>1 && $vid ~ /^xnonpar/ {
                if ($h != "NA" && $h+0 == $h) n++
            } END { print n+0 }' "${COLLATED}")
        assert_ge "${N_NONPAR_HET}" "1" \
            "collated TSV reports numeric chrX_female_het_ct for non-PAR chrX variants"
    fi
fi

echo ""

# =================================================================
# Test 6: Ancestry-stratified sex-chr QC flags populate for non-PAR
# =================================================================
#
# Regression test for the --ancestry-sex-chr-qc pathway: since
# load_sex_chr_qc() is shared between the full-cohort and per-ancestry
# invocations, the fix must cascade so that
# all_ancestries_chrX_female_hwe_pass is 1/0 (never silently-NA)
# for non-PAR chrX variants when ancestry-stratified sex-chr QC
# is provided.
echo "--- Test 6: ancestry-stratified chrX HWE pass flag ---"

# Re-use OUTPUT_DIR from Test 1 as a single "EUR-only" ancestry.
COLLATED_ANC="${TMP_DIR}/collated_with_anc.tsv"
if python3 "${REPO_DIR}/scripts/collate_variant_qc.py" \
        --sex-chr-qc-dir "${OUTPUT_DIR}" \
        --ancestry-sex-chr-qc "EUR:${OUTPUT_DIR}" \
        --output "${COLLATED_ANC}" > "${TMP_DIR}/collate_anc.log" 2>&1; then
    echo "  PASS: collate_variant_qc.py succeeded with --ancestry-sex-chr-qc"
    (( PASS++ )) || true
else
    echo "  FAIL: collate_variant_qc.py failed"
    cat "${TMP_DIR}/collate_anc.log"
    (( FAIL++ )) || true
fi

if [[ -s "${COLLATED_ANC}" ]]; then
    HDR=$(awk '!/^#/ { print; exit }' "${COLLATED_ANC}")
    FLAG_IDX=$(echo "${HDR}" | awk -F'\t' '{ for(i=1;i<=NF;i++) if($i=="all_ancestries_chrX_female_hwe_pass") { print i; exit } }')
    EUR_HWE_IDX=$(echo "${HDR}" | awk -F'\t' '{ for(i=1;i<=NF;i++) if($i=="EUR_chrX_female_hwe_p") { print i; exit } }')

    if [[ -z "${FLAG_IDX}" ]]; then
        echo "  FAIL: all_ancestries_chrX_female_hwe_pass column missing"
        (( FAIL++ )) || true
    else
        # At least one non-PAR chrX variant must have a 1 or 0 flag
        # (not 'NA').  Before the fix, every non-PAR variant's
        # EUR_chrX_female_hwe_p was 'NA', which made the flag 0
        # or silently NA for all of them.
        N_NONPAR_FLAG=$(awk -F'\t' -v f="${FLAG_IDX}" '
            !/^#/ && NR>1 && $1 ~ /^xnonpar/ {
                if ($f == "0" || $f == "1") n++
            } END { print n+0 }' "${COLLATED_ANC}")
        assert_ge "${N_NONPAR_FLAG}" "1" \
            "all_ancestries_chrX_female_hwe_pass populated (0/1) for non-PAR chrX"
    fi

    if [[ -n "${EUR_HWE_IDX}" ]]; then
        N_NONPAR_EUR=$(awk -F'\t' -v h="${EUR_HWE_IDX}" '
            !/^#/ && NR>1 && $1 ~ /^xnonpar/ {
                if ($h != "NA" && $h+0 == $h) n++
            } END { print n+0 }' "${COLLATED_ANC}")
        assert_ge "${N_NONPAR_EUR}" "1" \
            "EUR_chrX_female_hwe_p populated for non-PAR chrX variants"
    fi
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
