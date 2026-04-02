#!/usr/bin/env bash
#
# compute_sex_chr_variant_qc.sh
#
# Compute variant-level QC for sex chromosomes (chrX non-PAR and chrY).
#
# Best practice for sex chromosome QC (Anderson et al. 2010; Marees et al.
# 2018) requires sex-aware handling:
#
#   chrX non-PAR:
#     - HWE: tested in females only (males are hemizygous, HWE N/A)
#     - Missingness & frequency: computed on ALL samples (both sexes)
#
#   chrY:
#     - All metrics computed on males only
#
# PAR/XTR regions are excluded so that the metrics reflect only the
# non-pseudoautosomal portions of the sex chromosomes.
#
# The script requires a sample QC file with a `computed_gender` column
# to determine sample sex for subsetting.  It also requires plink2 and
# bcftools.
#
# References:
#   - Anderson et al. 2010, Nat Protoc 5:1564-73
#   - Marees et al. 2018, Int J Methods Psychiatr Res 27:e1608
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/utils.sh"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Compute sex chromosome variant QC (chrX non-PAR females-only HWE, chrY males-only).

Required:
  --vcf FILE            Input VCF/BCF file
  --sample-qc FILE      Sample QC TSV with computed_gender column (M/F/U)
  --output-dir DIR      Output directory for sex chromosome QC results

Options:
  --genome STRING       Genome build for PAR/XTR exclusion: CHM13, GRCh38,
                        or GRCh37 (default: CHM13)
  --threads INT         Number of threads (default: all available)
  --help                Show this help message

Output files:
  chrx_variant_qc.vmiss     chrX non-PAR per-variant missingness (all samples)
  chrx_variant_qc.hardy     chrX non-PAR HWE test (females only)
  chrx_variant_qc.afreq     chrX non-PAR allele frequencies (all samples)
  chry_variant_qc.vmiss     chrY per-variant missingness (males only)
  chry_variant_qc.afreq     chrY allele frequencies (males only)
  sex_chr_variant_qc_summary.txt  Summary statistics
EOF
    exit 0
}

VCF=""
SAMPLE_QC=""
OUTPUT_DIR=""
GENOME="CHM13"
THREADS=$(nproc 2>/dev/null || echo 1)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)        VCF="$2"; shift 2 ;;
        --sample-qc)  SAMPLE_QC="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --genome)     GENOME="$2"; shift 2 ;;
        --threads)    THREADS="$2"; shift 2 ;;
        --help)       usage ;;
        *)            echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "${VCF}" ]]; then
    echo "Error: --vcf is required" >&2
    exit 1
fi
if [[ -z "${SAMPLE_QC}" ]]; then
    echo "Error: --sample-qc is required" >&2
    exit 1
fi
if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "Error: --output-dir is required" >&2
    exit 1
fi
if [[ ! -f "${VCF}" ]]; then
    echo "Error: VCF file not found: ${VCF}" >&2
    exit 1
fi
if [[ ! -f "${SAMPLE_QC}" ]]; then
    echo "Error: Sample QC file not found: ${SAMPLE_QC}" >&2
    exit 1
fi

# Check required tools
if ! command -v plink2 &>/dev/null; then
    echo "Warning: plink2 not found. Skipping sex chromosome variant QC." >&2
    exit 0
fi
if ! command -v bcftools &>/dev/null; then
    echo "Warning: bcftools not found. Skipping sex chromosome variant QC." >&2
    exit 0
fi

mkdir -p "${OUTPUT_DIR}"
TMP_DIR="${OUTPUT_DIR}/tmp/sex_chr_qc"
mkdir -p "${TMP_DIR}"
trap 'rm -rf "${TMP_DIR}"' EXIT

echo "Computing sex chromosome variant QC..."
echo "  VCF:     ${VCF}"
echo "  Genome:  ${GENOME}"
echo "  Threads: ${THREADS}"
SEX_QC_START=${SECONDS}

# -----------------------------------------------------------------
# Step 1: Extract male/female sample lists from the sample QC file
# -----------------------------------------------------------------
# Use header-aware column lookup for computed_gender
gender_col=$(head -1 "${SAMPLE_QC}" | tr '\t' '\n' \
    | grep -n '^computed_gender$' | cut -d: -f1)

if [[ -z "${gender_col}" ]]; then
    echo "Error: computed_gender column not found in ${SAMPLE_QC}" >&2
    exit 1
fi

# Extract female sample IDs (M/F/U classification)
FEMALES="${TMP_DIR}/females.txt"
MALES="${TMP_DIR}/males.txt"
awk -F'\t' -v col="${gender_col}" 'NR>1 && $col == "F" {print $1}' \
    "${SAMPLE_QC}" > "${FEMALES}"
awk -F'\t' -v col="${gender_col}" 'NR>1 && $col == "M" {print $1}' \
    "${SAMPLE_QC}" > "${MALES}"

N_FEMALES=$(wc -l < "${FEMALES}" | tr -d ' ')
N_MALES=$(wc -l < "${MALES}" | tr -d ' ')
echo "  Females: ${N_FEMALES}"
echo "  Males:   ${N_MALES}"

if [[ "${N_FEMALES}" -eq 0 && "${N_MALES}" -eq 0 ]]; then
    echo "Warning: No male or female samples found. Skipping sex chromosome QC." >&2
    exit 0
fi

# -----------------------------------------------------------------
# Step 2: Generate PAR/XTR exclusion BED
# -----------------------------------------------------------------
PAR_BED="${TMP_DIR}/par_xtr.bed"
get_par_xtr_bed "${GENOME}" "${PAR_BED}"

# Determine chrX/chrY contig names from the VCF
CHRX_NAME=""
CHRY_NAME=""
for candidate in chrX X; do
    if bcftools index -s "${VCF}" 2>/dev/null | cut -f1 | grep -qx "${candidate}"; then
        CHRX_NAME="${candidate}"
        break
    fi
done
for candidate in chrY Y; do
    if bcftools index -s "${VCF}" 2>/dev/null | cut -f1 | grep -qx "${candidate}"; then
        CHRY_NAME="${candidate}"
        break
    fi
done

echo "  chrX contig: ${CHRX_NAME:-not found}"
echo "  chrY contig: ${CHRY_NAME:-not found}"

# -----------------------------------------------------------------
# Step 3: Filter VCF to remove intensity-only probes
# -----------------------------------------------------------------
FILTERED_VCF="${TMP_DIR}/filtered.bcf"
bcftools view -e 'INFO/INTENSITY_ONLY=1' --threads "${THREADS}" \
    -Ob -o "${FILTERED_VCF}" "${VCF}" 2>/dev/null
bcftools index --threads "${THREADS}" "${FILTERED_VCF}" 2>/dev/null

# -----------------------------------------------------------------
# Step 4: chrX non-PAR variant QC
# -----------------------------------------------------------------
if [[ -n "${CHRX_NAME}" ]]; then
    echo ""
    echo "  --- chrX non-PAR variant QC ---"

    CHRX_PREFIX="${OUTPUT_DIR}/chrx_variant_qc"

    # 4a. Missingness + frequency on ALL samples (chrX non-PAR)
    echo "    Missingness & frequency (all samples)..."
    plink2 \
        --bcf "${FILTERED_VCF}" \
        --chr "${CHRX_NAME}" \
        --exclude-range "${PAR_BED}" \
        --allow-extra-chr \
        --threads "${THREADS}" \
        --missing variant-only \
        --freq \
        --out "${CHRX_PREFIX}" 2>&1 | tail -3 | sed 's/^/    /'

    # 4b. HWE on females only (chrX non-PAR)
    # Males are hemizygous — HWE is only meaningful for diploid (female) samples.
    if [[ "${N_FEMALES}" -gt 0 ]]; then
        echo "    HWE test (females only, n=${N_FEMALES})..."
        CHRX_HWE_PREFIX="${TMP_DIR}/chrx_hwe_females"
        plink2 \
            --bcf "${FILTERED_VCF}" \
            --chr "${CHRX_NAME}" \
            --exclude-range "${PAR_BED}" \
            --keep "${FEMALES}" \
            --allow-extra-chr \
            --threads "${THREADS}" \
            --hardy midp \
            --out "${CHRX_HWE_PREFIX}" 2>&1 | tail -3 | sed 's/^/    /'

        # Copy the .hardy file to the output with the standard name
        if [[ -f "${CHRX_HWE_PREFIX}.hardy" ]]; then
            cp "${CHRX_HWE_PREFIX}.hardy" "${CHRX_PREFIX}.hardy"
        fi
    else
        echo "    Skipping HWE (no female samples)."
    fi

    N_CHRX=$(awk 'NR>1' "${CHRX_PREFIX}.vmiss" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
    echo "    chrX non-PAR variants: ${N_CHRX}"
else
    echo "  Skipping chrX QC: chrX contig not found in VCF."
fi

# -----------------------------------------------------------------
# Step 5: chrY variant QC (males only)
# -----------------------------------------------------------------
if [[ -n "${CHRY_NAME}" && "${N_MALES}" -gt 0 ]]; then
    echo ""
    echo "  --- chrY variant QC (males only, n=${N_MALES}) ---"

    CHRY_PREFIX="${OUTPUT_DIR}/chry_variant_qc"

    plink2 \
        --bcf "${FILTERED_VCF}" \
        --chr "${CHRY_NAME}" \
        --exclude-range "${PAR_BED}" \
        --keep "${MALES}" \
        --allow-extra-chr \
        --threads "${THREADS}" \
        --missing variant-only \
        --freq \
        --out "${CHRY_PREFIX}" 2>&1 | tail -3 | sed 's/^/    /'

    N_CHRY=$(awk 'NR>1' "${CHRY_PREFIX}.vmiss" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
    echo "    chrY variants: ${N_CHRY}"
else
    if [[ -z "${CHRY_NAME}" ]]; then
        echo "  Skipping chrY QC: chrY contig not found in VCF."
    else
        echo "  Skipping chrY QC: no male samples."
    fi
fi

# -----------------------------------------------------------------
# Step 6: Summary
# -----------------------------------------------------------------
SUMMARY="${OUTPUT_DIR}/sex_chr_variant_qc_summary.txt"
{
    echo "======================================"
    echo "  Sex Chromosome Variant QC Summary"
    echo "======================================"
    echo ""
    echo "Genome build: ${GENOME}"
    echo "Females: ${N_FEMALES}"
    echo "Males:   ${N_MALES}"
    echo ""

    if [[ -n "${CHRX_NAME}" && -f "${OUTPUT_DIR}/chrx_variant_qc.vmiss" ]]; then
        n_chrx=$(awk 'NR>1' "${OUTPUT_DIR}/chrx_variant_qc.vmiss" | wc -l | tr -d ' ')
        echo "chrX non-PAR variants: ${n_chrx}"
        echo "  Missingness & frequency: all samples"
        echo "  HWE: females only (n=${N_FEMALES})"
        if [[ -f "${OUTPUT_DIR}/chrx_variant_qc.vmiss" ]]; then
            awk 'NR>1 {
                n++; miss += $NF
                if ($NF+0 > 0.05) n_05++
            } END {
                if (n>0) {
                    printf "  Mean missingness: %.6f\n", miss/n
                    printf "  Variants with missingness >5%%: %d (%.2f%%)\n", n_05+0, (n_05+0)*100/n
                }
            }' "${OUTPUT_DIR}/chrx_variant_qc.vmiss"
        fi
    fi

    echo ""

    if [[ -n "${CHRY_NAME}" && -f "${OUTPUT_DIR}/chry_variant_qc.vmiss" ]]; then
        n_chry=$(awk 'NR>1' "${OUTPUT_DIR}/chry_variant_qc.vmiss" | wc -l | tr -d ' ')
        echo "chrY variants: ${n_chry}"
        echo "  All metrics: males only (n=${N_MALES})"
        if [[ -f "${OUTPUT_DIR}/chry_variant_qc.vmiss" ]]; then
            awk 'NR>1 {
                n++; miss += $NF
                if ($NF+0 > 0.05) n_05++
            } END {
                if (n>0) {
                    printf "  Mean missingness: %.6f\n", miss/n
                    printf "  Variants with missingness >5%%: %d (%.2f%%)\n", n_05+0, (n_05+0)*100/n
                }
            }' "${OUTPUT_DIR}/chry_variant_qc.vmiss"
        fi
    fi

    echo ""
    echo "======================================"
} | tee "${SUMMARY}"

SEX_QC_ELAPSED=$(( SECONDS - SEX_QC_START ))
echo ""
echo "Sex chromosome variant QC files:"
for f in chrx_variant_qc.vmiss chrx_variant_qc.hardy chrx_variant_qc.afreq \
         chry_variant_qc.vmiss chry_variant_qc.afreq; do
    [[ -f "${OUTPUT_DIR}/${f}" ]] && echo "  ${OUTPUT_DIR}/${f}"
done
echo "  ${SUMMARY}"
echo "  [diag] Total sex chr QC time: ${SEX_QC_ELAPSED}s"
