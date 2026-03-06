#!/usr/bin/env bash
#
# compute_variant_qc.sh
#
# Compute variant-level QC metrics using PLINK2 on autosomal variants.
# Produces per-variant missingness, Hardy-Weinberg equilibrium tests,
# and allele frequency statistics. These are standard best-practice
# metrics for genotyping array quality control.
#
# Can be run standalone or invoked by stage1/stage2 scripts.
#
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Compute variant-level QC metrics using PLINK2 on autosomal variants.

Required:
  --vcf FILE          Input VCF/BCF file
  --output-dir DIR    Output directory for QC results

Options:
  --prefix STRING     Output file prefix (default: variant_qc)
  --threads INT       Number of threads (default: all available)
  --help              Show this help message

Output files:
  <prefix>.vmiss      Per-variant missingness
  <prefix>.hardy      Hardy-Weinberg equilibrium test
  <prefix>.afreq      Allele frequency statistics
  <prefix>.log        PLINK2 log
  variant_qc_summary.txt   Summary statistics
EOF
    exit 0
}

VCF=""
OUTPUT_DIR=""
PREFIX_NAME="variant_qc"
THREADS=$(nproc 2>/dev/null || echo 1)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)        VCF="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --prefix)     PREFIX_NAME="$2"; shift 2 ;;
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
if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "Error: --output-dir is required" >&2
    exit 1
fi
if [[ ! -f "${VCF}" ]]; then
    echo "Error: VCF file not found: ${VCF}" >&2
    exit 1
fi

# Check plink2 is available
if ! command -v plink2 &>/dev/null; then
    echo "Warning: plink2 not found. Skipping variant-level QC." >&2
    echo "  Install plink2 from https://www.cog-genomics.org/plink/2.0/" >&2
    exit 0
fi

mkdir -p "${OUTPUT_DIR}"
TMP_DIR="${OUTPUT_DIR}/tmp/compute_variant_qc"
mkdir -p "${TMP_DIR}"

PREFIX="${OUTPUT_DIR}/${PREFIX_NAME}"

echo "Computing variant QC with PLINK2..."
echo "  VCF:     ${VCF}"
echo "  Output:  ${OUTPUT_DIR}"
echo "  Threads: ${THREADS}"
VAR_QC_START=${SECONDS}

# Filter out intensity-only probes before passing to plink2
# These probes lack genotype calls and would inflate missingness
FILTERED_VCF=$(mktemp -p "${TMP_DIR}" --suffix=.bcf filtered_variant_qc.XXXXXX)
trap 'rm -f "${FILTERED_VCF}" "${FILTERED_VCF}.csi"' EXIT

bcftools view -e 'INFO/INTENSITY_ONLY=1' --threads "${THREADS}" \
    -Ob -o "${FILTERED_VCF}" "${VCF}" 2>/dev/null
bcftools index --threads "${THREADS}" "${FILTERED_VCF}" 2>/dev/null
N_FILTERED=$(bcftools index -n "${FILTERED_VCF}" 2>/dev/null || echo "NA")
echo "  [diag] Filtered autosomal-input BCF ready: ${FILTERED_VCF} (variants=${N_FILTERED})"

# Run plink2 to compute variant-level QC on autosomes
# --missing variant-only: per-variant missingness
# --hardy midp: HWE test with mid-p adjustment (recommended for small samples)
# --freq: allele frequency distribution
# --het: per-sample inbreeding coefficient (F statistic)
plink2 \
    --bcf "${FILTERED_VCF}" \
    --autosome \
    --allow-extra-chr \
    --threads "${THREADS}" \
    --missing variant-only \
    --hardy midp \
    --freq \
    --het \
    --out "${PREFIX}" 2>&1 | tail -5
echo "  [diag] PLINK2 variant QC metrics generated for prefix: ${PREFIX}"

# Compute Ti/Tv ratio via bcftools stats (lightweight, single pass)
echo "Computing Ti/Tv ratio..."
TSTV_FILE="${OUTPUT_DIR}/tstv_stats.txt"
bcftools stats --threads "${THREADS}" "${FILTERED_VCF}" 2>/dev/null | \
    awk '/^SN/ || /^TSTV/' > "${TSTV_FILE}" 2>/dev/null || true
echo "  [diag] Ti/Tv stats written: ${TSTV_FILE}"

# Generate summary statistics
SUMMARY="${OUTPUT_DIR}/variant_qc_summary.txt"
{
    echo "======================================"
    echo "  Variant QC Summary"
    echo "======================================"
    echo ""

    # Variant count
    if [[ -f "${PREFIX}.vmiss" ]]; then
        n_variants=$(awk 'NR>1' "${PREFIX}.vmiss" | wc -l | tr -d ' ')
        echo "Total autosomal variants: ${n_variants}"
        echo ""

        # Missingness summary
        echo "--- Variant Missingness ---"
        awk 'NR>1 {
            n++; miss += $NF
            if ($NF+0 > 0.05) n_05++
            if ($NF+0 > 0.02) n_02++
            if ($NF+0 > 0.01) n_01++
            if (n==1 || $NF+0 < min_m) min_m=$NF+0
            if (n==1 || $NF+0 > max_m) max_m=$NF+0
        } END {
            if (n>0) {
                printf "Variants with missingness >5%%: %d (%.2f%%)\n", n_05+0, (n_05+0)*100/n
                printf "Variants with missingness >2%%: %d (%.2f%%)\n", n_02+0, (n_02+0)*100/n
                printf "Variants with missingness >1%%: %d (%.2f%%)\n", n_01+0, (n_01+0)*100/n
                printf "Mean variant missingness: %.6f\n", miss/n
                printf "Min variant missingness:  %.6f\n", min_m
                printf "Max variant missingness:  %.6f\n", max_m
            }
        }' "${PREFIX}.vmiss"
    fi

    echo ""

    # HWE summary
    if [[ -f "${PREFIX}.hardy" ]]; then
        echo "--- Hardy-Weinberg Equilibrium ---"
        awk 'NR>1 {
            n++; p = $NF
            if (p+0 < 1e-6) n_6++
            if (p+0 < 1e-10) n_10++
            if (p+0 < 1e-20) n_20++
        } END {
            if (n>0) {
                printf "Variants failing HWE (p<1e-6):  %d (%.2f%%)\n", n_6+0, (n_6+0)*100/n
                printf "Variants failing HWE (p<1e-10): %d (%.2f%%)\n", n_10+0, (n_10+0)*100/n
                printf "Variants failing HWE (p<1e-20): %d (%.2f%%)\n", n_20+0, (n_20+0)*100/n
            }
        }' "${PREFIX}.hardy"
    fi

    echo ""

    # MAF summary
    if [[ -f "${PREFIX}.afreq" ]]; then
        echo "--- Allele Frequency Distribution ---"
        awk 'NR>1 {
            n++
            af = $5+0
            maf = (af > 0.5) ? 1-af : af
            maf_sum += maf
            if (maf < 0.001) n_mono++
            if (maf < 0.01) n_rare++
            if (maf >= 0.01 && maf < 0.05) n_low++
            if (maf >= 0.05) n_common++
        } END {
            if (n>0) {
                printf "Monomorphic (MAF<0.1%%):     %d (%.2f%%)\n", n_mono+0, (n_mono+0)*100/n
                printf "Rare (MAF<1%%):              %d (%.2f%%)\n", n_rare+0, (n_rare+0)*100/n
                printf "Low frequency (1-5%%):       %d (%.2f%%)\n", n_low+0, (n_low+0)*100/n
                printf "Common (MAF>=5%%):           %d (%.2f%%)\n", n_common+0, (n_common+0)*100/n
                printf "Mean MAF: %.4f\n", maf_sum/n
            }
        }' "${PREFIX}.afreq"
    fi

    echo ""

    # Heterozygosity / inbreeding coefficient
    if [[ -f "${PREFIX}.het" ]]; then
        echo "--- Per-Sample Heterozygosity (Inbreeding Coefficient F) ---"
        awk 'NR>1 {
            n++; f = $NF+0; fsum += f
            if (n==1 || f < min_f) min_f = f
            if (n==1 || f > max_f) max_f = f
            if (f > 0.05 || f < -0.05) n_outlier++
        } END {
            if (n>0) {
                printf "Samples:                    %d\n", n
                printf "Mean F (inbreeding coeff):  %.4f\n", fsum/n
                printf "Min F:                      %.4f\n", min_f
                printf "Max F:                      %.4f\n", max_f
                printf "Outliers (|F|>0.05):        %d (%.1f%%)\n", n_outlier+0, (n_outlier+0)*100/n
            }
        }' "${PREFIX}.het"
    fi

    echo ""

    # Ti/Tv ratio
    if [[ -f "${TSTV_FILE}" ]]; then
        echo "--- Transition/Transversion Ratio ---"
        awk -F'\t' '/^TSTV/ {
            printf "Transitions:    %s\n", $3
            printf "Transversions:  %s\n", $4
            printf "Ti/Tv ratio:    %s\n", $5
        }' "${TSTV_FILE}" 2>/dev/null || echo "  (Could not parse Ti/Tv ratio)"
    fi

    echo ""
    echo "======================================"
} | tee "${SUMMARY}"

echo ""
echo "Variant QC files:"
for ext in vmiss hardy afreq het log; do
    [[ -f "${PREFIX}.${ext}" ]] && echo "  ${PREFIX}.${ext}"
done
echo "  ${SUMMARY}"
if [[ -f "${TSTV_FILE}" ]]; then
    echo "  ${TSTV_FILE}"
fi
VAR_QC_ELAPSED=$(( SECONDS - VAR_QC_START ))
echo "  [diag] Total variant QC stage time: ${VAR_QC_ELAPSED}s"
