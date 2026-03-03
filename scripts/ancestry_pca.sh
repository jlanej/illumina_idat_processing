#!/usr/bin/env bash
#
# ancestry_pca.sh
#
# Perform stringent best-practice variant and sample QC, compute ancestry
# principal components using flashpca2 on QC'd data, and project PCs to
# all samples.
#
# Steps:
#   1. Stringent variant QC (plink2): missingness, HWE, MAF
#   2. Stringent sample QC (plink2): call rate
#   3. LD pruning (plink2)
#   4. PCA on QC'd variants and samples (flashpca2)
#   5. Project PCs to all samples (flashpca2 --project)
#
# Usage:
#   ./ancestry_pca.sh --vcf input.bcf --output-dir pca_output
#
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Perform stringent QC and compute ancestry PCs using flashpca2.

Required:
  --vcf FILE          Input VCF/BCF file
  --output-dir DIR    Output directory

Options:
  --n-pcs INT         Number of principal components (default: 20)
  --max-missing FLOAT Max variant missingness rate (default: 0.02)
  --hwe-p FLOAT       HWE p-value threshold (default: 1e-6)
  --min-maf FLOAT     Minimum minor allele frequency (default: 0.05)
  --min-call-rate FLOAT Min sample call rate (default: 0.98)
  --ld-window INT     LD pruning window in kb (default: 1000)
  --ld-step INT       LD pruning step size (default: 100)
  --ld-r2 FLOAT       LD pruning r-squared threshold (default: 0.1)
  --threads INT       Number of threads (default: all available)
  --help              Show this help message

Output files:
  pca_projections.tsv      All-sample PC projections (default 20 PCs)
  qc_summary.txt           QC filtering summary
  pca_eigenvalues.txt      PCA eigenvalues
  pca_qc_samples.txt       QC'd sample set used for PCA
EOF
    exit 0
}

VCF=""
OUTPUT_DIR=""
N_PCS=20
MAX_MISSING=0.02
HWE_P="1e-6"
MIN_MAF=0.05
MIN_CALL_RATE=0.98
LD_WINDOW=1000
LD_STEP=100
LD_R2=0.1
THREADS=$(nproc 2>/dev/null || echo 1)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)            VCF="$2"; shift 2 ;;
        --output-dir)     OUTPUT_DIR="$2"; shift 2 ;;
        --n-pcs)          N_PCS="$2"; shift 2 ;;
        --max-missing)    MAX_MISSING="$2"; shift 2 ;;
        --hwe-p)          HWE_P="$2"; shift 2 ;;
        --min-maf)        MIN_MAF="$2"; shift 2 ;;
        --min-call-rate)  MIN_CALL_RATE="$2"; shift 2 ;;
        --ld-window)      LD_WINDOW="$2"; shift 2 ;;
        --ld-step)        LD_STEP="$2"; shift 2 ;;
        --ld-r2)          LD_R2="$2"; shift 2 ;;
        --threads)        THREADS="$2"; shift 2 ;;
        --help)           usage ;;
        *)                echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

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

# Verify tools
if ! command -v plink2 &>/dev/null; then
    echo "Error: plink2 not found." >&2
    exit 1
fi
if ! command -v flashpca &>/dev/null; then
    echo "Warning: flashpca not found. PCA computation will be skipped." >&2
    echo "  Install flashpca2 from https://github.com/gabraham/flashpca" >&2
    HAVE_FLASHPCA="false"
else
    HAVE_FLASHPCA="true"
fi

mkdir -p "${OUTPUT_DIR}"

PREFIX="${OUTPUT_DIR}/ancestry_pca"
SUMMARY="${OUTPUT_DIR}/qc_summary.txt"

echo "======================================================"
echo "  Ancestry PCA: Stringent QC and PC Computation"
echo "======================================================"
echo ""
echo "  VCF:            ${VCF}"
echo "  Output:         ${OUTPUT_DIR}"
echo "  N PCs:          ${N_PCS}"
echo "  Max missing:    ${MAX_MISSING}"
echo "  HWE p:          ${HWE_P}"
echo "  Min MAF:        ${MIN_MAF}"
echo "  Min call rate:  ${MIN_CALL_RATE}"
echo "  LD pruning:     window=${LD_WINDOW}kb step=${LD_STEP} r2=${LD_R2}"
echo "  Threads:        ${THREADS}"
echo ""

# -----------------------------------------------------------------
# Step 1: Filter intensity-only probes
# -----------------------------------------------------------------
echo "Step 1: Filtering intensity-only probes..."
FILTERED_VCF=$(mktemp --suffix=.bcf)
trap 'rm -f "${FILTERED_VCF}" "${FILTERED_VCF}.csi"' EXIT

bcftools view -e 'INFO/INTENSITY_ONLY=1' --threads "${THREADS}" \
    -Ob -o "${FILTERED_VCF}" "${VCF}" 2>/dev/null
bcftools index "${FILTERED_VCF}" 2>/dev/null

# -----------------------------------------------------------------
# Step 2: Stringent variant and sample QC with plink2
# -----------------------------------------------------------------
echo "Step 2: Applying stringent variant and sample QC..."

# First pass: variant QC
plink2 \
    --bcf "${FILTERED_VCF}" \
    --autosome \
    --allow-extra-chr \
    --threads "${THREADS}" \
    --geno "${MAX_MISSING}" \
    --hwe "${HWE_P}" \
    --maf "${MIN_MAF}" \
    --mind "$(awk -v cr="${MIN_CALL_RATE}" 'BEGIN{printf "%.4f", 1-cr}')" \
    --make-bed \
    --out "${PREFIX}_qc" 2>&1 | tail -10

N_VARIANTS_QC=0
N_SAMPLES_QC=0
if [[ -f "${PREFIX}_qc.bim" ]]; then
    N_VARIANTS_QC=$(wc -l < "${PREFIX}_qc.bim" | tr -d ' ')
fi
if [[ -f "${PREFIX}_qc.fam" ]]; then
    N_SAMPLES_QC=$(wc -l < "${PREFIX}_qc.fam" | tr -d ' ')
fi

echo "  QC'd variants: ${N_VARIANTS_QC}"
echo "  QC'd samples:  ${N_SAMPLES_QC}"

# Save QC'd sample list
if [[ -f "${PREFIX}_qc.fam" ]]; then
    awk '{print $2}' "${PREFIX}_qc.fam" > "${OUTPUT_DIR}/pca_qc_samples.txt"
fi

# -----------------------------------------------------------------
# Step 3: LD pruning
# -----------------------------------------------------------------
echo "Step 3: LD pruning..."

plink2 \
    --bfile "${PREFIX}_qc" \
    --threads "${THREADS}" \
    --indep-pairwise "${LD_WINDOW}kb" "${LD_STEP}" "${LD_R2}" \
    --out "${PREFIX}_ldprune" 2>&1 | tail -5

N_LD_PRUNED=0
if [[ -f "${PREFIX}_ldprune.prune.in" ]]; then
    N_LD_PRUNED=$(wc -l < "${PREFIX}_ldprune.prune.in" | tr -d ' ')
fi
echo "  Variants after LD pruning: ${N_LD_PRUNED}"

# Create LD-pruned dataset
plink2 \
    --bfile "${PREFIX}_qc" \
    --extract "${PREFIX}_ldprune.prune.in" \
    --threads "${THREADS}" \
    --make-bed \
    --out "${PREFIX}_ldpruned" 2>&1 | tail -3

# -----------------------------------------------------------------
# Step 4: PCA computation with flashpca2
# -----------------------------------------------------------------
if [[ "${HAVE_FLASHPCA}" == "true" && -f "${PREFIX}_ldpruned.bed" ]]; then
    echo "Step 4: Computing ${N_PCS} ancestry PCs with flashpca2..."

    flashpca \
        --bfile "${PREFIX}_ldpruned" \
        --ndim "${N_PCS}" \
        --outpc "${PREFIX}_pcs.txt" \
        --outval "${OUTPUT_DIR}/pca_eigenvalues.txt" \
        --outvec "${PREFIX}_loadings.txt" \
        --outmeansd "${PREFIX}_meansd.txt" \
        --numthreads "${THREADS}" 2>&1 | tail -5

    echo "  PCA complete."

    # -----------------------------------------------------------------
    # Step 5: Project PCs to all samples
    # -----------------------------------------------------------------
    echo "Step 5: Projecting PCs to all samples..."

    # Create all-sample LD-pruned dataset (no sample QC filter)
    plink2 \
        --bcf "${FILTERED_VCF}" \
        --autosome \
        --allow-extra-chr \
        --extract "${PREFIX}_ldprune.prune.in" \
        --threads "${THREADS}" \
        --make-bed \
        --out "${PREFIX}_all_ldpruned" 2>&1 | tail -3

    if [[ -f "${PREFIX}_all_ldpruned.bed" ]]; then
        flashpca \
            --bfile "${PREFIX}_all_ldpruned" \
            --project \
            --inmeansd "${PREFIX}_meansd.txt" \
            --inload "${PREFIX}_loadings.txt" \
            --outproj "${PREFIX}_projections.txt" \
            --numthreads "${THREADS}" 2>&1 | tail -5

        # Format projections into a clean TSV
        if [[ -f "${PREFIX}_projections.txt" ]]; then
            {
                printf "FID\tIID"
                for i in $(seq 1 "${N_PCS}"); do printf "\tPC%d" "$i"; done
                printf "\n"
                # Normalize whitespace to tabs for consistent TSV output
                sed 's/[[:space:]]\+/\t/g' "${PREFIX}_projections.txt"
            } > "${OUTPUT_DIR}/pca_projections.tsv"
            echo "  Projections saved: ${OUTPUT_DIR}/pca_projections.tsv"
        fi
    else
        echo "  Warning: Could not create all-sample dataset for projection."
    fi
else
    echo "Step 4-5: Skipping PCA (flashpca not available or no QC'd data)."

    # Fallback: use plink2 --pca for PCA computation
    if [[ -f "${PREFIX}_ldpruned.bed" ]]; then
        echo "  Using plink2 --pca as fallback..."
        plink2 \
            --bfile "${PREFIX}_ldpruned" \
            --pca "${N_PCS}" \
            --threads "${THREADS}" \
            --out "${PREFIX}_plink_pca" 2>&1 | tail -5

        if [[ -f "${PREFIX}_plink_pca.eigenvec" ]]; then
            # plink2 eigenvec already has header; rename to standard output
            cp "${PREFIX}_plink_pca.eigenvec" "${OUTPUT_DIR}/pca_projections.tsv"
            echo "  PCA projections (plink2): ${OUTPUT_DIR}/pca_projections.tsv"
        fi
        if [[ -f "${PREFIX}_plink_pca.eigenval" ]]; then
            cp "${PREFIX}_plink_pca.eigenval" "${OUTPUT_DIR}/pca_eigenvalues.txt"
        fi
    fi
fi

# -----------------------------------------------------------------
# QC Summary
# -----------------------------------------------------------------
{
    echo "======================================================"
    echo "  Ancestry PCA QC Summary"
    echo "======================================================"
    echo ""
    echo "Input VCF: ${VCF}"
    echo ""
    echo "QC Parameters:"
    echo "  Max variant missingness: ${MAX_MISSING}"
    echo "  HWE p-value threshold:  ${HWE_P}"
    echo "  Min MAF:                ${MIN_MAF}"
    echo "  Min sample call rate:   ${MIN_CALL_RATE}"
    echo "  LD pruning:             window=${LD_WINDOW}kb step=${LD_STEP} r2=${LD_R2}"
    echo ""
    echo "Results:"
    echo "  Variants after QC:      ${N_VARIANTS_QC}"
    echo "  Samples after QC:       ${N_SAMPLES_QC}"
    echo "  Variants after LD prune: ${N_LD_PRUNED}"
    echo "  PCs computed:           ${N_PCS}"
    if [[ "${HAVE_FLASHPCA}" == "true" ]]; then
        echo "  PCA method:             flashpca2"
    else
        echo "  PCA method:             plink2 --pca (flashpca2 not available)"
    fi
    echo ""
    echo "Output files:"
    for f in "${OUTPUT_DIR}/pca_projections.tsv" "${OUTPUT_DIR}/pca_eigenvalues.txt" \
             "${OUTPUT_DIR}/pca_qc_samples.txt" "${OUTPUT_DIR}/qc_summary.txt"; do
        [[ -f "${f}" ]] && echo "  ${f}"
    done
    echo ""
    echo "======================================================"
} | tee "${SUMMARY}"

echo ""
echo "Ancestry PCA complete."
