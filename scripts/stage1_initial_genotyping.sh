#!/usr/bin/env bash
#
# stage1_initial_genotyping.sh
#
# Stage 1: Initial genotyping of Illumina IDAT data.
#
# Converts IDAT files to GTC (genotype call) files using the idat2gtc
# bcftools plugin (equivalent to Illumina GenCall), then converts GTC
# files to VCF format with BAF and LRR intensities. Computes per-sample
# QC metrics including call rate and LRR standard deviation.
#
# This produces a raw estimate that is used to identify high-quality
# samples for reclustering in Stage 2.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Defaults
THREADS=1
BATCH_SIZE=256

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Stage 1: Initial genotyping - convert IDAT to VCF with QC metrics.

Required:
  --idat-dir DIR         Directory containing IDAT files
  --bpm FILE             Illumina BPM manifest file
  --egt FILE             Illumina EGT cluster file
  --csv FILE             Illumina CSV manifest file
  --ref-fasta FILE       Reference genome FASTA file
  --output-dir DIR       Output directory

Options:
  --sample-sheet FILE    Tab-delimited file mapping sample_id to IDAT files
                         (auto-generated from IDAT dir if not provided)
  --threads INT          Number of threads (default: ${THREADS})
  --batch-size INT       Batch size for IDAT to GTC conversion (default: ${BATCH_SIZE})
  --help                 Show this help message
EOF
    exit 0
}

IDAT_DIR=""
BPM=""
EGT=""
CSV=""
REF_FASTA=""
OUTPUT_DIR=""
SAMPLE_SHEET=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --idat-dir)      IDAT_DIR="$2"; shift 2 ;;
        --bpm)           BPM="$2"; shift 2 ;;
        --egt)           EGT="$2"; shift 2 ;;
        --csv)           CSV="$2"; shift 2 ;;
        --ref-fasta)     REF_FASTA="$2"; shift 2 ;;
        --output-dir)    OUTPUT_DIR="$2"; shift 2 ;;
        --sample-sheet)  SAMPLE_SHEET="$2"; shift 2 ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --batch-size)    BATCH_SIZE="$2"; shift 2 ;;
        --help)          usage ;;
        *)               echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Validate required arguments
for var in IDAT_DIR BPM EGT CSV REF_FASTA OUTPUT_DIR; do
    if [[ -z "${!var}" ]]; then
        echo "Error: --$(echo ${var} | tr '[:upper:]' '[:lower:]' | tr '_' '-') is required" >&2
        exit 1
    fi
done

# Validate files exist
for f in "${BPM}" "${EGT}" "${CSV}" "${REF_FASTA}"; do
    if [[ ! -f "${f}" ]]; then
        echo "Error: File not found: ${f}" >&2
        exit 1
    fi
done

if [[ ! -d "${IDAT_DIR}" ]]; then
    echo "Error: IDAT directory not found: ${IDAT_DIR}" >&2
    exit 1
fi

GTC_DIR="${OUTPUT_DIR}/gtc"
VCF_DIR="${OUTPUT_DIR}/vcf"
QC_DIR="${OUTPUT_DIR}/qc"
mkdir -p "${GTC_DIR}" "${VCF_DIR}" "${QC_DIR}"

echo "============================================"
echo "Stage 1: Initial Genotyping"
echo "============================================"
echo "IDAT dir:   ${IDAT_DIR}"
echo "BPM:        ${BPM}"
echo "EGT:        ${EGT}"
echo "CSV:        ${CSV}"
echo "Reference:  ${REF_FASTA}"
echo "Output dir: ${OUTPUT_DIR}"
echo "Threads:    ${THREADS}"
echo "============================================"
echo ""

# ---------------------------------------------------------------
# Step 1: Generate sample sheet if not provided
# ---------------------------------------------------------------
if [[ -z "${SAMPLE_SHEET}" ]]; then
    SAMPLE_SHEET="${OUTPUT_DIR}/sample_sheet.tsv"
    echo "Generating sample sheet from IDAT files..."
    {
        echo -e "sample_id\tgreen_idat\tred_idat"
        find "${IDAT_DIR}" -name "*_Grn.idat" -o -name "*_Grn.idat.gz" | sort | while read -r grn; do
            red="${grn/_Grn.idat/_Red.idat}"
            if [[ -f "${red}" ]] || [[ -f "${red}.gz" ]]; then
                # Extract sample ID from barcode_position pattern
                base=$(basename "${grn}" | sed 's/_Grn\.idat\(\.gz\)\?$//')
                echo -e "${base}\t$(basename "${grn}")\t$(basename "${red}")"
            fi
        done
    } > "${SAMPLE_SHEET}"
    n_samples=$(( $(wc -l < "${SAMPLE_SHEET}") - 1 ))
    echo "Found ${n_samples} samples"
fi

# ---------------------------------------------------------------
# Step 2: Convert IDAT files to GTC files
# ---------------------------------------------------------------
echo ""
echo "--- Step 1/3: Converting IDAT files to GTC ---"
STEP_START=${SECONDS}

# Collect all green and red IDAT files
GREEN_IDATS=()
RED_IDATS=()
while IFS=$'\t' read -r sample_id grn red; do
    [[ "${sample_id}" == "sample_id" ]] && continue
    grn_path="${IDAT_DIR}/${grn}"
    red_path="${IDAT_DIR}/${red}"
    if [[ -f "${grn_path}" ]] && [[ -f "${red_path}" ]]; then
        GREEN_IDATS+=("${grn_path}")
        RED_IDATS+=("${red_path}")
    elif [[ -f "${grn_path}.gz" ]] && [[ -f "${red_path}.gz" ]]; then
        GREEN_IDATS+=("${grn_path}.gz")
        RED_IDATS+=("${red_path}.gz")
    else
        echo "Warning: Skipping ${sample_id} - IDAT files not found" >&2
    fi
done < "${SAMPLE_SHEET}"

n_total=${#GREEN_IDATS[@]}
echo "Processing ${n_total} samples in batches of ${BATCH_SIZE}..."

for (( i=0; i<n_total; i+=BATCH_SIZE )); do
    batch_end=$(( i + BATCH_SIZE ))
    if (( batch_end > n_total )); then
        batch_end=${n_total}
    fi
    batch_num=$(( i / BATCH_SIZE + 1 ))
    echo "  Batch ${batch_num}: samples $((i+1))-${batch_end}"

    # Build file arguments for this batch
    GRN_ARGS=()
    RED_ARGS=()
    for (( j=i; j<batch_end; j++ )); do
        GRN_ARGS+=("${GREEN_IDATS[j]}")
        RED_ARGS+=("${RED_IDATS[j]}")
    done

    bcftools +idat2gtc \
        --bpm "${BPM}" \
        --egt "${EGT}" \
        --output "${GTC_DIR}" \
        "${GRN_ARGS[@]}" 2>&1 | tail -1
done

echo "GTC conversion complete. Files in: ${GTC_DIR}"
STEP_ELAPSED=$(( SECONDS - STEP_START ))
echo "  Step 1/3 completed in ${STEP_ELAPSED}s"

# ---------------------------------------------------------------
# Step 3: Convert GTC files to VCF
# ---------------------------------------------------------------
echo ""
echo "--- Step 2/3: Converting GTC files to VCF ---"
STEP_START=${SECONDS}

VCF_OUTPUT="${VCF_DIR}/stage1_initial.bcf"
EXTRA_TSV="${QC_DIR}/gtc_metadata.tsv"

bcftools +gtc2vcf \
    --no-version -Ou \
    --bpm "${BPM}" \
    --csv "${CSV}" \
    --egt "${EGT}" \
    --fasta-ref "${REF_FASTA}" \
    --gtcs "${GTC_DIR}" \
    --extra "${EXTRA_TSV}" \
    --threads "${THREADS}" | \
bcftools sort -Ou -T "${OUTPUT_DIR}/bcftools." | \
bcftools norm --no-version -Ob -c x -f "${REF_FASTA}" \
    -o "${VCF_OUTPUT}" --write-index

echo "VCF created: ${VCF_OUTPUT}"
STEP_ELAPSED=$(( SECONDS - STEP_START ))
echo "  Step 2/3 completed in ${STEP_ELAPSED}s"

# ---------------------------------------------------------------
# Step 4: Compute per-sample QC metrics
# ---------------------------------------------------------------
echo ""
echo "--- Step 3/3: Computing QC metrics ---"
STEP_START=${SECONDS}

source "${SCRIPT_DIR}/collect_qc_metrics.sh"
collect_qc_metrics "${VCF_OUTPUT}" "${EXTRA_TSV}" "${QC_DIR}/stage1_sample_qc.tsv"
STEP_ELAPSED=$(( SECONDS - STEP_START ))
echo "  Step 3/3 completed in ${STEP_ELAPSED}s"

echo ""
echo "============================================"
echo "Stage 1 complete!"
echo "============================================"
echo "  VCF:           ${VCF_OUTPUT}"
echo "  QC metrics:    ${QC_DIR}/stage1_sample_qc.tsv"
echo "  GTC metadata:  ${EXTRA_TSV}"
echo ""
echo "Review QC metrics before proceeding to Stage 2."
