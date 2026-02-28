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
# Idempotent: each step uses a .done sentinel file so it can be safely
# re-run after interruption without repeating completed work.
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
  --force                Force re-run of all steps, ignoring checkpoints
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
FORCE="false"

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
        --force)         FORCE="true"; shift ;;
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
CHECKPOINT_DIR="${OUTPUT_DIR}/.checkpoints"
mkdir -p "${GTC_DIR}" "${VCF_DIR}" "${QC_DIR}" "${CHECKPOINT_DIR}"

# Helper: check if a step is already complete
step_done() {
    local step_name="$1"
    [[ "${FORCE}" == "true" ]] && return 1
    [[ -f "${CHECKPOINT_DIR}/${step_name}.done" ]]
}

# Helper: mark a step as complete
mark_done() {
    local step_name="$1"
    date -u "+%Y-%m-%dT%H:%M:%SZ" > "${CHECKPOINT_DIR}/${step_name}.done"
}

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
    if [[ -f "${SAMPLE_SHEET}" && "${FORCE}" != "true" ]]; then
        n_samples=$(( $(wc -l < "${SAMPLE_SHEET}") - 1 ))
        echo "Sample sheet already exists (${n_samples} samples): ${SAMPLE_SHEET}"
    else
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
        } > "${SAMPLE_SHEET}.tmp"
        mv "${SAMPLE_SHEET}.tmp" "${SAMPLE_SHEET}"
        n_samples=$(( $(wc -l < "${SAMPLE_SHEET}") - 1 ))
        echo "Found ${n_samples} samples"
    fi
fi

# ---------------------------------------------------------------
# Step 2: Convert IDAT files to GTC files
# ---------------------------------------------------------------
echo ""
echo "--- Step 1/3: Converting IDAT files to GTC ---"

if step_done "stage1_idat2gtc"; then
    n_existing=$(find "${GTC_DIR}" -name '*.gtc' 2>/dev/null | wc -l)
    echo "  [checkpoint] GTC conversion already complete (${n_existing} files). Skipping."
else
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

        # Build interleaved green/red file arguments for this batch.
        # idat2gtc expects pairs: <green1.idat> <red1.idat> [<green2.idat> <red2.idat> ...]
        IDAT_ARGS=()
        for (( j=i; j<batch_end; j++ )); do
            IDAT_ARGS+=("${GREEN_IDATS[j]}")
            IDAT_ARGS+=("${RED_IDATS[j]}")
        done

        bcftools +idat2gtc \
            --bpm "${BPM}" \
            --egt "${EGT}" \
            --output "${GTC_DIR}" \
            "${IDAT_ARGS[@]}" 2>&1 | tail -1
    done

    echo "GTC conversion complete. Files in: ${GTC_DIR}"
    STEP_ELAPSED=$(( SECONDS - STEP_START ))
    echo "  Step 1/3 completed in ${STEP_ELAPSED}s"
    mark_done "stage1_idat2gtc"
fi

# ---------------------------------------------------------------
# Step 3: Convert GTC files to VCF
# ---------------------------------------------------------------
echo ""
echo "--- Step 2/3: Converting GTC files to VCF ---"

VCF_OUTPUT="${VCF_DIR}/stage1_initial.bcf"
EXTRA_TSV="${QC_DIR}/gtc_metadata.tsv"

if step_done "stage1_gtc2vcf" && [[ -f "${VCF_OUTPUT}" ]]; then
    echo "  [checkpoint] VCF conversion already complete. Skipping."
    echo "  VCF: ${VCF_OUTPUT}"
else
    STEP_START=${SECONDS}

    # Debug: show manifest and GTC file details for troubleshooting
    echo "  BPM file: ${BPM} ($(wc -c < "${BPM}") bytes)"
    echo "  CSV file: ${CSV} ($(wc -c < "${CSV}") bytes)"
    echo "  EGT file: ${EGT} ($(wc -c < "${EGT}") bytes)"
    echo "  GTC dir:  ${GTC_DIR} ($(find "${GTC_DIR}" -name '*.gtc' | wc -l) files)"
    echo "  Reference: ${REF_FASTA}"

    # Convert GTC to VCF (pre-normalization)
    # Note: --do-not-check-bpm is required because the realigned CSV may have
    # coordinates on a different build than the BPM (e.g., GRCh37 BPM with
    # GRCh38-realigned CSV), so the coordinate check would falsely fail.
    PRE_NORM_BCF="${VCF_DIR}/stage1_pre_norm.bcf"

    bcftools +gtc2vcf \
        --no-version -Ou \
        --do-not-check-bpm \
        --bpm "${BPM}" \
        --csv "${CSV}" \
        --egt "${EGT}" \
        --fasta-ref "${REF_FASTA}" \
        --gtcs "${GTC_DIR}" \
        --extra "${EXTRA_TSV}" \
        --threads "${THREADS}" | \
    bcftools sort -Ob -T "${OUTPUT_DIR}/bcftools." \
        -o "${PRE_NORM_BCF}" --write-index

    # Diagnostic: count variants before normalization
    N_PRE_NORM=$(bcftools index -n "${PRE_NORM_BCF}" 2>/dev/null || \
        bcftools view -H "${PRE_NORM_BCF}" | wc -l)
    echo "  Variants before normalization: ${N_PRE_NORM}"

    # Normalize with -c ws (warn and swap REF/ALT to match reference).
    # CRITICAL: Do NOT use -c x here. When processing cross-build data
    # (e.g., GRCh37 manifests on GRCh38 reference), -c x silently sets
    # genotypes to missing for all REF-mismatched probes, which can
    # deflate call rates from >97% to ~65-75% and inflate LRR SD.
    # The -c ws mode swaps REF/ALT alleles to match the reference,
    # preserving genotype calls while fixing allele orientation.
    NORM_LOG="${QC_DIR}/norm_warnings.log"

    bcftools norm --no-version -Ob -c ws -f "${REF_FASTA}" \
        "${PRE_NORM_BCF}" \
        -o "${VCF_OUTPUT}.tmp" --write-index 2>"${NORM_LOG}"

    # Atomic move: only replace final output after successful completion
    mv "${VCF_OUTPUT}.tmp" "${VCF_OUTPUT}"
    mv "${VCF_OUTPUT}.tmp.csi" "${VCF_OUTPUT}.csi" 2>/dev/null || true

    # Diagnostic: count variants after normalization and REF swaps
    N_POST_NORM=$(bcftools index -n "${VCF_OUTPUT}" 2>/dev/null || \
        bcftools view -H "${VCF_OUTPUT}" | wc -l)
    N_REF_SWAPS=$(grep -c "REF_MISMATCH" "${NORM_LOG}" 2>/dev/null || true)
    echo "  Variants after normalization:  ${N_POST_NORM}"
    echo "  REF/ALT swaps (REF mismatch):  ${N_REF_SWAPS}"
    if [[ "${N_REF_SWAPS}" -gt 0 ]]; then
        PCT_SWAPS=$(awk -v s="${N_REF_SWAPS}" -v t="${N_PRE_NORM}" 'BEGIN { printf "%.1f", (t>0) ? s*100.0/t : 0 }')
        echo "  REF swap percentage:           ${PCT_SWAPS}%"
        echo "  (REF/ALT swapped to match reference — genotype calls preserved)"
    fi
    echo "  Norm warnings log: ${NORM_LOG}"

    # Clean up pre-norm BCF
    rm -f "${PRE_NORM_BCF}" "${PRE_NORM_BCF}.csi"

    echo "VCF created: ${VCF_OUTPUT}"
    STEP_ELAPSED=$(( SECONDS - STEP_START ))
    echo "  Step 2/3 completed in ${STEP_ELAPSED}s"
    mark_done "stage1_gtc2vcf"
fi

# ---------------------------------------------------------------
# Step 4: Compute per-sample QC metrics
# ---------------------------------------------------------------
echo ""
echo "--- Step 3/3: Computing QC metrics ---"

QC_OUTPUT="${QC_DIR}/stage1_sample_qc.tsv"

if step_done "stage1_qc" && [[ -f "${QC_OUTPUT}" ]]; then
    n_qc=$(( $(wc -l < "${QC_OUTPUT}") - 1 ))
    echo "  [checkpoint] QC metrics already computed (${n_qc} samples). Skipping."
else
    STEP_START=${SECONDS}

    source "${SCRIPT_DIR}/collect_qc_metrics.sh"
    collect_qc_metrics "${VCF_OUTPUT}" "${EXTRA_TSV}" "${QC_OUTPUT}"
    STEP_ELAPSED=$(( SECONDS - STEP_START ))
    echo "  Step 3/3 completed in ${STEP_ELAPSED}s"
    mark_done "stage1_qc"
fi

echo ""
echo "============================================"
echo "Stage 1 complete!"
echo "============================================"
echo "  VCF:           ${VCF_OUTPUT}"
echo "  QC metrics:    ${QC_DIR}/stage1_sample_qc.tsv"
echo "  GTC metadata:  ${EXTRA_TSV}"
echo ""
echo "Review QC metrics before proceeding to Stage 2."
