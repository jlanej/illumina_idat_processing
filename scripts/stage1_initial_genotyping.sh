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
  --skip-failures        Continue past corrupt/truncated IDAT files instead of
                         halting.  Failed samples are logged to
                         <output-dir>/qc/failed_idat2gtc.tsv and excluded from
                         downstream outputs.  Without this flag the pipeline
                         stops immediately on the first IDAT read error.
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
SKIP_FAILURES="false"

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
        --skip-failures) SKIP_FAILURES="true"; shift ;;
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
echo "--- Step 1/4: Converting IDAT files to GTC ---"

if step_done "stage1_idat2gtc"; then
    n_existing=$(find "${GTC_DIR}" -name '*.gtc' 2>/dev/null | wc -l)
    echo "  [checkpoint] GTC conversion already complete (${n_existing} files). Skipping."
else
    STEP_START=${SECONDS}

    # Failed-sample report — records every sample that could not be converted
    # (e.g. truncated/corrupt IDAT files).  Downstream steps use this to
    # understand why certain samples are absent from the final VCF.
    FAILED_LOG="${QC_DIR}/failed_idat2gtc.tsv"
    echo -e "sample_id\tgreen_idat\tred_idat\terror\ttimestamp" > "${FAILED_LOG}"

    # Collect all green and red IDAT files, skipping samples with existing GTC
    GREEN_IDATS=()
    RED_IDATS=()
    SAMPLE_IDS=()
    n_skipped=0
    while IFS=$'\t' read -r sample_id grn red; do
        [[ "${sample_id}" == "sample_id" ]] && continue
        grn_path="${IDAT_DIR}/${grn}"
        red_path="${IDAT_DIR}/${red}"

        # Check if GTC output already exists for this sample
        gtc_file="${GTC_DIR}/${sample_id}.gtc"
        if [[ -f "${gtc_file}" && -s "${gtc_file}" ]]; then
            (( n_skipped++ )) || true
            continue
        fi

        if [[ -f "${grn_path}" ]] && [[ -f "${red_path}" ]]; then
            GREEN_IDATS+=("${grn_path}")
            RED_IDATS+=("${red_path}")
            SAMPLE_IDS+=("${sample_id}")
        elif [[ -f "${grn_path}.gz" ]] && [[ -f "${red_path}.gz" ]]; then
            GREEN_IDATS+=("${grn_path}.gz")
            RED_IDATS+=("${red_path}.gz")
            SAMPLE_IDS+=("${sample_id}")
        else
            echo "Warning: Skipping ${sample_id} - IDAT files not found" >&2
        fi
    done < "${SAMPLE_SHEET}"

    n_total=${#GREEN_IDATS[@]}
    n_failed=0
    if [[ "${n_skipped}" -gt 0 ]]; then
        echo "  Skipping ${n_skipped} samples with existing GTC files."
    fi

    if [[ "${n_total}" -eq 0 ]]; then
        echo "  All samples already have GTC files. Nothing to convert."
    else
        echo "Processing ${n_total} samples in batches of ${BATCH_SIZE}..."

        for (( i=0; i<n_total; i+=BATCH_SIZE )); do
            batch_end=$(( i + BATCH_SIZE ))
            if (( batch_end > n_total )); then
                batch_end=${n_total}
            fi
            batch_num=$(( i / BATCH_SIZE + 1 ))
            echo "  Batch ${batch_num}: samples $((i+1))-${batch_end}"

            # Build interleaved green/red file arguments for this batch.
            # idat2gtc expects pairs: <green1.idat> <red1.idat> …
            IDAT_ARGS=()
            for (( j=i; j<batch_end; j++ )); do
                IDAT_ARGS+=("${GREEN_IDATS[j]}")
                IDAT_ARGS+=("${RED_IDATS[j]}")
            done

            BATCH_LOG=$(mktemp)
            if bcftools +idat2gtc \
                --bpm "${BPM}" \
                --egt "${EGT}" \
                --output "${GTC_DIR}" \
                "${IDAT_ARGS[@]}" > "${BATCH_LOG}" 2>&1; then
                tail -1 "${BATCH_LOG}"
            else
                BATCH_RC=$?
                # Batch failed — retry each sample individually to isolate
                # the problematic file(s).
                echo "  Batch ${batch_num} failed (exit ${BATCH_RC}). Retrying samples individually..."
                for (( j=i; j<batch_end; j++ )); do
                    SAMPLE_LOG=$(mktemp)
                    if bcftools +idat2gtc \
                        --bpm "${BPM}" \
                        --egt "${EGT}" \
                        --output "${GTC_DIR}" \
                        "${GREEN_IDATS[j]}" "${RED_IDATS[j]}" > "${SAMPLE_LOG}" 2>&1; then
                        : # success
                    else
                        SAMPLE_RC=$?
                        ERROR_MSG=$(tail -1 "${SAMPLE_LOG}" | tr '\t' ' ')
                        (( n_failed++ )) || true

                        if [[ "${SKIP_FAILURES}" == "true" ]]; then
                            echo "  WARNING: Skipping ${SAMPLE_IDS[j]} (exit ${SAMPLE_RC}): ${ERROR_MSG}" >&2
                            echo -e "${SAMPLE_IDS[j]}\t$(basename "${GREEN_IDATS[j]}")\t$(basename "${RED_IDATS[j]}")\t${ERROR_MSG}\t$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
                                >> "${FAILED_LOG}"
                            # Remove any partial/empty GTC left behind
                            rm -f "${GTC_DIR}/${SAMPLE_IDS[j]}.gtc"
                        else
                            echo "  ERROR: Failed to convert ${SAMPLE_IDS[j]} (exit ${SAMPLE_RC}): ${ERROR_MSG}" >&2
                            echo "  The IDAT file may be corrupt or truncated." >&2
                            echo "  To skip failures and continue, re-run with --skip-failures" >&2
                            rm -f "${SAMPLE_LOG}" "${BATCH_LOG}"
                            exit 1
                        fi
                    fi
                    rm -f "${SAMPLE_LOG}"
                done
            fi
            rm -f "${BATCH_LOG}"
        done
    fi

    # Report failed samples
    n_failed_logged=$(( $(wc -l < "${FAILED_LOG}") - 1 ))
    if [[ "${n_failed_logged}" -gt 0 ]]; then
        echo ""
        echo "  ============================================"
        echo "  WARNING: ${n_failed_logged} sample(s) failed IDAT-to-GTC conversion"
        echo "  ============================================"
        echo "  These samples had corrupt or truncated IDAT files and were skipped."
        echo "  They will be absent from downstream VCF and QC outputs."
        echo "  Failed sample report: ${FAILED_LOG}"
        echo ""
        echo "  Failed samples:"
        tail -n +2 "${FAILED_LOG}" | while IFS=$'\t' read -r sid grn red err ts; do
            echo "    - ${sid}: ${err}"
        done
        echo ""
    else
        # No failures — remove the empty header-only log
        rm -f "${FAILED_LOG}"
    fi

    n_gtc_created=$(find "${GTC_DIR}" -name '*.gtc' 2>/dev/null | wc -l)
    echo "GTC conversion complete. ${n_gtc_created} files in: ${GTC_DIR}"
    STEP_ELAPSED=$(( SECONDS - STEP_START ))
    echo "  Step 1/4 completed in ${STEP_ELAPSED}s"
    mark_done "stage1_idat2gtc"
fi

# ---------------------------------------------------------------
# Step 3: Convert GTC files to VCF
# ---------------------------------------------------------------
echo ""
echo "--- Step 2/4: Converting GTC files to VCF ---"

VCF_OUTPUT="${VCF_DIR}/stage1_initial.bcf"
EXTRA_TSV="${QC_DIR}/gtc_metadata.tsv"

if step_done "stage1_gtc2vcf" && [[ -f "${VCF_OUTPUT}" ]]; then
    echo "  [checkpoint] VCF conversion already complete. Skipping."
    echo "  VCF: ${VCF_OUTPUT}"
else
    STEP_START=${SECONDS}

    # Clean up any leftover temp files from a previous interrupted run
    rm -f "${VCF_OUTPUT}.tmp" "${VCF_OUTPUT}.tmp.csi"

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
    echo "  Step 2/4 completed in ${STEP_ELAPSED}s"
    mark_done "stage1_gtc2vcf"
fi

# ---------------------------------------------------------------
# Step 4: Compute per-sample QC metrics
# ---------------------------------------------------------------
echo ""
echo "--- Step 3/4: Computing QC metrics ---"

QC_OUTPUT="${QC_DIR}/stage1_sample_qc.tsv"

if step_done "stage1_qc" && [[ -f "${QC_OUTPUT}" ]]; then
    n_qc=$(( $(wc -l < "${QC_OUTPUT}") - 1 ))
    echo "  [checkpoint] QC metrics already computed (${n_qc} samples). Skipping."
else
    STEP_START=${SECONDS}

    source "${SCRIPT_DIR}/collect_qc_metrics.sh"
    collect_qc_metrics "${VCF_OUTPUT}" "${EXTRA_TSV}" "${QC_OUTPUT}" "${THREADS}"
    STEP_ELAPSED=$(( SECONDS - STEP_START ))
    echo "  Step 3/4 completed in ${STEP_ELAPSED}s"
    mark_done "stage1_qc"
fi

# ---------------------------------------------------------------
# Step 4: Compute variant-level QC metrics
# ---------------------------------------------------------------
echo ""
echo "--- Step 4/4: Computing variant QC metrics ---"

VARIANT_QC_DIR="${QC_DIR}/variant_qc"

if step_done "stage1_variant_qc" && [[ -d "${VARIANT_QC_DIR}" ]]; then
    echo "  [checkpoint] Variant QC already computed. Skipping."
else
    STEP_START=${SECONDS}

    bash "${SCRIPT_DIR}/compute_variant_qc.sh" \
        --vcf "${VCF_OUTPUT}" \
        --output-dir "${VARIANT_QC_DIR}" \
        --threads "${THREADS}"

    STEP_ELAPSED=$(( SECONDS - STEP_START ))
    echo "  Step 4/4 completed in ${STEP_ELAPSED}s"
    mark_done "stage1_variant_qc"
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
