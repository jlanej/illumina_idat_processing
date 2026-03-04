#!/usr/bin/env bash
#
# stage2_recluster.sh
#
# Stage 2: Recluster and reprocess Illumina array data.
#
# Uses QC metrics from Stage 1 to identify high-quality samples (by call
# rate and LRR standard deviation), then recomputes the EGT genotype
# cluster file from the high-quality samples' intensity data. This new
# EGT is used to re-call genotypes (via idat2gtc) for all study samples,
# producing more accurate genotype calls and intensity estimates.
#
# This approach is fundamentally different from just using --adjust-clusters
# in gtc2vcf (which only median-adjusts BAF/LRR post-hoc). By recomputing
# the actual EGT cluster definitions, we improve the underlying genotype
# calls themselves, since the GenCall algorithm uses the EGT to assign
# each sample's intensities to the correct genotype cluster (AA/AB/BB).
#
# Inspired by the MoChA/mochawdl pipeline approach of iterative QC.
#
# Idempotent: each step uses a .done sentinel file so it can be safely
# re-run after interruption without repeating completed work.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Defaults
MIN_CALL_RATE=0.97
MAX_LRR_SD=0.35
THREADS=1
MIN_CLUSTER_SAMPLES=5
BATCH_SIZE=256

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Stage 2: Recompute EGT clusters from high-quality samples and reprocess
all samples with improved genotype calling.

Required:
  --idat-dir DIR              Directory containing IDAT files
  --bpm FILE                  Illumina BPM manifest file
  --egt FILE                  Illumina EGT cluster file (original)
  --csv FILE                  Illumina CSV manifest file
  --ref-fasta FILE            Reference genome FASTA file
  --stage1-dir DIR            Stage 1 output directory (contains gtc/ and qc/)
  --output-dir DIR            Output directory for Stage 2

Options:
  --min-call-rate FLOAT       Minimum call rate for HQ samples (default: ${MIN_CALL_RATE})
  --max-lrr-sd FLOAT          Maximum LRR SD for HQ samples (default: ${MAX_LRR_SD})
  --min-cluster-samples INT   Minimum samples per genotype to recompute cluster (default: ${MIN_CLUSTER_SAMPLES})
  --sample-name-map FILE      Two-column tab-delimited file mapping IDAT root
                              names to desired sample names (renames GTC files)
  --threads INT               Number of threads (default: ${THREADS})
  --batch-size INT            Batch size for IDAT to GTC conversion (default: ${BATCH_SIZE})
  --skip-failures             Continue past corrupt/truncated IDAT files instead of
                              halting.  Failed samples are logged to
                              <output-dir>/qc/failed_idat2gtc_stage2.tsv.
  --force                     Force re-run of all steps, ignoring checkpoints
  --help                      Show this help message
EOF
    exit 0
}

IDAT_DIR=""
BPM=""
EGT=""
CSV=""
REF_FASTA=""
STAGE1_DIR=""
OUTPUT_DIR=""
SAMPLE_NAME_MAP=""
FORCE="false"
SKIP_FAILURES="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --idat-dir)       IDAT_DIR="$2"; shift 2 ;;
        --bpm)            BPM="$2"; shift 2 ;;
        --egt)            EGT="$2"; shift 2 ;;
        --csv)            CSV="$2"; shift 2 ;;
        --ref-fasta)      REF_FASTA="$2"; shift 2 ;;
        --stage1-dir)     STAGE1_DIR="$2"; shift 2 ;;
        --output-dir)     OUTPUT_DIR="$2"; shift 2 ;;
        --sample-name-map) SAMPLE_NAME_MAP="$2"; shift 2 ;;
        --min-call-rate)  MIN_CALL_RATE="$2"; shift 2 ;;
        --max-lrr-sd)     MAX_LRR_SD="$2"; shift 2 ;;
        --min-cluster-samples) MIN_CLUSTER_SAMPLES="$2"; shift 2 ;;
        --threads)        THREADS="$2"; shift 2 ;;
        --batch-size)     BATCH_SIZE="$2"; shift 2 ;;
        --skip-failures)  SKIP_FAILURES="true"; shift ;;
        --force)          FORCE="true"; shift ;;
        --help)           usage ;;
        *)                echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Validate required arguments
for var in IDAT_DIR BPM EGT CSV REF_FASTA STAGE1_DIR OUTPUT_DIR; do
    if [[ -z "${!var}" ]]; then
        echo "Error: --$(echo ${var} | tr '[:upper:]' '[:lower:]' | tr '_' '-') is required" >&2
        exit 1
    fi
done

STAGE1_QC="${STAGE1_DIR}/qc/stage1_sample_qc.tsv"
STAGE1_GTC_DIR="${STAGE1_DIR}/gtc"
if [[ ! -f "${STAGE1_QC}" ]]; then
    echo "Error: Stage 1 QC file not found: ${STAGE1_QC}" >&2
    echo "Run stage1_initial_genotyping.sh first." >&2
    exit 1
fi
if [[ ! -d "${STAGE1_GTC_DIR}" ]]; then
    echo "Error: Stage 1 GTC directory not found: ${STAGE1_GTC_DIR}" >&2
    echo "Run stage1_initial_genotyping.sh first." >&2
    exit 1
fi

GTC_DIR="${OUTPUT_DIR}/gtc"
VCF_DIR="${OUTPUT_DIR}/vcf"
QC_DIR="${OUTPUT_DIR}/qc"
CLUSTER_DIR="${OUTPUT_DIR}/clusters"
CHECKPOINT_DIR="${OUTPUT_DIR}/.checkpoints"
mkdir -p "${GTC_DIR}" "${VCF_DIR}" "${QC_DIR}" "${CLUSTER_DIR}" "${CHECKPOINT_DIR}"

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
echo "Stage 2: Recluster and Reprocess"
echo "============================================"
echo "IDAT dir:              ${IDAT_DIR}"
echo "Stage 1 dir:           ${STAGE1_DIR}"
echo "Output dir:            ${OUTPUT_DIR}"
echo "Min call rate:         ${MIN_CALL_RATE}"
echo "Max LRR SD:            ${MAX_LRR_SD}"
echo "Min cluster samples:   ${MIN_CLUSTER_SAMPLES}"
echo "Threads:               ${THREADS}"
echo "============================================"
echo ""

# ---------------------------------------------------------------
# Step 1: Identify high-quality samples from Stage 1 QC metrics
# ---------------------------------------------------------------
echo "--- Step 1/7: Identifying high-quality samples ---"

HQ_SAMPLES="${QC_DIR}/high_quality_samples.txt"
EXCLUDED_SAMPLES="${QC_DIR}/excluded_samples.txt"

if step_done "stage2_hq_select" && [[ -f "${HQ_SAMPLES}" ]]; then
    n_hq=$(wc -l < "${HQ_SAMPLES}")
    echo "  [checkpoint] HQ sample selection already complete (${n_hq} samples). Skipping."
    echo ""
else
    awk -F'\t' -v min_cr="${MIN_CALL_RATE}" -v max_sd="${MAX_LRR_SD}" '
        NR == 1 { next }
        {
            sample = $1
            call_rate = $2
            lrr_sd = $3
            if (call_rate+0 >= min_cr+0 && (lrr_sd == "NA" || lrr_sd+0 <= max_sd+0)) {
                print sample > "/dev/fd/3"
            } else {
                print sample > "/dev/fd/4"
            }
        }
    ' "${STAGE1_QC}" 3>"${HQ_SAMPLES}" 4>"${EXCLUDED_SAMPLES}"

    n_total=$(( $(wc -l < "${STAGE1_QC}") - 1 ))
    n_hq=$(wc -l < "${HQ_SAMPLES}")
    n_excluded=$(wc -l < "${EXCLUDED_SAMPLES}")

    echo "  Total samples:            ${n_total}"
    echo "  High-quality samples:     ${n_hq}"
    echo "  Excluded samples:         ${n_excluded}"

    # Diagnostic: show QC threshold evaluation details
    echo "  [diag] QC thresholds: call_rate >= ${MIN_CALL_RATE}, lrr_sd <= ${MAX_LRR_SD}"
    echo "  [diag] Stage 1 QC file first 5 lines:"
    head -5 "${STAGE1_QC}" | sed 's/^/    [diag]   /'
    if [[ "${n_excluded}" -gt 0 ]]; then
        echo "  [diag] First 5 excluded samples with values:"
        awk -F'\t' -v min_cr="${MIN_CALL_RATE}" -v max_sd="${MAX_LRR_SD}" '
            NR == 1 { next }
            {
                cr = $2; sd = $3
                fail_cr = (cr+0 < min_cr+0) ? "FAIL" : "ok"
                fail_sd = (sd != "NA" && sd+0 > max_sd+0) ? "FAIL" : "ok"
                if (fail_cr == "FAIL" || fail_sd == "FAIL") {
                    n++
                    if (n <= 5) printf "    [diag]   %s: call_rate=%s(%s) lrr_sd=%s(%s)\n", $1, cr, fail_cr, sd, fail_sd
                }
            }' "${STAGE1_QC}"
    fi
    echo ""

    if [[ "${n_hq}" -lt 10 ]]; then
        echo "Warning: Only ${n_hq} high-quality samples found. Reclustering may not be reliable." >&2
        echo "Consider relaxing --min-call-rate or --max-lrr-sd thresholds." >&2
    fi

    mark_done "stage2_hq_select"
fi

# ---------------------------------------------------------------
# Step 2: Recompute EGT cluster file from high-quality samples
# ---------------------------------------------------------------
echo "--- Step 2/7: Recomputing EGT clusters from high-quality samples ---"

RECLUSTERED_EGT="${CLUSTER_DIR}/reclustered.egt"

if step_done "stage2_recluster" && [[ -f "${RECLUSTERED_EGT}" ]]; then
    echo "  [checkpoint] EGT reclustering already complete. Skipping."
else
    STEP_START=${SECONDS}

    python3 "${SCRIPT_DIR}/recluster_egt.py" \
        --egt "${EGT}" \
        --bpm "${BPM}" \
        --gtc-dir "${STAGE1_GTC_DIR}" \
        --hq-samples "${HQ_SAMPLES}" \
        --output-egt "${RECLUSTERED_EGT}" \
        --min-cluster-samples "${MIN_CLUSTER_SAMPLES}" \
        --workers "${THREADS}"

    echo "  Reclustered EGT: ${RECLUSTERED_EGT}"
    STEP_ELAPSED=$(( SECONDS - STEP_START ))
    echo "  Step 2/7 completed in ${STEP_ELAPSED}s"
    mark_done "stage2_recluster"
fi

# ---------------------------------------------------------------
# Step 3: Re-run idat2gtc with the new EGT clusters
# ---------------------------------------------------------------
echo ""
echo "--- Step 3/7: Re-genotyping all samples with new clusters ---"

if step_done "stage2_idat2gtc"; then
    n_existing=$(find "${GTC_DIR}" -name '*.gtc' 2>/dev/null | wc -l)
    echo "  [checkpoint] Re-genotyping already complete (${n_existing} files). Skipping."
else
    STEP_START=${SECONDS}

    # Failed-sample report for stage 2
    FAILED_LOG="${QC_DIR}/failed_idat2gtc_stage2.tsv"
    echo -e "sample_id\tgreen_idat\tred_idat\terror\ttimestamp" > "${FAILED_LOG}"

    # Build interleaved green/red IDAT pairs, skipping samples with existing GTC.
    # idat2gtc expects pairs: <green1.idat> <red1.idat> …
    ALL_GREEN=()
    ALL_RED=()
    ALL_IDS=()
    n_skipped=0
    while IFS= read -r -d '' grn; do
        red="${grn/_Grn.idat/_Red.idat}"
        if [[ -f "${red}" ]]; then
            # Check if GTC output already exists for this sample
            sample_id=$(basename "${grn}" | sed 's/_Grn\.idat\(\.gz\)\?$//')
            gtc_file="${GTC_DIR}/${sample_id}.gtc"
            if [[ -f "${gtc_file}" && -s "${gtc_file}" ]]; then
                (( n_skipped++ )) || true
                continue
            fi
            ALL_GREEN+=("${grn}")
            ALL_RED+=("${red}")
            ALL_IDS+=("${sample_id}")
        else
            echo "Warning: Red IDAT not found for $(basename "${grn}"), skipping" >&2
        fi
    done < <(find "${IDAT_DIR}" \( -name "*_Grn.idat" -o -name "*_Grn.idat.gz" \) -print0 | sort -z)

    n_pairs=${#ALL_GREEN[@]}
    n_failed=0
    if [[ "${n_skipped}" -gt 0 ]]; then
        echo "  Skipping ${n_skipped} samples with existing GTC files."
    fi

    if [[ "${n_pairs}" -eq 0 ]]; then
        echo "  All samples already have GTC files. Nothing to convert."
    else
        echo "  Re-calling genotypes for ${n_pairs} samples with reclustered EGT..."
        echo "  Processing in batches of ${BATCH_SIZE}..."

        for (( i=0; i<n_pairs; i+=BATCH_SIZE )); do
            batch_end=$(( i + BATCH_SIZE ))
            if (( batch_end > n_pairs )); then
                batch_end=${n_pairs}
            fi
            batch_num=$(( i / BATCH_SIZE + 1 ))
            echo "  Batch ${batch_num}: samples $((i+1))-${batch_end}"

            IDAT_ARGS=()
            for (( j=i; j<batch_end; j++ )); do
                IDAT_ARGS+=("${ALL_GREEN[j]}")
                IDAT_ARGS+=("${ALL_RED[j]}")
            done

            BATCH_LOG=$(mktemp)
            if bcftools +idat2gtc \
                --bpm "${BPM}" \
                --egt "${RECLUSTERED_EGT}" \
                --output "${GTC_DIR}" \
                "${IDAT_ARGS[@]}" > "${BATCH_LOG}" 2>&1; then
                tail -1 "${BATCH_LOG}"
            else
                BATCH_RC=$?
                echo "  Batch ${batch_num} failed (exit ${BATCH_RC}). Retrying samples individually..."
                for (( j=i; j<batch_end; j++ )); do
                    SAMPLE_LOG=$(mktemp)
                    if bcftools +idat2gtc \
                        --bpm "${BPM}" \
                        --egt "${RECLUSTERED_EGT}" \
                        --output "${GTC_DIR}" \
                        "${ALL_GREEN[j]}" "${ALL_RED[j]}" > "${SAMPLE_LOG}" 2>&1; then
                        : # success
                    else
                        SAMPLE_RC=$?
                        ERROR_MSG=$(tail -1 "${SAMPLE_LOG}" | tr '\t' ' ')
                        (( n_failed++ )) || true

                        if [[ "${SKIP_FAILURES}" == "true" ]]; then
                            echo "  WARNING: Skipping ${ALL_IDS[j]} (exit ${SAMPLE_RC}): ${ERROR_MSG}" >&2
                            echo -e "${ALL_IDS[j]}\t$(basename "${ALL_GREEN[j]}")\t$(basename "${ALL_RED[j]}")\t${ERROR_MSG}\t$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
                                >> "${FAILED_LOG}"
                            rm -f "${GTC_DIR}/${ALL_IDS[j]}.gtc"
                        else
                            echo "  ERROR: Failed to convert ${ALL_IDS[j]} (exit ${SAMPLE_RC}): ${ERROR_MSG}" >&2
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
        echo "  WARNING: ${n_failed_logged} sample(s) failed IDAT-to-GTC conversion (stage 2)"
        echo "  ============================================"
        echo "  These samples had corrupt or truncated IDAT files and were skipped."
        echo "  They will be absent from downstream VCF and QC outputs."
        echo "  Failed sample report: ${FAILED_LOG}"
        echo ""
        echo "  Failed samples:"
        tail -n +2 "${FAILED_LOG}" | while IFS=$'\t' read -r sid grn red err _; do
            echo "    - ${sid}: ${err}"
        done
        echo ""
    else
        rm -f "${FAILED_LOG}"
    fi

    n_gtc_created=$(find "${GTC_DIR}" -name '*.gtc' 2>/dev/null | wc -l)
    echo "  GTC re-calling complete. ${n_gtc_created} files in: ${GTC_DIR}"
    STEP_ELAPSED=$(( SECONDS - STEP_START ))
    echo "  Step 3/7 completed in ${STEP_ELAPSED}s"
    mark_done "stage2_idat2gtc"
fi

# ---------------------------------------------------------------
# Rename GTC files using sample name map (if provided)
# ---------------------------------------------------------------
STAGE2_NAME_MAP_FILE="${OUTPUT_DIR}/sample_name_mapping.tsv"

if [[ -n "${SAMPLE_NAME_MAP}" && -f "${SAMPLE_NAME_MAP}" ]]; then
    if step_done "stage2_rename_gtc" && [[ -f "${STAGE2_NAME_MAP_FILE}" ]]; then
        n_renamed=$(awk -F'\t' '$1 != $2' "${STAGE2_NAME_MAP_FILE}" | wc -l)
        echo "  [checkpoint] GTC rename already complete (${n_renamed} renamed). Skipping."
    else
        echo ""
        echo "--- Applying sample name map (stage 2) ---"

        declare -A S2_MAP_OLD_TO_NEW
        while IFS=$'\t' read -r old_name new_name; do
            [[ -z "${old_name}" || "${old_name}" == "#"* ]] && continue
            S2_MAP_OLD_TO_NEW["${old_name}"]="${new_name}"
        done < "${SAMPLE_NAME_MAP}"

        # Build list of GTC files to process
        S2_GTC_SAMPLES=()
        while IFS= read -r gtc_path; do
            S2_GTC_SAMPLES+=("$(basename "${gtc_path}" .gtc)")
        done < <(find "${GTC_DIR}" -name '*.gtc' | sort)

        {
            echo -e "sample_id\toriginal_name"
            n_renamed=0
            for gtc_id in "${S2_GTC_SAMPLES[@]}"; do
                if [[ -n "${S2_MAP_OLD_TO_NEW[${gtc_id}]+x}" ]]; then
                    new_name="${S2_MAP_OLD_TO_NEW[${gtc_id}]}"
                    if [[ "${new_name}" != "${gtc_id}" ]]; then
                        mv "${GTC_DIR}/${gtc_id}.gtc" "${GTC_DIR}/${new_name}.gtc"
                        (( n_renamed++ )) || true
                    fi
                    echo -e "${new_name}\t${gtc_id}"
                else
                    echo -e "${gtc_id}\t${gtc_id}"
                fi
            done
        } > "${STAGE2_NAME_MAP_FILE}"

        echo "  Renamed ${n_renamed} GTC files"
        echo "  Name mapping saved: ${STAGE2_NAME_MAP_FILE}"
        mark_done "stage2_rename_gtc"
    fi
else
    # No name map: write identity mapping for downstream consistency
    if [[ ! -f "${STAGE2_NAME_MAP_FILE}" ]]; then
        {
            echo -e "sample_id\toriginal_name"
            find "${GTC_DIR}" -name '*.gtc' | sort | while read -r gtc_path; do
                sid=$(basename "${gtc_path}" .gtc)
                echo -e "${sid}\t${sid}"
            done
        } > "${STAGE2_NAME_MAP_FILE}"
    fi
fi

# ---------------------------------------------------------------
# Step 4: Convert reclustered GTC files to VCF
# ---------------------------------------------------------------
echo ""
echo "--- Step 4/7: Converting reclustered GTC files to VCF ---"

VCF_OUTPUT="${VCF_DIR}/stage2_reclustered.bcf"
EXTRA_TSV="${QC_DIR}/gtc_metadata_stage2.tsv"

if step_done "stage2_gtc2vcf" && [[ -f "${VCF_OUTPUT}" ]]; then
    echo "  [checkpoint] VCF conversion already complete. Skipping."
else
    STEP_START=${SECONDS}

    # Clean up any leftover temp files from a previous interrupted run
    rm -f "${VCF_OUTPUT}.tmp" "${VCF_OUTPUT}.tmp.csi"

    # Convert reclustered GTC to VCF (pre-normalization)
    PRE_NORM_BCF="${VCF_DIR}/stage2_pre_norm.bcf"

    # Dynamically allocate sort memory: use ~50% of available RAM (capped at 64G)
    SORT_MEM="4G"
    if [[ -f /proc/meminfo ]]; then
        SORT_MEM=$(awk '/MemAvailable/{
            mem_gb = int($2 / 1024 / 1024 * 0.5)
            if (mem_gb < 4) mem_gb = 4
            if (mem_gb > 64) mem_gb = 64
            printf "%dG", mem_gb
        }' /proc/meminfo)
    fi
    echo "  Sort memory: ${SORT_MEM}"

    bcftools +gtc2vcf \
        --no-version -Ou \
        --do-not-check-bpm \
        --bpm "${BPM}" \
        --csv "${CSV}" \
        --egt "${RECLUSTERED_EGT}" \
        --fasta-ref "${REF_FASTA}" \
        --gtcs "${GTC_DIR}" \
        --extra "${EXTRA_TSV}" \
        --adjust-clusters \
        --threads "${THREADS}" | \
    bcftools sort -Ob -m "${SORT_MEM}" -T "${OUTPUT_DIR}/bcftools." \
        -o "${PRE_NORM_BCF}" --write-index

    # Diagnostic: count variants before normalization
    N_PRE_NORM=$(bcftools index -n "${PRE_NORM_BCF}" 2>/dev/null || \
        bcftools view -H "${PRE_NORM_BCF}" | wc -l)
    echo "  Variants before normalization: ${N_PRE_NORM}"

    # Normalize with -c ws (warn and swap REF/ALT to match reference).
    # CRITICAL: Do NOT use -c x here. See stage1_initial_genotyping.sh for details.
    NORM_LOG="${QC_DIR}/norm_warnings_stage2.log"

    bcftools norm --no-version -Ob -c ws -f "${REF_FASTA}" \
        --threads "${THREADS}" \
        "${PRE_NORM_BCF}" \
        -o "${VCF_OUTPUT}.tmp" --write-index 2>"${NORM_LOG}"

    # Atomic move
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
    fi
    echo "  Norm warnings log: ${NORM_LOG}"

    # Clean up pre-norm BCF
    rm -f "${PRE_NORM_BCF}" "${PRE_NORM_BCF}.csi"

    echo "  Reclustered VCF created: ${VCF_OUTPUT}"
    STEP_ELAPSED=$(( SECONDS - STEP_START ))
    echo "  Step 4/7 completed in ${STEP_ELAPSED}s"
    mark_done "stage2_gtc2vcf"
fi

# ---------------------------------------------------------------
# Step 5: Recompute QC metrics after reclustering
# ---------------------------------------------------------------
echo ""
echo "--- Step 5/7: Recomputing QC metrics ---"

STAGE2_QC="${QC_DIR}/stage2_sample_qc.tsv"

if step_done "stage2_qc" && [[ -f "${STAGE2_QC}" ]]; then
    n_qc=$(( $(wc -l < "${STAGE2_QC}") - 1 ))
    echo "  [checkpoint] QC metrics already computed (${n_qc} samples). Skipping."
else
    STEP_START=${SECONDS}

    source "${SCRIPT_DIR}/collect_qc_metrics.sh"
    collect_qc_metrics "${VCF_OUTPUT}" "${EXTRA_TSV}" "${STAGE2_QC}" "${THREADS}"

    # Add original_name column from name mapping
    if [[ -f "${STAGE2_NAME_MAP_FILE}" ]]; then
        awk -F'\t' '
            NR==FNR { if (FNR > 1) map[$1] = $2; next }
            FNR==1 { print "sample_id\toriginal_name\t" substr($0, index($0,$2)); next }
            { printf "%s\t%s", $1, ($1 in map ? map[$1] : $1)
              for (i=2; i<=NF; i++) printf "\t%s", $i
              printf "\n" }
        ' "${STAGE2_NAME_MAP_FILE}" "${STAGE2_QC}" > "${STAGE2_QC}.tmp"
        mv "${STAGE2_QC}.tmp" "${STAGE2_QC}"
    fi

    STEP_ELAPSED=$(( SECONDS - STEP_START ))
    echo "  Step 5/7 completed in ${STEP_ELAPSED}s"
    mark_done "stage2_qc"
fi

# ---------------------------------------------------------------
# Step 6: Compute variant-level QC metrics after reclustering
# ---------------------------------------------------------------
echo ""
echo "--- Step 6/7: Computing variant QC metrics ---"

VARIANT_QC_DIR="${QC_DIR}/variant_qc"

if step_done "stage2_variant_qc" && [[ -d "${VARIANT_QC_DIR}" ]]; then
    echo "  [checkpoint] Variant QC already computed. Skipping."
else
    STEP_START=${SECONDS}

    bash "${SCRIPT_DIR}/compute_variant_qc.sh" \
        --vcf "${VCF_OUTPUT}" \
        --output-dir "${VARIANT_QC_DIR}" \
        --threads "${THREADS}"

    STEP_ELAPSED=$(( SECONDS - STEP_START ))
    echo "  Step 6/7 completed in ${STEP_ELAPSED}s"
    mark_done "stage2_variant_qc"
fi

# ---------------------------------------------------------------
# Step 7: Generate QC comparison (Stage 1 vs Stage 2)
# ---------------------------------------------------------------
echo ""
echo "--- Step 7/7: Generating QC comparison ---"

COMPARISON_DIR="${QC_DIR}/comparison"
STAGE1_VARIANT_QC="${STAGE1_DIR}/qc/variant_qc"

if step_done "stage2_comparison" && [[ -d "${COMPARISON_DIR}" ]]; then
    echo "  [checkpoint] QC comparison already generated. Skipping."
else
    STEP_START=${SECONDS}

    COMPARISON_ARGS=(
        --stage1-sample-qc "${STAGE1_QC}"
        --stage2-sample-qc "${STAGE2_QC}"
        --output-dir "${COMPARISON_DIR}"
    )
    if [[ -d "${STAGE1_VARIANT_QC}" ]]; then
        COMPARISON_ARGS+=(--stage1-variant-qc-dir "${STAGE1_VARIANT_QC}")
    fi
    if [[ -d "${VARIANT_QC_DIR}" ]]; then
        COMPARISON_ARGS+=(--stage2-variant-qc-dir "${VARIANT_QC_DIR}")
    fi

    if command -v python3 &>/dev/null; then
        if ! python3 "${SCRIPT_DIR}/plot_qc_comparison.py" "${COMPARISON_ARGS[@]}" 2>&1; then
            echo "  Warning: QC comparison script failed. Check logs for details." >&2
        fi
    else
        echo "  Warning: python3 not found; skipping QC comparison plots."
    fi

    STEP_ELAPSED=$(( SECONDS - STEP_START ))
    echo "  Step 7/7 completed in ${STEP_ELAPSED}s"
    mark_done "stage2_comparison"
fi

# Generate comparison report
echo ""
echo "--- QC Comparison (Stage 1 vs Stage 2) ---"

# Use header-based column lookup (robust to column order changes)
python3 -c "
import sys
def read_qc(path):
    d = {}
    with open(path) as f:
        hdr = f.readline().strip().split('\t')
        cr_i = hdr.index('call_rate') if 'call_rate' in hdr else 1
        sd_i = hdr.index('lrr_sd') if 'lrr_sd' in hdr else 2
        for line in f:
            fs = line.strip().split('\t')
            try:
                cr = float(fs[cr_i]) if fs[cr_i] != 'NA' else None
                sd = float(fs[sd_i]) if fs[sd_i] != 'NA' else None
                d[fs[0]] = (cr, sd)
            except (ValueError, IndexError):
                pass
    return d

s1 = read_qc('${STAGE1_QC}')
s2 = read_qc('${STAGE2_QC}')
common = set(s1) & set(s2)
cr_diffs = [s2[s][0]-s1[s][0] for s in common if s1[s][0] is not None and s2[s][0] is not None]
sd_diffs = [s2[s][1]-s1[s][1] for s in common if s1[s][1] is not None and s2[s][1] is not None]
if cr_diffs:
    print(f'  Avg call rate change:  {sum(cr_diffs)/len(cr_diffs):+.4f}')
if sd_diffs:
    print(f'  Avg LRR SD change:     {sum(sd_diffs)/len(sd_diffs):+.4f}')
" 2>/dev/null || echo "  (Could not compute comparison - check output files)"

echo ""
echo "============================================"
echo "Stage 2 complete!"
echo "============================================"
echo "  Reclustered EGT:   ${RECLUSTERED_EGT}"
echo "  VCF:               ${VCF_OUTPUT}"
echo "  QC metrics:        ${QC_DIR}/stage2_sample_qc.tsv"
echo "  Variant QC:        ${VARIANT_QC_DIR}/"
echo "  QC comparison:     ${COMPARISON_DIR}/"
echo "  HQ sample list:    ${HQ_SAMPLES}"
echo "  Excluded samples:  ${EXCLUDED_SAMPLES}"
echo ""
echo "The reclustered VCF uses study-specific genotype clusters"
echo "and is ready for downstream analysis (phasing, MoChA, imputation)."
