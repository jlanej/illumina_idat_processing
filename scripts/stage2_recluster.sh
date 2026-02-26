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
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Defaults
MIN_CALL_RATE=0.97
MAX_LRR_SD=0.35
THREADS=1
MIN_CLUSTER_SAMPLES=5

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
  --threads INT               Number of threads (default: ${THREADS})
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

while [[ $# -gt 0 ]]; do
    case "$1" in
        --idat-dir)       IDAT_DIR="$2"; shift 2 ;;
        --bpm)            BPM="$2"; shift 2 ;;
        --egt)            EGT="$2"; shift 2 ;;
        --csv)            CSV="$2"; shift 2 ;;
        --ref-fasta)      REF_FASTA="$2"; shift 2 ;;
        --stage1-dir)     STAGE1_DIR="$2"; shift 2 ;;
        --output-dir)     OUTPUT_DIR="$2"; shift 2 ;;
        --min-call-rate)  MIN_CALL_RATE="$2"; shift 2 ;;
        --max-lrr-sd)     MAX_LRR_SD="$2"; shift 2 ;;
        --min-cluster-samples) MIN_CLUSTER_SAMPLES="$2"; shift 2 ;;
        --threads)        THREADS="$2"; shift 2 ;;
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
mkdir -p "${GTC_DIR}" "${VCF_DIR}" "${QC_DIR}" "${CLUSTER_DIR}"

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
echo "--- Step 1/5: Identifying high-quality samples ---"

HQ_SAMPLES="${QC_DIR}/high_quality_samples.txt"
EXCLUDED_SAMPLES="${QC_DIR}/excluded_samples.txt"

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

# ---------------------------------------------------------------
# Step 2: Recompute EGT cluster file from high-quality samples
# ---------------------------------------------------------------
echo "--- Step 2/5: Recomputing EGT clusters from high-quality samples ---"
STEP_START=${SECONDS}

RECLUSTERED_EGT="${CLUSTER_DIR}/reclustered.egt"

python3 "${SCRIPT_DIR}/recluster_egt.py" \
    --egt "${EGT}" \
    --bpm "${BPM}" \
    --gtc-dir "${STAGE1_GTC_DIR}" \
    --hq-samples "${HQ_SAMPLES}" \
    --output-egt "${RECLUSTERED_EGT}" \
    --min-cluster-samples "${MIN_CLUSTER_SAMPLES}"

echo "  Reclustered EGT: ${RECLUSTERED_EGT}"
STEP_ELAPSED=$(( SECONDS - STEP_START ))
echo "  Step 2/5 completed in ${STEP_ELAPSED}s"

# ---------------------------------------------------------------
# Step 3: Re-run idat2gtc with the new EGT clusters
# ---------------------------------------------------------------
echo ""
echo "--- Step 3/5: Re-genotyping all samples with new clusters ---"
STEP_START=${SECONDS}

ALL_GRN_IDATS=()
while IFS= read -r -d '' grn; do
    ALL_GRN_IDATS+=("${grn}")
done < <(find "${IDAT_DIR}" -name "*_Grn.idat" -print0 | sort -z)

echo "  Re-calling genotypes for ${#ALL_GRN_IDATS[@]} samples with reclustered EGT..."

bcftools +idat2gtc \
    --bpm "${BPM}" \
    --egt "${RECLUSTERED_EGT}" \
    --output "${GTC_DIR}" \
    "${ALL_GRN_IDATS[@]}" 2>&1 | tail -3

echo "  GTC re-calling complete."
STEP_ELAPSED=$(( SECONDS - STEP_START ))
echo "  Step 3/5 completed in ${STEP_ELAPSED}s"

# ---------------------------------------------------------------
# Step 4: Convert reclustered GTC files to VCF
# ---------------------------------------------------------------
echo ""
echo "--- Step 4/5: Converting reclustered GTC files to VCF ---"
STEP_START=${SECONDS}

VCF_OUTPUT="${VCF_DIR}/stage2_reclustered.bcf"
EXTRA_TSV="${QC_DIR}/gtc_metadata_stage2.tsv"

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
bcftools sort -Ou -T "${OUTPUT_DIR}/bcftools." | \
bcftools norm --no-version -Ob -c x -f "${REF_FASTA}" \
    -o "${VCF_OUTPUT}" --write-index

echo "  Reclustered VCF created: ${VCF_OUTPUT}"
STEP_ELAPSED=$(( SECONDS - STEP_START ))
echo "  Step 4/5 completed in ${STEP_ELAPSED}s"

# ---------------------------------------------------------------
# Step 5: Recompute QC metrics after reclustering
# ---------------------------------------------------------------
echo ""
echo "--- Step 5/5: Recomputing QC metrics ---"
STEP_START=${SECONDS}

source "${SCRIPT_DIR}/collect_qc_metrics.sh"
collect_qc_metrics "${VCF_OUTPUT}" "${EXTRA_TSV}" "${QC_DIR}/stage2_sample_qc.tsv"
STEP_ELAPSED=$(( SECONDS - STEP_START ))
echo "  Step 5/5 completed in ${STEP_ELAPSED}s"

# Generate comparison report
echo ""
echo "--- QC Comparison (Stage 1 vs Stage 2) ---"

paste "${STAGE1_QC}" "${QC_DIR}/stage2_sample_qc.tsv" | \
    awk -F'\t' 'NR==1 {next}
    {
        s1_cr = $2; s1_sd = $3
        s2_cr = $6; s2_sd = $7
        if (s1_cr != "NA" && s2_cr != "NA") {
            cr_diff = s2_cr - s1_cr
            cr_sum += cr_diff; cr_n++
        }
        if (s1_sd != "NA" && s2_sd != "NA") {
            sd_diff = s2_sd - s1_sd
            sd_sum += sd_diff; sd_n++
        }
    }
    END {
        if (cr_n > 0) printf "  Avg call rate change:  %+.4f\n", cr_sum/cr_n
        if (sd_n > 0) printf "  Avg LRR SD change:     %+.4f\n", sd_sum/sd_n
    }' 2>/dev/null || echo "  (Could not compute comparison - check output files)"

echo ""
echo "============================================"
echo "Stage 2 complete!"
echo "============================================"
echo "  Reclustered EGT:   ${RECLUSTERED_EGT}"
echo "  VCF:               ${VCF_OUTPUT}"
echo "  QC metrics:        ${QC_DIR}/stage2_sample_qc.tsv"
echo "  HQ sample list:    ${HQ_SAMPLES}"
echo "  Excluded samples:  ${EXCLUDED_SAMPLES}"
echo ""
echo "The reclustered VCF uses study-specific genotype clusters"
echo "and is ready for downstream analysis (phasing, MoChA, imputation)."
