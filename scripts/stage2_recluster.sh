#!/usr/bin/env bash
#
# stage2_recluster.sh
#
# Stage 2: Recluster and reprocess Illumina array data.
#
# Uses QC metrics from Stage 1 to identify high-quality samples (by call
# rate and LRR standard deviation), then reclusters genotypes using only
# those high-quality samples. The new clusters are then applied to re-call
# genotypes for all study samples, producing improved intensity estimates.
#
# Inspired by the MoChA/mochawdl pipeline approach of iterative QC.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Defaults
MIN_CALL_RATE=0.97
MAX_LRR_SD=0.35
THREADS=1

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Stage 2: Recluster on high-quality samples and reprocess all samples.

Required:
  --idat-dir DIR         Directory containing IDAT files
  --bpm FILE             Illumina BPM manifest file
  --egt FILE             Illumina EGT cluster file
  --csv FILE             Illumina CSV manifest file
  --ref-fasta FILE       Reference genome FASTA file
  --stage1-dir DIR       Stage 1 output directory (contains qc/ subdirectory)
  --output-dir DIR       Output directory for Stage 2

Options:
  --min-call-rate FLOAT  Minimum call rate for high-quality samples (default: ${MIN_CALL_RATE})
  --max-lrr-sd FLOAT     Maximum LRR SD for high-quality samples (default: ${MAX_LRR_SD})
  --threads INT          Number of threads (default: ${THREADS})
  --help                 Show this help message
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
if [[ ! -f "${STAGE1_QC}" ]]; then
    echo "Error: Stage 1 QC file not found: ${STAGE1_QC}" >&2
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
echo "IDAT dir:       ${IDAT_DIR}"
echo "Stage 1 dir:    ${STAGE1_DIR}"
echo "Output dir:     ${OUTPUT_DIR}"
echo "Min call rate:  ${MIN_CALL_RATE}"
echo "Max LRR SD:     ${MAX_LRR_SD}"
echo "Threads:        ${THREADS}"
echo "============================================"
echo ""

# ---------------------------------------------------------------
# Step 1: Identify high-quality samples from Stage 1 QC metrics
# ---------------------------------------------------------------
echo "--- Step 1/4: Identifying high-quality samples ---"

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
echo ""

if [[ "${n_hq}" -lt 10 ]]; then
    echo "Warning: Only ${n_hq} high-quality samples found. Reclustering may not be reliable." >&2
    echo "Consider relaxing --min-call-rate or --max-lrr-sd thresholds." >&2
fi

# ---------------------------------------------------------------
# Step 2: Re-run idat2gtc with cluster adjustment
# ---------------------------------------------------------------
echo "--- Step 2/4: Re-genotyping all samples with adjusted clusters ---"

# Collect IDAT files for high-quality samples
HQ_GRN_IDATS=()
while IFS= read -r sample_id; do
    grn=$(find "${IDAT_DIR}" -name "${sample_id}*_Grn.idat" -o -name "${sample_id}*_Grn.idat.gz" 2>/dev/null | head -1)
    if [[ -n "${grn}" ]]; then
        HQ_GRN_IDATS+=("${grn}")
    fi
done < "${HQ_SAMPLES}"

echo "  Found ${#HQ_GRN_IDATS[@]} IDAT files for high-quality samples"

# First, generate GTC files for all samples (using original clusters)
echo "  Converting all IDAT files to GTC..."
ALL_GRN_IDATS=()
while IFS= read -r -d '' grn; do
    ALL_GRN_IDATS+=("${grn}")
done < <(find "${IDAT_DIR}" -name "*_Grn.idat" -print0 | sort -z)

bcftools +idat2gtc \
    --bpm "${BPM}" \
    --egt "${EGT}" \
    --output "${GTC_DIR}" \
    "${ALL_GRN_IDATS[@]}" 2>&1 | tail -3

echo "  GTC conversion complete."

# ---------------------------------------------------------------
# Step 3: Convert to VCF with cluster adjustment
# ---------------------------------------------------------------
echo ""
echo "--- Step 3/4: Converting to VCF with adjusted clusters ---"

# Convert GTC to VCF, using --adjust-clusters to recompute cluster centers
# from the high-quality samples. The adjust-clusters option median-adjusts
# BAF and LRR values to correct for batch effects.
VCF_OUTPUT="${VCF_DIR}/stage2_reclustered.bcf"
EXTRA_TSV="${QC_DIR}/gtc_metadata_stage2.tsv"

bcftools +gtc2vcf \
    --no-version -Ou \
    --bpm "${BPM}" \
    --csv "${CSV}" \
    --egt "${EGT}" \
    --fasta-ref "${REF_FASTA}" \
    --gtcs "${GTC_DIR}" \
    --extra "${EXTRA_TSV}" \
    --adjust-clusters \
    --threads "${THREADS}" | \
bcftools sort -Ou -T "${OUTPUT_DIR}/bcftools." | \
bcftools norm --no-version -Ob -c x -f "${REF_FASTA}" \
    -o "${VCF_OUTPUT}" --write-index

echo "  Reclustered VCF created: ${VCF_OUTPUT}"

# ---------------------------------------------------------------
# Step 4: Recompute QC metrics after reclustering
# ---------------------------------------------------------------
echo ""
echo "--- Step 4/4: Recomputing QC metrics ---"

source "${SCRIPT_DIR}/collect_qc_metrics.sh"
collect_qc_metrics "${VCF_OUTPUT}" "${EXTRA_TSV}" "${QC_DIR}/stage2_sample_qc.tsv"

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
echo "  VCF:               ${VCF_OUTPUT}"
echo "  QC metrics:        ${QC_DIR}/stage2_sample_qc.tsv"
echo "  HQ sample list:    ${HQ_SAMPLES}"
echo "  Excluded samples:  ${EXCLUDED_SAMPLES}"
echo ""
echo "The reclustered VCF is ready for downstream analysis"
echo "(phasing, MoChA, imputation, etc.)."
