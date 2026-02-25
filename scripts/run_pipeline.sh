#!/usr/bin/env bash
#
# run_pipeline.sh
#
# Main entry point for the Illumina IDAT processing pipeline.
# Processes raw Illumina array IDAT data in two stages:
#
#   Stage 1: Initial genotyping - IDAT → GTC → VCF with QC metrics
#   Stage 2: Reclustering - filter high-quality samples, recluster, reprocess
#
# Designed to require minimal arguments and automate manifest/resource
# acquisition where possible.
#
# Inspired by mocha (https://github.com/freeseek/mocha) and mochawdl
# (https://github.com/freeseek/mochawdl).
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Defaults
GENOME="GRCh38"
THREADS=1
MIN_CALL_RATE=0.97
MAX_LRR_SD=0.35
SKIP_STAGE2="false"
SKIP_DOWNLOAD="false"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Process raw Illumina IDAT files to genotyped VCF with two-stage QC.

Required (minimal):
  --idat-dir DIR         Directory containing IDAT files
  --output-dir DIR       Output directory

Manifest options (one of the following):
  --array-name NAME      Known array name (auto-downloads manifests)
  --manifest-dir DIR     Directory containing BPM, EGT, and CSV files
  --bpm FILE             BPM manifest file (with --egt and --csv)
  --egt FILE             EGT cluster file (with --bpm and --csv)
  --csv FILE             CSV manifest file (with --bpm and --egt)

Reference options:
  --ref-dir DIR          Directory with reference genome (auto-downloads if missing)
  --ref-fasta FILE       Reference FASTA file (overrides --ref-dir)
  --genome NAME          Genome build: GRCh37 or GRCh38 (default: ${GENOME})

Processing options:
  --threads INT          Number of threads (default: ${THREADS})
  --min-call-rate FLOAT  Min call rate for Stage 2 HQ samples (default: ${MIN_CALL_RATE})
  --max-lrr-sd FLOAT     Max LRR SD for Stage 2 HQ samples (default: ${MAX_LRR_SD})
  --skip-stage2          Skip Stage 2 (reclustering)
  --skip-download        Do not auto-download manifests or reference

  --help                 Show this help message

Examples:
  # Minimal: auto-download manifests and reference for a known array
  $(basename "$0") --idat-dir /data/idats --array-name GSA-24v3-0_A1 --output-dir /data/output

  # With pre-downloaded manifests and reference
  $(basename "$0") --idat-dir /data/idats --manifest-dir /data/manifests \\
      --ref-fasta /data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \\
      --output-dir /data/output

  # Provide manifest files directly
  $(basename "$0") --idat-dir /data/idats --bpm array.bpm --egt array.egt \\
      --csv array.csv --ref-dir /data/ref --output-dir /data/output
EOF
    exit 0
}

IDAT_DIR=""
OUTPUT_DIR=""
ARRAY_NAME=""
MANIFEST_DIR=""
BPM=""
EGT=""
CSV=""
REF_DIR=""
REF_FASTA=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --idat-dir)       IDAT_DIR="$2"; shift 2 ;;
        --output-dir)     OUTPUT_DIR="$2"; shift 2 ;;
        --array-name)     ARRAY_NAME="$2"; shift 2 ;;
        --manifest-dir)   MANIFEST_DIR="$2"; shift 2 ;;
        --bpm)            BPM="$2"; shift 2 ;;
        --egt)            EGT="$2"; shift 2 ;;
        --csv)            CSV="$2"; shift 2 ;;
        --ref-dir)        REF_DIR="$2"; shift 2 ;;
        --ref-fasta)      REF_FASTA="$2"; shift 2 ;;
        --genome)         GENOME="$2"; shift 2 ;;
        --threads)        THREADS="$2"; shift 2 ;;
        --min-call-rate)  MIN_CALL_RATE="$2"; shift 2 ;;
        --max-lrr-sd)     MAX_LRR_SD="$2"; shift 2 ;;
        --skip-stage2)    SKIP_STAGE2="true"; shift ;;
        --skip-download)  SKIP_DOWNLOAD="true"; shift ;;
        --help)           usage ;;
        *)                echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "${IDAT_DIR}" ]]; then
    echo "Error: --idat-dir is required" >&2
    exit 1
fi
if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "Error: --output-dir is required" >&2
    exit 1
fi
if [[ ! -d "${IDAT_DIR}" ]]; then
    echo "Error: IDAT directory not found: ${IDAT_DIR}" >&2
    exit 1
fi

# Check bcftools is available
if ! command -v bcftools &>/dev/null; then
    echo "Error: bcftools not found. Run scripts/install_dependencies.sh first." >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

echo "======================================================"
echo "  Illumina IDAT Processing Pipeline"
echo "======================================================"
echo ""

# ---------------------------------------------------------------
# Resolve manifest files
# ---------------------------------------------------------------
if [[ -n "${MANIFEST_DIR}" ]]; then
    # Auto-detect BPM, EGT, CSV from manifest directory
    [[ -z "${BPM}" ]] && BPM=$(find "${MANIFEST_DIR}" -name "*.bpm" -print -quit 2>/dev/null)
    [[ -z "${EGT}" ]] && EGT=$(find "${MANIFEST_DIR}" -name "*.egt" -print -quit 2>/dev/null)
    [[ -z "${CSV}" ]] && CSV=$(find "${MANIFEST_DIR}" -name "*.csv" -print -quit 2>/dev/null)
fi

if [[ -z "${BPM}" || -z "${EGT}" || -z "${CSV}" ]]; then
    if [[ "${SKIP_DOWNLOAD}" == "true" ]]; then
        echo "Error: Manifest files not found and --skip-download is set." >&2
        echo "Provide --bpm, --egt, --csv or --manifest-dir." >&2
        exit 1
    fi

    MANIFEST_DIR="${OUTPUT_DIR}/manifests"
    echo "Downloading manifest files..."
    DOWNLOAD_ARGS=(--output-dir "${MANIFEST_DIR}")

    if [[ -n "${ARRAY_NAME}" ]]; then
        DOWNLOAD_ARGS+=(--array-name "${ARRAY_NAME}")
    else
        echo "  Attempting to auto-detect array type from IDAT files..."
        DOWNLOAD_ARGS+=(--detect-from "${IDAT_DIR}")
    fi

    bash "${SCRIPT_DIR}/download_manifests.sh" "${DOWNLOAD_ARGS[@]}"
    echo ""

    [[ -z "${BPM}" ]] && BPM=$(find "${MANIFEST_DIR}" -name "*.bpm" -print -quit 2>/dev/null)
    [[ -z "${EGT}" ]] && EGT=$(find "${MANIFEST_DIR}" -name "*.egt" -print -quit 2>/dev/null)
    [[ -z "${CSV}" ]] && CSV=$(find "${MANIFEST_DIR}" -name "*.csv" -print -quit 2>/dev/null)

    if [[ -z "${BPM}" || -z "${EGT}" || -z "${CSV}" ]]; then
        echo "Error: Could not find all manifest files after download." >&2
        echo "Please provide --bpm, --egt, and --csv directly." >&2
        exit 1
    fi
fi

# ---------------------------------------------------------------
# Resolve reference genome
# ---------------------------------------------------------------
if [[ -z "${REF_FASTA}" ]]; then
    if [[ -z "${REF_DIR}" ]]; then
        REF_DIR="${OUTPUT_DIR}/reference"
    fi

    # Try to find existing reference
    if [[ "${GENOME}" == "GRCh38" ]]; then
        REF_FASTA="${REF_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    else
        REF_FASTA="${REF_DIR}/human_g1k_v37.fasta"
    fi

    if [[ ! -f "${REF_FASTA}" ]]; then
        if [[ "${SKIP_DOWNLOAD}" == "true" ]]; then
            echo "Error: Reference genome not found and --skip-download is set." >&2
            exit 1
        fi
        echo "Downloading reference genome (${GENOME})..."
        bash "${SCRIPT_DIR}/download_reference.sh" \
            --genome "${GENOME}" \
            --output-dir "${REF_DIR}"
        echo ""
    fi
fi

if [[ ! -f "${REF_FASTA}" ]]; then
    echo "Error: Reference FASTA not found: ${REF_FASTA}" >&2
    exit 1
fi

# ---------------------------------------------------------------
# Print configuration
# ---------------------------------------------------------------
echo "Configuration:"
echo "  IDAT dir:       ${IDAT_DIR}"
echo "  Output dir:     ${OUTPUT_DIR}"
echo "  BPM:            ${BPM}"
echo "  EGT:            ${EGT}"
echo "  CSV:            ${CSV}"
echo "  Reference:      ${REF_FASTA}"
echo "  Genome:         ${GENOME}"
echo "  Threads:        ${THREADS}"
echo "  Min call rate:  ${MIN_CALL_RATE}"
echo "  Max LRR SD:     ${MAX_LRR_SD}"
echo ""

# ---------------------------------------------------------------
# Stage 1: Initial genotyping
# ---------------------------------------------------------------
STAGE1_DIR="${OUTPUT_DIR}/stage1"

echo "======================================================"
echo "  Running Stage 1: Initial Genotyping"
echo "======================================================"
echo ""

bash "${SCRIPT_DIR}/stage1_initial_genotyping.sh" \
    --idat-dir "${IDAT_DIR}" \
    --bpm "${BPM}" \
    --egt "${EGT}" \
    --csv "${CSV}" \
    --ref-fasta "${REF_FASTA}" \
    --output-dir "${STAGE1_DIR}" \
    --threads "${THREADS}"

echo ""

# ---------------------------------------------------------------
# Stage 2: Recluster and reprocess
# ---------------------------------------------------------------
if [[ "${SKIP_STAGE2}" == "true" ]]; then
    echo "Skipping Stage 2 (--skip-stage2 specified)"
    FINAL_VCF="${STAGE1_DIR}/vcf/stage1_initial.bcf"
    FINAL_QC="${STAGE1_DIR}/qc/stage1_sample_qc.tsv"
else
    STAGE2_DIR="${OUTPUT_DIR}/stage2"

    echo "======================================================"
    echo "  Running Stage 2: Recluster and Reprocess"
    echo "======================================================"
    echo ""

    bash "${SCRIPT_DIR}/stage2_recluster.sh" \
        --idat-dir "${IDAT_DIR}" \
        --bpm "${BPM}" \
        --egt "${EGT}" \
        --csv "${CSV}" \
        --ref-fasta "${REF_FASTA}" \
        --stage1-dir "${STAGE1_DIR}" \
        --output-dir "${STAGE2_DIR}" \
        --min-call-rate "${MIN_CALL_RATE}" \
        --max-lrr-sd "${MAX_LRR_SD}" \
        --threads "${THREADS}"

    FINAL_VCF="${STAGE2_DIR}/vcf/stage2_reclustered.bcf"
    FINAL_QC="${STAGE2_DIR}/qc/stage2_sample_qc.tsv"
fi

echo ""
echo "======================================================"
echo "  Pipeline Complete"
echo "======================================================"
echo ""
echo "Final outputs:"
echo "  VCF:        ${FINAL_VCF}"
echo "  QC metrics: ${FINAL_QC}"
echo ""
echo "The VCF contains genotypes with BAF and LRR intensities"
echo "ready for downstream analysis (phasing, MoChA, imputation)."
