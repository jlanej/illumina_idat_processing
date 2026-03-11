#!/usr/bin/env bash
#
# ancestry_stratified_qc.sh
#
# Perform ancestry-stratified variant QC and PCA using peddy ancestry
# predictions. For each ancestry group meeting a minimum sample count,
# this script:
#   1. Subsets the VCF to samples of that ancestry
#   2. Computes ancestry-specific variant QC (missingness, HWE, MAF)
#   3. Computes ancestry-specific PCA
#   4. Collates all per-ancestry variant QC into a single unified file
#
# Scientific rationale:
#   GWAS best practice recommends performing variant QC within genetically
#   homogeneous groups. Hardy-Weinberg equilibrium (HWE) tests are
#   sensitive to population structure — running HWE on a mixed-ancestry
#   sample inflates deviation due to Wahlund effect, leading to
#   inappropriate exclusion of valid variants. Similarly, allele
#   frequency estimates and missingness patterns can differ across
#   ancestries. Computing PCA within a homogeneous subgroup resolves
#   finer-grained population structure that is masked when PCA is
#   performed on the full multi-ancestry cohort.
#
# References:
#   - Anderson et al. 2010, Nat Protoc 5:1564-73 (ancestry-stratified QC)
#   - Marees et al. 2018, Int J Methods Psychiatr Res 27:e1608
#   - Peterson et al. 2019, AJHG 105:921-935 (within-ancestry PCA)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Perform ancestry-stratified variant QC and PCA using peddy predictions.

Required:
  --vcf FILE              Input VCF/BCF file
  --peddy-het-check FILE  Peddy het_check.csv with ancestry predictions
  --output-dir DIR        Output directory for stratified QC results

Options:
  --min-samples INT       Minimum samples per ancestry group (default: 100)
  --threads INT           Number of threads (default: all available)
  --force                 Re-run even if outputs exist
  --help                  Show this help message

Output files:
  collated_variant_qc.tsv          Unified variant QC across all ancestries
  ancestry_stratified_summary.txt  Summary of stratified analyses
  <ANCESTRY>/variant_qc/           Per-ancestry variant QC
  <ANCESTRY>/pca/                  Per-ancestry PCA projections
EOF
    exit 0
}

VCF=""
PEDDY_HET_CHECK=""
OUTPUT_DIR=""
MIN_SAMPLES=100
THREADS=$(nproc 2>/dev/null || echo 1)
FORCE="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)              VCF="$2"; shift 2 ;;
        --peddy-het-check)  PEDDY_HET_CHECK="$2"; shift 2 ;;
        --output-dir)       OUTPUT_DIR="$2"; shift 2 ;;
        --min-samples)      MIN_SAMPLES="$2"; shift 2 ;;
        --threads)          THREADS="$2"; shift 2 ;;
        --force)            FORCE="true"; shift ;;
        --help)             usage ;;
        *)                  echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "${VCF}" ]]; then
    echo "Error: --vcf is required" >&2
    exit 1
fi
if [[ -z "${PEDDY_HET_CHECK}" ]]; then
    echo "Error: --peddy-het-check is required" >&2
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
if [[ ! -f "${PEDDY_HET_CHECK}" ]]; then
    echo "Error: Peddy het_check.csv not found: ${PEDDY_HET_CHECK}" >&2
    exit 1
fi

# Check required tools
for tool in bcftools; do
    if ! command -v "${tool}" &>/dev/null; then
        echo "Error: ${tool} not found" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}"

SUMMARY="${OUTPUT_DIR}/ancestry_stratified_summary.txt"
COLLATED="${OUTPUT_DIR}/collated_variant_qc.tsv"

# Check idempotency
if [[ "${FORCE}" != "true" && -f "${SUMMARY}" && -f "${COLLATED}" ]]; then
    echo "Ancestry stratified QC outputs already exist. Skipping (use --force to re-run)."
    exit 0
fi

echo "======================================================"
echo "  Ancestry-Stratified QC"
echo "======================================================"
echo ""
echo "  VCF:              ${VCF}"
echo "  Peddy het_check:  ${PEDDY_HET_CHECK}"
echo "  Output:           ${OUTPUT_DIR}"
echo "  Min samples:      ${MIN_SAMPLES}"
echo "  Threads:          ${THREADS}"
echo ""

# -----------------------------------------------------------------
# Step 1: Parse peddy ancestry predictions and identify groups
# -----------------------------------------------------------------
echo "Step 1: Parsing peddy ancestry predictions..."

# Extract sample_id and ancestry-prediction columns from peddy het_check.csv
# peddy het_check.csv is a proper CSV with headers including:
#   sample_id, ancestry-prediction, ancestry-prob, etc.
declare -A ANCESTRY_SAMPLES
declare -A ANCESTRY_COUNTS

# Use awk to parse the CSV and extract ancestry assignments
# We need to find the column indices for sample_id and ancestry-prediction
ANCESTRY_MAP="${OUTPUT_DIR}/ancestry_assignments.tsv"
python3 -c "
import csv, sys
with open('${PEDDY_HET_CHECK}') as f:
    reader = csv.DictReader(f)
    print('sample_id\tancestry')
    for row in reader:
        sid = row.get('sample_id', '')
        anc = row.get('ancestry-prediction', 'UNKNOWN')
        if sid and anc:
            print(f'{sid}\t{anc}')
" > "${ANCESTRY_MAP}"

# Count samples per ancestry
echo "  Ancestry group counts:"
ANCESTRIES=()
while IFS=$'\t' read -r ancestry count; do
    echo "    ${ancestry}: ${count} samples"
    if [[ "${count}" -ge "${MIN_SAMPLES}" ]]; then
        ANCESTRIES+=("${ancestry}")
        echo "      → Meets minimum (${MIN_SAMPLES}), will run stratified QC"
    else
        echo "      → Below minimum (${MIN_SAMPLES}), skipping stratified QC"
    fi
done < <(awk -F'\t' 'NR>1 {print $2}' "${ANCESTRY_MAP}" | sort | uniq -c | sort -rn | awk '{print $2"\t"$1}')

echo ""
echo "  Ancestry groups for stratified analysis: ${ANCESTRIES[*]:-none}"
echo ""

if [[ ${#ANCESTRIES[@]} -eq 0 ]]; then
    echo "No ancestry groups meet the minimum sample count (${MIN_SAMPLES})."
    echo "Skipping ancestry-stratified QC."
    {
        echo "======================================================"
        echo "  Ancestry-Stratified QC Summary"
        echo "======================================================"
        echo ""
        echo "No ancestry groups met the minimum sample threshold (${MIN_SAMPLES})."
        echo "Ancestry-stratified analyses were not performed."
        echo ""
        echo "======================================================"
    } > "${SUMMARY}"
    exit 0
fi

# -----------------------------------------------------------------
# Step 2: For each qualifying ancestry, run stratified analyses
# -----------------------------------------------------------------
PROCESSED_ANCESTRIES=()
for ANCESTRY in "${ANCESTRIES[@]}"; do
    echo "======================================================"
    echo "  Processing ancestry group: ${ANCESTRY}"
    echo "======================================================"
    echo ""

    ANC_DIR="${OUTPUT_DIR}/${ANCESTRY}"
    mkdir -p "${ANC_DIR}"

    # Create sample list for this ancestry
    SAMPLE_LIST="${ANC_DIR}/samples.txt"
    awk -F'\t' -v anc="${ANCESTRY}" 'NR>1 && $2==anc {print $1}' "${ANCESTRY_MAP}" > "${SAMPLE_LIST}"
    N_SAMPLES=$(wc -l < "${SAMPLE_LIST}" | tr -d ' ')
    echo "  Samples: ${N_SAMPLES}"

    # Subset VCF to this ancestry group
    ANC_VCF="${ANC_DIR}/${ANCESTRY}_subset.bcf"
    if [[ "${FORCE}" == "true" || ! -f "${ANC_VCF}" ]]; then
        echo "  Subsetting VCF to ${ANCESTRY} samples..."
        bcftools view -S "${SAMPLE_LIST}" --threads "${THREADS}" \
            -Ob -o "${ANC_VCF}" "${VCF}" 2>/dev/null
        bcftools index --threads "${THREADS}" "${ANC_VCF}" 2>/dev/null
    else
        echo "  ${ANCESTRY} subset VCF already exists."
    fi

    # ----- Ancestry-specific variant QC -----
    ANC_VQC_DIR="${ANC_DIR}/variant_qc"
    if [[ "${FORCE}" == "true" || ! -f "${ANC_VQC_DIR}/variant_qc_summary.txt" ]]; then
        echo "  Running ancestry-specific variant QC for ${ANCESTRY}..."
        bash "${SCRIPT_DIR}/compute_variant_qc.sh" \
            --vcf "${ANC_VCF}" \
            --output-dir "${ANC_VQC_DIR}" \
            --threads "${THREADS}" 2>&1 | sed 's/^/    /'
    else
        echo "  ${ANCESTRY} variant QC already exists."
    fi

    # ----- Ancestry-specific PCA -----
    ANC_PCA_DIR="${ANC_DIR}/pca"
    if [[ "${FORCE}" == "true" || ! -f "${ANC_PCA_DIR}/pca_projections.tsv" ]]; then
        echo "  Running ancestry-specific PCA for ${ANCESTRY}..."
        bash "${SCRIPT_DIR}/ancestry_pca.sh" \
            --vcf "${ANC_VCF}" \
            --output-dir "${ANC_PCA_DIR}" \
            --threads "${THREADS}" 2>&1 | sed 's/^/    /' || true
    else
        echo "  ${ANCESTRY} PCA already exists."
    fi

    PROCESSED_ANCESTRIES+=("${ANCESTRY}")
    echo ""
done

# -----------------------------------------------------------------
# Step 3: Collate variant QC across all ancestries
# -----------------------------------------------------------------
echo "======================================================"
echo "  Collating Variant QC Across Ancestries"
echo "======================================================"
echo ""

# Build arguments for the collation script
COLLATE_ARGS=(
    --output "${COLLATED}"
)

# Add the full-cohort variant QC if available (check stage2 first, then stage1)
FULL_VQC_DIR=""
for vqc_candidate in \
    "$(dirname "${OUTPUT_DIR}")/stage2/qc/variant_qc" \
    "$(dirname "${OUTPUT_DIR}")/stage1/qc/variant_qc"; do
    if [[ -d "${vqc_candidate}" && -f "${vqc_candidate}/variant_qc.vmiss" ]]; then
        FULL_VQC_DIR="${vqc_candidate}"
        break
    fi
done

if [[ -n "${FULL_VQC_DIR}" ]]; then
    COLLATE_ARGS+=(--all-variant-qc-dir "${FULL_VQC_DIR}")
fi

# Add per-ancestry variant QC directories
for ANCESTRY in "${PROCESSED_ANCESTRIES[@]}"; do
    ANC_VQC_DIR="${OUTPUT_DIR}/${ANCESTRY}/variant_qc"
    if [[ -d "${ANC_VQC_DIR}" && -f "${ANC_VQC_DIR}/variant_qc.vmiss" ]]; then
        COLLATE_ARGS+=(--ancestry-variant-qc "${ANCESTRY}:${ANC_VQC_DIR}")
    fi
done

if command -v python3 &>/dev/null; then
    python3 "${SCRIPT_DIR}/collate_variant_qc.py" "${COLLATE_ARGS[@]}" 2>&1 || true
else
    echo "Warning: python3 not found. Skipping variant QC collation." >&2
fi

# -----------------------------------------------------------------
# Step 4: Write summary
# -----------------------------------------------------------------
{
    echo "======================================================"
    echo "  Ancestry-Stratified QC Summary"
    echo "======================================================"
    echo ""
    echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo ""
    echo "Parameters:"
    echo "  Minimum samples per ancestry: ${MIN_SAMPLES}"
    echo "  VCF:                          ${VCF}"
    echo "  Peddy het_check:              ${PEDDY_HET_CHECK}"
    echo ""
    echo "Ancestry Groups Analyzed:"
    for ANCESTRY in "${PROCESSED_ANCESTRIES[@]}"; do
        N_SAMPLES=$(wc -l < "${OUTPUT_DIR}/${ANCESTRY}/samples.txt" | tr -d ' ')
        echo "  ${ANCESTRY}: ${N_SAMPLES} samples"
        ANC_VQC="${OUTPUT_DIR}/${ANCESTRY}/variant_qc/variant_qc_summary.txt"
        if [[ -f "${ANC_VQC}" ]]; then
            echo "    Variant QC: available"
        fi
        ANC_PCA="${OUTPUT_DIR}/${ANCESTRY}/pca/pca_projections.tsv"
        if [[ -f "${ANC_PCA}" ]]; then
            echo "    PCA:        available"
        fi
    done
    echo ""
    echo "Collated variant QC: ${COLLATED}"
    echo ""
    echo "Scientific rationale:"
    echo "  Ancestry-stratified QC is performed because Hardy-Weinberg"
    echo "  equilibrium tests assume panmixia. Running HWE on a mixed-"
    echo "  ancestry cohort inflates deviation due to the Wahlund effect,"
    echo "  which can cause inappropriate exclusion of valid variants."
    echo "  Similarly, allele frequencies and genotyping performance may"
    echo "  differ across ancestries. Within-ancestry PCA resolves finer"
    echo "  population structure that is masked in multi-ancestry PCA."
    echo ""
    echo "References:"
    echo "  Anderson et al. 2010, Nat Protoc 5:1564-73"
    echo "  Marees et al. 2018, Int J Methods Psychiatr Res 27:e1608"
    echo "  Peterson et al. 2019, AJHG 105:921-935"
    echo ""
    echo "======================================================"
} | tee "${SUMMARY}"

echo ""
echo "Ancestry-stratified QC complete."
