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
  --exclude-samples FILE  File of sample IDs (one per line) to exclude from
                          HWE testing and PCA training sets.  These samples
                          are removed prior to variant QC and within-ancestry
                          PCA.  Used for related and het-outlier samples.
  --tstv-file FILE        bcftools stats tstv output file for project-level
                          Ti/Tv ratio in collated output
  --sample-qc FILE        Sample QC TSV with computed_gender column for
                          sex-chromosome variant QC (chrX HWE on females,
                          chrY on males)
  --genome STRING         Genome build for PAR/XTR exclusion in sex-chr QC
                          (CHM13, GRCh38, or GRCh37; default: CHM13)
  --threads INT           Number of threads (default: all available)
  --force                 Re-run even if outputs exist
  --help                  Show this help message

Output files:
  collated_variant_qc.tsv          Unified variant QC across all ancestries
  hwe_passing_variants.txt         Variant IDs passing HWE in all ancestry groups
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
EXCLUDE_SAMPLES=""
TSTV_FILE=""
SAMPLE_QC=""
GENOME="CHM13"
THREADS=$(nproc 2>/dev/null || echo 1)
FORCE="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)              VCF="$2"; shift 2 ;;
        --peddy-het-check)  PEDDY_HET_CHECK="$2"; shift 2 ;;
        --output-dir)       OUTPUT_DIR="$2"; shift 2 ;;
        --min-samples)      MIN_SAMPLES="$2"; shift 2 ;;
        --exclude-samples)  EXCLUDE_SAMPLES="$2"; shift 2 ;;
        --tstv-file)        TSTV_FILE="$2"; shift 2 ;;
        --sample-qc)        SAMPLE_QC="$2"; shift 2 ;;
        --genome)           GENOME="$2"; shift 2 ;;
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
if ! command -v bcftools &>/dev/null; then
    echo "Error: bcftools not found" >&2
    exit 1
fi

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
if [[ -n "${EXCLUDE_SAMPLES}" && -s "${EXCLUDE_SAMPLES}" ]]; then
    N_PRE_EXCLUDE=$(wc -l < "${EXCLUDE_SAMPLES}" | tr -d ' ')
    echo "  Pre-QC excluded:  ${N_PRE_EXCLUDE} samples (${EXCLUDE_SAMPLES})"
fi
echo ""

# -----------------------------------------------------------------
# Step 1: Parse peddy ancestry predictions and identify groups
# -----------------------------------------------------------------
echo "Step 1: Parsing peddy ancestry predictions..."

# Extract sample_id and ancestry-prediction columns from peddy het_check.csv
# peddy het_check.csv is a proper CSV with headers including:
#   sample_id, ancestry-prediction, ancestry-prob, etc.
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
    echo "Per-ancestry analyses will be skipped."
    echo "Will still collate full-cohort variant QC (if available) for report plots."
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

    # Create unrelated sample list (excluding pre-QC samples)
    UNRELATED_SAMPLE_LIST="${ANC_DIR}/unrelated_samples.txt"
    if [[ -n "${EXCLUDE_SAMPLES}" && -s "${EXCLUDE_SAMPLES}" ]]; then
        grep -vxFf "${EXCLUDE_SAMPLES}" "${SAMPLE_LIST}" > "${UNRELATED_SAMPLE_LIST}" || true
        N_UNRELATED=$(wc -l < "${UNRELATED_SAMPLE_LIST}" | tr -d ' ')
        echo "  Unrelated samples (after exclusions): ${N_UNRELATED}"
    else
        cp "${SAMPLE_LIST}" "${UNRELATED_SAMPLE_LIST}"
        N_UNRELATED="${N_SAMPLES}"
    fi

    # Skip this ancestry if all samples were excluded
    if [[ "${N_UNRELATED}" -eq 0 ]]; then
        echo "  WARNING: All ${N_SAMPLES} samples in ${ANCESTRY} were excluded."
        echo "           Skipping variant QC and PCA for this group."
        echo ""
        continue
    fi

    # Subset VCF to unrelated samples of this ancestry for variant QC
    ANC_VCF="${ANC_DIR}/${ANCESTRY}_subset.bcf"
    if [[ "${FORCE}" == "true" || ! -f "${ANC_VCF}" ]]; then
        echo "  Subsetting VCF to ${ANCESTRY} unrelated samples..."
        bcftools view -S "${UNRELATED_SAMPLE_LIST}" --threads "${THREADS}" \
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
    # PCA training is on unrelated samples; all ancestry samples are projected
    ANC_PCA_DIR="${ANC_DIR}/pca"
    ANC_PCA_VCF="${ANC_DIR}/${ANCESTRY}_all_subset.bcf"
    if [[ "${FORCE}" == "true" || ! -f "${ANC_PCA_DIR}/pca_projections.tsv" ]]; then
        # Create full-ancestry VCF (including excluded samples) for projection
        if [[ -n "${EXCLUDE_SAMPLES}" && -s "${EXCLUDE_SAMPLES}" ]]; then
            bcftools view -S "${SAMPLE_LIST}" --threads "${THREADS}" \
                -Ob -o "${ANC_PCA_VCF}" "${VCF}" 2>/dev/null
            bcftools index --threads "${THREADS}" "${ANC_PCA_VCF}" 2>/dev/null
        else
            ANC_PCA_VCF="${ANC_VCF}"
        fi
        echo "  Running ancestry-specific PCA for ${ANCESTRY}..."
        ANC_PCA_ARGS=(
            --vcf "${ANC_PCA_VCF}"
            --output-dir "${ANC_PCA_DIR}"
            --threads "${THREADS}"
        )
        if [[ -n "${EXCLUDE_SAMPLES}" && -s "${EXCLUDE_SAMPLES}" ]]; then
            ANC_PCA_ARGS+=(--exclude-samples "${EXCLUDE_SAMPLES}")
        fi
        bash "${SCRIPT_DIR}/ancestry_pca.sh" "${ANC_PCA_ARGS[@]}" 2>&1 | sed 's/^/    /' || true
    else
        echo "  ${ANCESTRY} PCA already exists."
    fi

    PROCESSED_ANCESTRIES+=("${ANCESTRY}")
    echo ""
done

# -----------------------------------------------------------------
# Step 2b: Sex-chromosome variant QC (chrX and chrY)
# -----------------------------------------------------------------
SEX_CHR_QC_DIR="${OUTPUT_DIR}/sex_chr_qc"

echo "======================================================"
echo "  Sex-Chromosome Variant QC (chrX/chrY)"
echo "======================================================"
echo ""

if [[ -n "${SAMPLE_QC}" && -f "${SAMPLE_QC}" ]]; then
    if [[ "${FORCE}" == "true" || ! -d "${SEX_CHR_QC_DIR}/chrX" ]]; then
        echo "  Running sex-chromosome variant QC..."
        bash "${SCRIPT_DIR}/compute_sex_chr_variant_qc.sh" \
            --vcf "${VCF}" \
            --sample-qc "${SAMPLE_QC}" \
            --output-dir "${SEX_CHR_QC_DIR}" \
            --genome "${GENOME}" \
            --threads "${THREADS}" 2>&1 | sed 's/^/    /'
    else
        echo "  Sex-chromosome variant QC already exists."
    fi
else
    echo "  WARNING: Skipping sex-chromosome variant QC." >&2
    echo "    --sample-qc was not provided or file does not exist." >&2
    echo "    chrX and chrY variants will be absent from collated_variant_qc.tsv." >&2
    echo "    To include sex-chr QC, re-run with --sample-qc <sample_qc.tsv>." >&2
fi
echo ""

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

# Include Ti/Tv stats: prefer explicit --tstv-file arg, fallback to
# standard location relative to variant QC directory
TSTV_RESOLVED=""
if [[ -n "${TSTV_FILE}" && -f "${TSTV_FILE}" ]]; then
    TSTV_RESOLVED="${TSTV_FILE}"
elif [[ -n "${FULL_VQC_DIR}" ]]; then
    TSTV_CANDIDATE="${FULL_VQC_DIR}/../tstv_stats.txt"
    if [[ -f "${TSTV_CANDIDATE}" ]]; then
        TSTV_RESOLVED="${TSTV_CANDIDATE}"
    fi
fi
if [[ -n "${TSTV_RESOLVED}" ]]; then
    COLLATE_ARGS+=(--tstv-file "${TSTV_RESOLVED}")
fi

# Add per-ancestry variant QC directories
for ANCESTRY in "${PROCESSED_ANCESTRIES[@]}"; do
    ANC_VQC_DIR="${OUTPUT_DIR}/${ANCESTRY}/variant_qc"
    if [[ -d "${ANC_VQC_DIR}" && -f "${ANC_VQC_DIR}/variant_qc.vmiss" ]]; then
        COLLATE_ARGS+=(--ancestry-variant-qc "${ANCESTRY}:${ANC_VQC_DIR}")
    fi
done

# Add sex-chromosome QC directory if available
if [[ -d "${SEX_CHR_QC_DIR}" ]]; then
    COLLATE_ARGS+=(--sex-chr-qc-dir "${SEX_CHR_QC_DIR}")
fi

if command -v python3 &>/dev/null; then
    python3 "${SCRIPT_DIR}/collate_variant_qc.py" "${COLLATE_ARGS[@]}" 2>&1 || true
else
    echo "Warning: python3 not found. Skipping variant QC collation." >&2
fi

# -----------------------------------------------------------------
# Step 3b: Build HWE-passing variant list (intersection across ancestries)
# -----------------------------------------------------------------
HWE_PASSING="${OUTPUT_DIR}/hwe_passing_variants.txt"
HWE_P_THRESH="1e-6"

echo "======================================================"
echo "  Building HWE-Passing Variant List"
echo "======================================================"
echo ""

if [[ ${#PROCESSED_ANCESTRIES[@]} -gt 0 ]]; then
    # For each ancestry, extract variant IDs where HWE p >= threshold
    HWE_FIRST=""
    for ANCESTRY in "${PROCESSED_ANCESTRIES[@]}"; do
        HARDY_FILE="${OUTPUT_DIR}/${ANCESTRY}/variant_qc/variant_qc.hardy"
        ANC_HWE_PASS="${OUTPUT_DIR}/${ANCESTRY}/variant_qc/hwe_passing_ids.txt"
        if [[ -f "${HARDY_FILE}" ]]; then
            # Extract ID column for variants with HWE p >= threshold
            awk -v thresh="${HWE_P_THRESH}" \
                'NR>1 { p=$NF+0; if (p >= thresh) print $2 }' \
                "${HARDY_FILE}" | sort > "${ANC_HWE_PASS}"
            N_PASS=$(wc -l < "${ANC_HWE_PASS}" | tr -d ' ')
            echo "  ${ANCESTRY}: ${N_PASS} variants pass HWE (p >= ${HWE_P_THRESH})"
            if [[ -z "${HWE_FIRST}" ]]; then
                HWE_FIRST="${ANC_HWE_PASS}"
            fi
        else
            echo "  ${ANCESTRY}: No .hardy file found — skipping"
        fi
    done

    # Intersect across all ancestries: only keep variants passing HWE in ALL groups
    if [[ -n "${HWE_FIRST}" ]]; then
        cp "${HWE_FIRST}" "${HWE_PASSING}.tmp"
        for ANCESTRY in "${PROCESSED_ANCESTRIES[@]}"; do
            ANC_HWE_PASS="${OUTPUT_DIR}/${ANCESTRY}/variant_qc/hwe_passing_ids.txt"
            if [[ -f "${ANC_HWE_PASS}" ]]; then
                comm -12 "${HWE_PASSING}.tmp" "${ANC_HWE_PASS}" > "${HWE_PASSING}.tmp2"
                mv "${HWE_PASSING}.tmp2" "${HWE_PASSING}.tmp"
            fi
        done
        mv "${HWE_PASSING}.tmp" "${HWE_PASSING}"
        N_HWE_PASS=$(wc -l < "${HWE_PASSING}" | tr -d ' ')
        echo ""
        echo "  HWE-passing variants (intersection across all ancestries): ${N_HWE_PASS}"
    else
        touch "${HWE_PASSING}"
        echo "  No .hardy files found — HWE-passing list is empty"
    fi
else
    touch "${HWE_PASSING}"
    echo "  No ancestry groups processed — HWE-passing list is empty"
fi
echo ""

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
    if [[ ${#PROCESSED_ANCESTRIES[@]} -eq 0 ]]; then
        echo "  none (no groups met minimum sample threshold ${MIN_SAMPLES})"
    else
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
    fi
    echo ""
    echo "Collated variant QC: ${COLLATED}"
    echo ""
    if [[ -f "${HWE_PASSING}" ]]; then
        N_HWE_FINAL=$(wc -l < "${HWE_PASSING}" | tr -d ' ')
        echo "HWE-passing variant list: ${HWE_PASSING} (${N_HWE_FINAL} variants)"
        echo "  These variants pass HWE (p >= ${HWE_P_THRESH}) across ALL ancestry groups."
        echo "  Used as the variant set for full-cohort PCA (Price et al. 2006)."
    fi
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
