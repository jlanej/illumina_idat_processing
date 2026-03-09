#!/usr/bin/env bash
#
# run_peddy.sh
#
# Run peddy for pedigree/sex/ancestry QC on genotyped VCF output.
#
# Peddy validates sex, ancestry, and relatedness from VCF + PED input.
# If no pedigree file is provided, a minimal PED is generated (all
# samples as unrelated singletons) so that peddy can still predict
# sex, ancestry, and discover relationships.
#
# Genome-aware VCF preparation:
#   When the pipeline genome is not GRCh38, the stage2 BCF is lifted
#   over to GRCh38 coordinates via bcftools +liftover and filtered to
#   peddy GRCH38 site coordinate windows. This produces a small,
#   correctly-coordinated VCF that peddy can use regardless of the
#   source reference genome.
#
# After running, a final pedigree incorporating peddy's discovered
# relationships is written.
#
# Idempotent: skips execution when expected outputs already exist
# unless --force is specified.
#
set -euo pipefail

# --- Resource URLs ---
CHAIN_URL_CHM13="https://hgdownload.soe.ucsc.edu/gbdb/hs1/liftOver/chm13v2-hg38.over.chain.gz"
CHAIN_URL_GRCh37="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
PEDDY_SITES_URL="https://raw.githubusercontent.com/brentp/peddy/master/peddy/GRCH38.sites"
GRCH38_FASTA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
GRCH38_FASTA_NAME="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
PEDDY_MIN_OVERLAP_WARN_COUNT=500
PEDDY_MIN_OVERLAP_WARN_SITE_FRACTION=0.005
PEDDY_MIN_LIFTOVER_RETAIN_FRACTION=0.10
PEDDY_COORD_BUFFER_BP=100
PEDDY_COUNT_WAIT_ATTEMPTS=50

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Run peddy pedigree/sex/ancestry QC on pipeline VCF output.

Required:
  --vcf FILE             Input VCF/BCF file
  --output-dir DIR       Output directory for peddy results

Optional:
  --genome NAME          Genome build: CHM13, GRCh37, or GRCh38 (default: GRCh38)
                         When not GRCh38, coordinates are lifted over via
                         bcftools +liftover so peddy can match its GRCH38 sites.
  --src-fasta FILE       Source reference FASTA (required when genome != GRCh38)
  --ref-dir DIR          Reference directory (for caching GRCh38 target FASTA)
  --ped-file FILE        Pedigree file (.ped/.fam). If not provided, a minimal
                         PED with unrelated singletons is generated.
  --sample-qc FILE       Sample QC TSV for sex annotation in generated PED
  --threads INT          Number of threads (default: 4)
  --debug                Print expanded bcftools commands (default: enabled)
  --no-debug             Disable debug command logging
  --force                Re-run even if outputs exist
  --help                 Show this help message
EOF
    exit 0
}

VCF=""
OUTPUT_DIR=""
GENOME="GRCh38"
SRC_FASTA=""
REF_DIR=""
PED_FILE=""
SAMPLE_QC=""
THREADS=4
DEBUG="true"
FORCE="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)         VCF="$2"; shift 2 ;;
        --output-dir)  OUTPUT_DIR="$2"; shift 2 ;;
        --genome)      GENOME="$2"; shift 2 ;;
        --src-fasta)   SRC_FASTA="$2"; shift 2 ;;
        --ref-dir)     REF_DIR="$2"; shift 2 ;;
        --ped-file)    PED_FILE="$2"; shift 2 ;;
        --sample-qc)   SAMPLE_QC="$2"; shift 2 ;;
        --threads)     THREADS="$2"; shift 2 ;;
        --debug)       DEBUG="true"; shift ;;
        --no-debug)    DEBUG="false"; shift ;;
        --force)       FORCE="true"; shift ;;
        --help)        usage ;;
        *)             echo "Error: Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Validate required arguments
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

mkdir -p "${OUTPUT_DIR}"

PREFIX="${OUTPUT_DIR}/peddy"

# --- Idempotency check ---
EXPECTED_OUTPUTS=(
    "${PREFIX}.het_check.csv"
    "${PREFIX}.sex_check.csv"
    "${PREFIX}.ped_check.csv"
    "${PREFIX}.peddy.ped"
)

if [[ "${FORCE}" != "true" ]]; then
    all_exist=true
    for f in "${EXPECTED_OUTPUTS[@]}"; do
        if [[ ! -f "$f" ]]; then
            all_exist=false
            break
        fi
    done
    if [[ "${all_exist}" == "true" ]]; then
        echo "Peddy outputs already exist in ${OUTPUT_DIR}. Skipping (use --force to re-run)."
        exit 0
    fi
fi

# --- Prepare working directories ---
TMP_DIR="${OUTPUT_DIR}/tmp/run_peddy"
RESOURCE_DIR="${OUTPUT_DIR}/resources"
mkdir -p "${TMP_DIR}" "${RESOURCE_DIR}"
if [[ -n "${REF_DIR}" ]]; then
    mkdir -p "${REF_DIR}"
fi

# ======================================================================
#  Helper functions
# ======================================================================

# Download a file with wget or curl, with caching
_cached_download() {
    local url="$1" dest="$2" label="$3"
    if [[ -f "${dest}" ]]; then
        return 0
    fi
    echo "  Downloading ${label}..."
    if wget -q "${url}" -O "${dest}" 2>/dev/null || \
       curl -sfL "${url}" -o "${dest}" 2>/dev/null; then
        echo "  Cached: ${dest}"
    else
        echo "Error: Failed to download ${label} from ${url}" >&2
        rm -f "${dest}"
        return 1
    fi
}

# Ensure BCF/VCF is indexed
_ensure_index() {
    local vcf="$1"
    if [[ "${vcf}" == *.bcf && ! -f "${vcf}.csi" ]]; then
        bcftools index "${vcf}" 2>/dev/null || true
    elif [[ "${vcf}" == *.vcf.gz && ! -f "${vcf}.tbi" && ! -f "${vcf}.csi" ]]; then
        bcftools index -t "${vcf}" 2>/dev/null || true
    fi
}

# Create bare -> chr chromosome rename map
_create_add_chr_map() {
    local dest="$1"
    : > "${dest}"
    for i in $(seq 1 22); do echo "${i} chr${i}"; done >> "${dest}"
    echo "X chrX" >> "${dest}"
    echo "Y chrY" >> "${dest}"
    echo "MT chrM" >> "${dest}"
}

# Detect whether VCF chromosomes use "chr" prefix
_has_chr_prefix() {
    local vcf="$1"
    local first_chr
    first_chr=$(bcftools view -H "${vcf}" 2>/dev/null | head -1 | cut -f1)
    [[ "${first_chr}" == chr* ]]
}

_count_variants() {
    local vcf="$1"
    local indexed_count_raw
    local indexed_count_trimmed
    if indexed_count_raw=$(bcftools index -n "${vcf}" 2>/dev/null); then
        indexed_count_trimmed=$(tr -d '[:space:]' <<< "${indexed_count_raw}")
        if [[ "${indexed_count_trimmed}" =~ ^[0-9]+$ ]]; then
            echo "${indexed_count_trimmed}"
            return
        elif [[ "${DEBUG}" == "true" ]]; then
            echo "  [debug] bcftools index -n returned non-numeric count for ${vcf}: ${indexed_count_raw}" >&2
        fi
    elif [[ "${DEBUG}" == "true" ]]; then
        echo "  [debug] bcftools index -n failed for ${vcf}; falling back to streamed count" >&2
    fi
    bcftools view -H "${vcf}" 2>/dev/null | wc -l | tr -d ' '
}

_debug_log_command() {
    local cmd="$1"
    if [[ "${DEBUG}" == "true" ]]; then
        echo "  [debug] ${cmd}" >&2
    fi
}

_percent() {
    local numerator="$1"
    local denominator="$2"
    awk -v n="${numerator}" -v d="${denominator}" \
        'BEGIN{ if (d > 0) printf "%.2f", (100*n/d); else printf "0.00" }'
}

_report_peddy_overlap() {
    local context="$1"
    local candidate_label="$2"
    local candidate_count="$3"
    local sites_count="$4"
    local matched_count="$5"

    local candidate_pct sites_pct
    candidate_pct=$(_percent "${matched_count}" "${candidate_count}")
    sites_pct=$(_percent "${matched_count}" "${sites_count}")
    echo "  ${context}: ${matched_count} variants matched peddy sites (${candidate_pct}% of ${candidate_label}, ${sites_pct}% of GRCH38.sites)"

    local warn_overlap="false"
    if (( matched_count < PEDDY_MIN_OVERLAP_WARN_COUNT )); then
        warn_overlap="true"
    elif awk -v n="${matched_count}" -v d="${sites_count}" -v thr="${PEDDY_MIN_OVERLAP_WARN_SITE_FRACTION}" \
        'BEGIN{ exit !(d > 0 && (n/d) < thr) }'; then
        warn_overlap="true"
    fi

    if [[ "${warn_overlap}" == "true" ]]; then
        echo "Warning: Very low overlap with peddy GRCh38 sites (${matched_count} matched). Check genome build, liftover chain/source FASTA compatibility, and chromosome naming."
        if [[ "${GENOME}" == "GRCh38" ]]; then
            echo "Hint: For legacy GRCh37 manifests (e.g., 1000G Omni2.5), verify manifest realignment to GRCh38 completed before genotyping. If coordinates are still GRCh37, run with --genome GRCh37 and provide --src-fasta so peddy can liftover."
        fi
    fi
}

_prepare_peddy_site_windows() {
    local sites_file="$1"
    local windows_file="$2"
    local buffer_bp="$3"
    # peddy GRCH38.sites format is CHROM:POS:REF:ALT (one record per line).
    # End coordinates beyond chromosome bounds are tolerated by bcftools region
    # handling and are effectively clipped during querying.
    awk -F: -v b="${buffer_bp}" 'NF >= 2 {chr=$1; if (chr !~ /^chr/) chr="chr"chr; start=$2-b; if (start < 1) start=1; end=$2+b; print chr"\t"start"\t"end}' \
        "${sites_file}" > "${windows_file}"
}

# ======================================================================
#  Prepare peddy-compatible VCF
# ======================================================================
# Peddy requires a VCF.gz with tabix index.
#
# Strategy (following user specification):
#   bcftools +liftover  (source -> GRCh38, when genome != GRCh38)
#     | bcftools view -T GRCH38.sites.windows  (coordinate ±buffer subset)
#     | bcftools sort -Oz -o peddy_input.vcf.gz
# ======================================================================

_prepare_peddy_input() {
    local src_vcf="$1"

    _ensure_index "${src_vcf}"

    # Download peddy GRCH38 sites (used for filtering in all paths)
    _cached_download "${PEDDY_SITES_URL}" \
        "${RESOURCE_DIR}/GRCH38.sites" "peddy GRCH38 sites"
    _prepare_peddy_site_windows \
        "${RESOURCE_DIR}/GRCH38.sites" \
        "${RESOURCE_DIR}/GRCH38.sites.windows" \
        "${PEDDY_COORD_BUFFER_BP}"

    local sites_count
    sites_count=$(wc -l < "${RESOURCE_DIR}/GRCH38.sites" | tr -d ' ')
    echo "  Peddy GRCH38 sites: ${sites_count}"
    echo "  Coordinate match window: +/-${PEDDY_COORD_BUFFER_BP} bp"

    local source_variant_count
    source_variant_count=$(_count_variants "${src_vcf}")
    echo "  Source variants: ${source_variant_count}"

    if [[ "${GENOME}" == "GRCh38" ]]; then
        echo "Genome is GRCh38 — subsetting VCF to peddy sites..."
        _prepare_grch38_subset "${src_vcf}" "${source_variant_count}" "${sites_count}"
    elif [[ "${GENOME}" == "CHM13" || "${GENOME}" == "GRCh37" ]]; then
        echo "Genome is ${GENOME} — lifting over to GRCh38 via bcftools +liftover..."
        _prepare_liftover_subset "${src_vcf}" "${source_variant_count}" "${sites_count}"
    else
        echo "Warning: Unknown genome '${GENOME}'. Converting VCF as-is."
        _prepare_basic_convert "${src_vcf}"
    fi
}

# ------------------------------------------------------------------
# GRCh38: already correct positions, coordinate-filter and sort
# ------------------------------------------------------------------
_prepare_grch38_subset() {
    local src_vcf="$1"
    local source_variant_count="$2"
    local sites_count="$3"

    _debug_log_command "bcftools view \"${src_vcf}\" --threads \"${THREADS}\" -Ou | bcftools view -T \"${RESOURCE_DIR}/GRCH38.sites.windows\" -Ou | bcftools sort -Oz -o \"${TMP_DIR}/peddy_input.vcf.gz\""
    bcftools view "${src_vcf}" --threads "${THREADS}" -Ou | \
    bcftools view -T "${RESOURCE_DIR}/GRCH38.sites.windows" -Ou | \
    bcftools sort -Oz -o "${TMP_DIR}/peddy_input.vcf.gz"
    bcftools index -t "${TMP_DIR}/peddy_input.vcf.gz"
    INPUT_VCF="${TMP_DIR}/peddy_input.vcf.gz"

    local n_variants
    n_variants=$(_count_variants "${INPUT_VCF}")
    echo "  Peddy input: ${n_variants} variants (GRCh38 coordinates)"
    local source_pct
    source_pct=$(_percent "${n_variants}" "${source_variant_count}")
    echo "  Final peddy subset retained ${source_pct}% of source VCF variants"
    _report_peddy_overlap "GRCh38 coordinate-window overlap" "source VCF variants" "${source_variant_count}" "${sites_count}" "${n_variants}"
}

# ------------------------------------------------------------------
# CHM13 / GRCh37: liftover to GRCh38, normalize, and coordinate-filter
# ------------------------------------------------------------------
_prepare_liftover_subset() {
    local src_vcf="$1"
    local source_variant_count="$2"
    local sites_count="$3"

    # Validate source FASTA
    if [[ -z "${SRC_FASTA}" || ! -f "${SRC_FASTA}" ]]; then
        echo "Error: --src-fasta is required for liftover (genome=${GENOME})" >&2
        echo "Falling back to basic VCF conversion (peddy may not find all sites)." >&2
        _prepare_basic_convert "${src_vcf}"
        return
    fi

    # Download chain file (source -> GRCh38)
    local chain_url=""
    case "${GENOME}" in
        CHM13)  chain_url="${CHAIN_URL_CHM13}" ;;
        GRCh37) chain_url="${CHAIN_URL_GRCh37}" ;;
    esac
    _cached_download "${chain_url}" \
        "${RESOURCE_DIR}/liftover_chain.gz" \
        "liftOver chain (${GENOME} -> GRCh38)"

    # Download / locate GRCh38 target reference FASTA
    local target_fasta=""
    target_fasta=$(_resolve_grch38_fasta)
    if [[ -z "${target_fasta}" || ! -f "${target_fasta}" ]]; then
        echo "Error: GRCh38 reference FASTA not available for liftover." >&2
        echo "Falling back to basic VCF conversion." >&2
        _prepare_basic_convert "${src_vcf}"
        return
    fi

    # For GRCh37 (bare chr names), the UCSC chain file expects chr-
    # prefixed names.  Rename the VCF and create a chr-prefixed copy
    # of the source reference so everything is consistent.
    local liftover_vcf="${src_vcf}"
    local liftover_src_fasta="${SRC_FASTA}"

    if [[ "${GENOME}" == "GRCh37" ]] && ! _has_chr_prefix "${src_vcf}"; then
        echo "  GRCh37 detected with bare chromosome names — adding chr prefix..."
        _create_add_chr_map "${TMP_DIR}/add_chr.txt"

        _debug_log_command "bcftools annotate --rename-chrs \"${TMP_DIR}/add_chr.txt\" \"${src_vcf}\" --threads \"${THREADS}\" -Ob -o \"${TMP_DIR}/chr_prefixed.bcf\""
        bcftools annotate --rename-chrs "${TMP_DIR}/add_chr.txt" \
            "${src_vcf}" --threads "${THREADS}" \
            -Ob -o "${TMP_DIR}/chr_prefixed.bcf"
        bcftools index "${TMP_DIR}/chr_prefixed.bcf"
        liftover_vcf="${TMP_DIR}/chr_prefixed.bcf"

        # Create chr-prefixed source reference (cached)
        local chr_ref
        chr_ref="${RESOURCE_DIR}/chr_$(basename "${SRC_FASTA}")"
        if [[ ! -f "${chr_ref}" ]]; then
            echo "  Creating chr-prefixed source reference (one-time)..."
            awk '{if (/^>/) {sub(/^>/, ">chr"); print} else print}' \
                "${SRC_FASTA}" > "${chr_ref}"
            samtools faidx "${chr_ref}"
        fi
        liftover_src_fasta="${chr_ref}"
    fi

    # Run the pipe: liftover | normalize to GRCh38 REF
    #              -> inline pre-site-filter variant count + coordinate-window subset.
    # We intentionally avoid writing the full lifted VCF and only persist the
    # peddy marker-window subset.
    echo "  Running bcftools +liftover pipeline..."
    local peddy_vcf="${TMP_DIR}/peddy_input.vcf.gz"
    local lifted_count_file="${TMP_DIR}/lifted_variants.count"
    rm -f "${lifted_count_file}"
    _debug_log_command "bcftools +liftover \"${liftover_vcf}\" -- -s \"${liftover_src_fasta}\" -f \"${target_fasta}\" -c \"${RESOURCE_DIR}/liftover_chain.gz\" --no-tags-update 2>\"${TMP_DIR}/liftover.log\" | bcftools norm -f \"${target_fasta}\" -c s -Ou 2>\"${TMP_DIR}/liftover.norm.log\" | tee >(bcftools view -H | wc -l > \"${lifted_count_file}\") | bcftools view -T \"${RESOURCE_DIR}/GRCH38.sites.windows\" -Ou | bcftools sort -Oz -o \"${peddy_vcf}\""
    bcftools +liftover "${liftover_vcf}" -- \
        -s "${liftover_src_fasta}" \
        -f "${target_fasta}" \
        -c "${RESOURCE_DIR}/liftover_chain.gz" \
        --no-tags-update 2>"${TMP_DIR}/liftover.log" | \
    # -c s: fix REF mismatches against GRCh38 target FASTA after coordinate liftover.
    bcftools norm -f "${target_fasta}" -c s -Ou 2>"${TMP_DIR}/liftover.norm.log" | \
    tee >(bcftools view -H | wc -l > "${lifted_count_file}") | \
    bcftools view -T "${RESOURCE_DIR}/GRCH38.sites.windows" -Ou | \
    bcftools sort -Oz -o "${peddy_vcf}"

    local wait_attempts=0
    while (( wait_attempts < PEDDY_COUNT_WAIT_ATTEMPTS )); do
        if [[ -f "${lifted_count_file}" ]]; then
            break
        fi
        ((wait_attempts += 1))
        sleep 0.1
    done

    bcftools index -t "${peddy_vcf}"
    INPUT_VCF="${peddy_vcf}"

    local lifted_variants
    lifted_variants="0"
    if [[ -f "${lifted_count_file}" ]]; then
        lifted_variants=$(tr -d '[:space:]' < "${lifted_count_file}")
    else
        echo "Warning: Streamed liftover variant count file was not produced (liftover stream failed early or async count write did not complete)."
    fi
    if [[ ! "${lifted_variants}" =~ ^[0-9]+$ ]]; then
        echo "Warning: Failed to read streamed liftover variant count; defaulting to 0."
        lifted_variants="0"
    fi
    local lifted_pct
    lifted_pct=$(_percent "${lifted_variants}" "${source_variant_count}")
    echo "  Liftover output (pre-site-filter): ${lifted_variants} variants (${lifted_pct}% of source variants)"

    if (( lifted_variants == 0 )); then
        echo "Warning: Liftover produced zero variants. Verify source reference FASTA and chain file compatibility."
    elif awk -v n="${lifted_variants}" -v d="${source_variant_count}" -v thr="${PEDDY_MIN_LIFTOVER_RETAIN_FRACTION}" \
        'BEGIN{ exit !(d > 0 && (n/d) < thr) }'; then
        echo "Warning: Liftover retained only ${lifted_pct}% of source variants. This suggests a possible build or contig naming mismatch."
    fi

    # Report liftover stats
    if [[ -f "${TMP_DIR}/liftover.log" ]]; then
        echo "  Liftover log:"
        sed 's/^/    /' "${TMP_DIR}/liftover.log" | head -20
    fi
    if [[ -s "${TMP_DIR}/liftover.norm.log" ]]; then
        echo "  Liftover normalization log:"
        sed 's/^/    /' "${TMP_DIR}/liftover.norm.log" | head -10
    fi

    local n_variants
    n_variants=$(_count_variants "${INPUT_VCF}")
    echo "  Peddy input: ${n_variants} variants (lifted to GRCh38 coordinates)"
    _report_peddy_overlap "Liftover coordinate-window overlap" "lifted variants" "${lifted_variants}" "${sites_count}" "${n_variants}"
}

# ------------------------------------------------------------------
# Resolve path to GRCh38 target reference FASTA (download if needed)
# ------------------------------------------------------------------
_resolve_grch38_fasta() {
    # Check common locations
    local candidates=(
        "${REF_DIR}/${GRCH38_FASTA_NAME}"
        "${RESOURCE_DIR}/${GRCH38_FASTA_NAME}"
    )
    for f in "${candidates[@]}"; do
        if [[ -f "${f}" && -f "${f}.fai" ]]; then
            echo "${f}"
            return
        fi
    done

    # Download into ref-dir or resources
    local dest_dir="${REF_DIR:-${RESOURCE_DIR}}"
    local dest="${dest_dir}/${GRCH38_FASTA_NAME}"

    if [[ ! -f "${dest}" ]]; then
        echo "  Downloading GRCh38 reference FASTA for liftover target..."
        if wget -q "${GRCH38_FASTA_URL}" -O "${dest}.gz" 2>/dev/null || \
           curl -sfL "${GRCH38_FASTA_URL}" -o "${dest}.gz" 2>/dev/null; then
            echo "  Decompressing GRCh38 reference..."
            gunzip -f "${dest}.gz"
        else
            echo "Error: Failed to download GRCh38 reference." >&2
            rm -f "${dest}.gz"
            return 1
        fi
    fi

    # Index if needed
    if [[ ! -f "${dest}.fai" ]]; then
        echo "  Indexing GRCh38 reference..."
        samtools faidx "${dest}"
    fi

    echo "${dest}"
}

# ------------------------------------------------------------------
# Fallback: basic BCF -> VCF.gz (no liftover, no site filtering)
# ------------------------------------------------------------------
_prepare_basic_convert() {
    local src_vcf="$1"
    INPUT_VCF="${TMP_DIR}/peddy_input.vcf.gz"

    if [[ "${src_vcf}" == *.bcf ]]; then
        echo "  Converting BCF to VCF.gz for peddy..."
        _debug_log_command "bcftools view \"${src_vcf}\" -Oz -o \"${INPUT_VCF}\" --threads \"${THREADS}\""
        bcftools view "${src_vcf}" -Oz -o "${INPUT_VCF}" --threads "${THREADS}"
    elif [[ "${src_vcf}" == *.vcf && ! "${src_vcf}" == *.vcf.gz ]]; then
        echo "  Compressing VCF to VCF.gz for peddy..."
        _debug_log_command "bcftools view \"${src_vcf}\" -Oz -o \"${INPUT_VCF}\" --threads \"${THREADS}\""
        bcftools view "${src_vcf}" -Oz -o "${INPUT_VCF}" --threads "${THREADS}"
    elif [[ "${src_vcf}" == *.vcf.gz ]]; then
        INPUT_VCF="${src_vcf}"
        _ensure_index "${src_vcf}"
        return
    fi
    bcftools index -t "${INPUT_VCF}"
}

# ======================================================================
#  Main execution
# ======================================================================

INPUT_VCF=""
_prepare_peddy_input "${VCF}"

if [[ -z "${INPUT_VCF}" || ! -f "${INPUT_VCF}" ]]; then
    echo "Error: Failed to prepare peddy input VCF" >&2
    exit 1
fi

# --- Generate PED file if not provided ---
USED_PED="${PED_FILE}"
GENERATED_PED="${OUTPUT_DIR}/generated_input.ped"

if [[ -z "${PED_FILE}" || ! -f "${PED_FILE}" ]]; then
    echo "No pedigree file provided. Generating minimal PED from VCF samples..."

    bcftools query -l "${INPUT_VCF}" > "${TMP_DIR}/sample_list.txt"

    declare -A sex_map
    if [[ -n "${SAMPLE_QC}" && -f "${SAMPLE_QC}" ]]; then
        while IFS=$'\t' read -r sid _ _ _ _ _ _ gender _; do
            case "${gender}" in
                M|male|1)   sex_map["${sid}"]="1" ;;
                F|female|2) sex_map["${sid}"]="2" ;;
                *)          sex_map["${sid}"]="0" ;;
            esac
        done < <(tail -n +2 "${SAMPLE_QC}")
    fi

    : > "${GENERATED_PED}"
    while read -r sample; do
        sex="${sex_map[${sample}]:-0}"
        echo -e "${sample}\t${sample}\t0\t0\t${sex}\t-9" >> "${GENERATED_PED}"
    done < "${TMP_DIR}/sample_list.txt"

    USED_PED="${GENERATED_PED}"
    echo "  Generated PED: ${GENERATED_PED} ($(wc -l < "${GENERATED_PED}") samples)"
else
    echo "Using provided pedigree file: ${PED_FILE}"
    cp "${PED_FILE}" "${OUTPUT_DIR}/input.ped"
fi

# --- Run peddy ---
echo ""
echo "Running peddy..."
echo "  VCF:     ${INPUT_VCF}"
echo "  PED:     ${USED_PED}"
echo "  Genome:  ${GENOME}"
echo "  Prefix:  ${PREFIX}"
echo "  Threads: ${THREADS}"
echo ""

PEDDY_START=${SECONDS}

python3 -m peddy \
    -p "${THREADS}" \
    --plot \
    --prefix "${PREFIX}" \
    "${INPUT_VCF}" \
    "${USED_PED}" 2>&1 || {
        echo "Warning: peddy exited with non-zero status. Some outputs may be incomplete." >&2
    }

PEDDY_ELAPSED=$(( SECONDS - PEDDY_START ))
echo ""
echo "Peddy completed in $(( PEDDY_ELAPSED / 60 ))m $(( PEDDY_ELAPSED % 60 ))s"

# --- Create final pedigree with discovered relationships ---
FINAL_PED="${OUTPUT_DIR}/peddy_final.ped"

if [[ -f "${PREFIX}.peddy.ped" ]]; then
    echo "Creating final pedigree from peddy results..."
    cp "${PREFIX}.peddy.ped" "${FINAL_PED}"
    echo "  Final pedigree: ${FINAL_PED}"
fi

# --- Report output summary ---
echo ""
echo "Peddy output files:"
for f in "${EXPECTED_OUTPUTS[@]}"; do
    if [[ -f "$f" ]]; then
        row_count=$(wc -l < "$f" | tr -d ' ')
        echo "  $(basename "$f"): ${row_count} lines"
    else
        echo "  $(basename "$f"): MISSING"
    fi
done

if [[ -f "${PREFIX}.html" ]]; then
    echo "  peddy.html: interactive report"
fi
if [[ -f "${PREFIX}.background_pca.json" ]]; then
    echo "  background_pca.json: ancestry PCA data"
fi

# Cleanup temporary files
rm -rf "${TMP_DIR}"

echo ""
echo "Done."
