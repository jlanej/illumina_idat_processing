#!/usr/bin/env bash
#
# test_resource_links.sh
#
# Validate that array manifest and reference resource URLs are reachable.
# Checks HTTP status codes for all URLs in config/manifest_urls.tsv and
# key external resources used by the pipeline.
#
# NOTE: Illumina has retired their webdata.illumina.com direct download URLs.
# Manifest URLs now point to support.illumina.com product file pages.
# These pages require free registration to download files, but the pages
# themselves should be reachable (HTTP 200).
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
MANIFEST_FILE="${REPO_DIR}/config/manifest_urls.tsv"

PASS=0
FAIL=0
SKIP=0

check_url() {
    local url="$1"
    local description="$2"

    # Skip FTP URLs (cannot check with curl HEAD request reliably)
    if [[ "${url}" == ftp://* ]]; then
        echo "  SKIP: ${description} (FTP URL)"
        (( SKIP++ )) || true
        return 0
    fi

    local http_code
    http_code=$(curl -o /dev/null -s -w "%{http_code}" --head --max-time 30 \
        -L "${url}" 2>/dev/null) || http_code="000"

    if [[ "${http_code}" =~ ^(200|301|302|303|307|308)$ ]]; then
        echo "  PASS: ${description} (HTTP ${http_code})"
        (( PASS++ )) || true
    else
        echo "  FAIL: ${description} (HTTP ${http_code})"
        echo "        URL: ${url}"
        (( FAIL++ )) || true
    fi
}

echo "============================================"
echo "  Testing Resource URLs"
echo "============================================"
echo ""

# Test manifest URLs from config file
# Since Illumina has moved to support.illumina.com, BPM/EGT/CSV URLs now
# point to the same product file download page. We deduplicate and check
# each unique URL once per array.
echo "--- Manifest URLs (config/manifest_urls.tsv) ---"
if [[ -f "${MANIFEST_FILE}" ]]; then
    while IFS=$'\t' read -r array_name bpm_url egt_url csv_url; do
        [[ "${array_name}" =~ ^# ]] && continue
        [[ -z "${array_name}" ]] && continue

        echo ""
        echo "  Array: ${array_name}"

        # Collect unique URLs for this array
        declare -A seen_urls=()
        for url in "${bpm_url}" "${egt_url}" "${csv_url}"; do
            if [[ "${url}" == ftp://* ]]; then
                echo "  SKIP: ${array_name} (FTP URL)"
                (( SKIP++ )) || true
            elif [[ -z "${seen_urls[${url}]+_}" ]]; then
                seen_urls["${url}"]=1
                check_url "${url}" "${array_name} product files"
            fi
        done
        unset seen_urls
    done < "${MANIFEST_FILE}"
else
    echo "  FAIL: Manifest URL file not found: ${MANIFEST_FILE}"
    (( FAIL++ )) || true
fi

echo ""
echo "--- 1000 Genomes Resources ---"

# 1000G IDAT archive
check_url "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/broad_intensities/" \
    "1000G broad_intensities directory"

# 1000G Omni2.5 manifest files (HumanOmni2.5-4v1)
check_url "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/broad_intensities/HumanOmni2.5-4v1-Multi_B.bpm" \
    "1000G Omni2.5 BPM manifest"

check_url "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/broad_intensities/HumanOmni2.5-4v1-Multi_B.egt" \
    "1000G Omni2.5 EGT cluster file"

check_url "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/broad_intensities/HumanOmni2.5-4v1_B.csv" \
    "1000G Omni2.5 CSV manifest"

echo ""
echo "--- Reference Genome Resources ---"

check_url "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz" \
    "CHM13v2.0 reference genome (T2T, default)"

check_url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" \
    "GRCh38 reference genome"

check_url "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz" \
    "GRCh37 reference genome"

echo ""
echo "--- Dependency Sources ---"

check_url "https://github.com/samtools/bcftools/releases/download/1.23/bcftools-1.23.tar.bz2" \
    "bcftools 1.23 release"

check_url "https://raw.githubusercontent.com/freeseek/gtc2vcf/master/gtc2vcf.c" \
    "gtc2vcf plugin source"

check_url "https://raw.githubusercontent.com/freeseek/mocha/master/mocha.c" \
    "mocha plugin source"

echo ""
echo "============================================"
echo "  Results: ${PASS} passed, ${FAIL} failed, ${SKIP} skipped"
echo "============================================"

if [[ "${FAIL}" -gt 0 ]]; then
    echo ""
    echo "WARNING: ${FAIL} resource URL(s) are unreachable."
    echo "This may indicate broken links that need updating."
    exit 1
fi
