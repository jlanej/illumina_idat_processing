#!/usr/bin/env bash
#
# collect_qc_metrics.sh
#
# Helper functions to extract per-sample QC metrics from genotyping output.
# Computes call rate and LRR standard deviation per sample.
#
# Sourced by stage1 and stage2 scripts; not intended to be run directly.
#

collect_qc_metrics() {
    local vcf_file="$1"
    local metadata_tsv="$2"
    local output_file="$3"

    echo "Extracting per-sample QC metrics..."

    # Extract computed gender and call rate from GTC metadata if available
    local gender_file call_rate_file
    gender_file=$(mktemp)
    call_rate_file=$(mktemp)
    lrr_sd_file=$(mktemp)
    trap 'rm -f "${gender_file}" "${call_rate_file}" "${lrr_sd_file}"' RETURN

    # Get sample names from VCF
    local samples_file
    samples_file=$(mktemp)
    bcftools query -l "${vcf_file}" > "${samples_file}"

    # -----------------------------------------------------------------
    # Diagnostic: Variant counts at each filtering stage
    # -----------------------------------------------------------------
    echo "  [diag] === VCF variant diagnostics ==="
    local n_total_vars n_intensity_only n_genotypeable n_missing_ref
    n_total_vars=$(bcftools view -H "${vcf_file}" 2>/dev/null | wc -l | tr -d ' ')
    n_intensity_only=$(bcftools view -i 'INFO/INTENSITY_ONLY=1' -H "${vcf_file}" 2>/dev/null | wc -l | tr -d ' ')
    n_genotypeable=$(( n_total_vars - n_intensity_only ))
    echo "  [diag] Total variants in VCF:          ${n_total_vars}"
    echo "  [diag] Intensity-only probes:           ${n_intensity_only}"
    echo "  [diag] Genotypeable variants:            ${n_genotypeable}"

    # Count variants with all-missing genotypes (indicator of norm -c x damage)
    n_all_missing=$(bcftools view -e 'INFO/INTENSITY_ONLY=1' "${vcf_file}" 2>/dev/null | \
        bcftools view -e 'COUNT(GT!="mis")>0' -H 2>/dev/null | wc -l | tr -d ' ') || n_all_missing="N/A"
    echo "  [diag] Variants with ALL samples missing: ${n_all_missing}"

    # Check for first sample: how many genotypes are called vs missing
    local first_sample
    first_sample=$(head -1 "${samples_file}")
    if [[ -n "${first_sample}" ]]; then
        local n_called_first n_missing_first
        n_called_first=$(bcftools view -e 'INFO/INTENSITY_ONLY=1' -s "${first_sample}" "${vcf_file}" 2>/dev/null | \
            bcftools query -f '[%GT]\n' 2>/dev/null | grep -cv '^\./\.\|^\.\|^\./\.' || true)
        n_missing_first=$(bcftools view -e 'INFO/INTENSITY_ONLY=1' -s "${first_sample}" "${vcf_file}" 2>/dev/null | \
            bcftools query -f '[%GT]\n' 2>/dev/null | grep -c '^\./\.\|^\.\|^\./\.' || true)
        echo "  [diag] Sample '${first_sample}' genotype breakdown:"
        echo "  [diag]   Called:  ${n_called_first}"
        echo "  [diag]   Missing: ${n_missing_first}"
        if [[ "${n_genotypeable}" -gt 0 ]]; then
            local pct_called
            pct_called=$(awk -v c="${n_called_first}" -v t="${n_genotypeable}" 'BEGIN { printf "%.4f", (t>0) ? c/t : 0 }')
            echo "  [diag]   Call rate: ${pct_called}"
        fi
    fi
    echo "  [diag] ==================================="

    # Compute call rate per sample: fraction of non-missing genotypes
    # Exclude intensity-only probes (no genotype, only BAF/LRR) from the
    # denominator so call rate reflects genotypeable loci only.
    echo "  Computing call rate per sample..."
    bcftools view -e 'INFO/INTENSITY_ONLY=1' "${vcf_file}" | \
    bcftools query -f '[%SAMPLE\t%GT\n]' | \
        awk -F'\t' '{
            total[$1]++
            if ($2 != "./." && $2 != "." && $2 != ".|.") called[$1]++
        } END {
            for (s in total) {
                cr = (total[s] > 0) ? called[s] / total[s] : 0
                print s "\t" cr
            }
        }' | sort -k1,1 > "${call_rate_file}" 2>/dev/null || true

    # Diagnostic: log call rate file stats
    if [[ -s "${call_rate_file}" ]]; then
        n_cr=$(wc -l < "${call_rate_file}")
        echo "  [diag] Call rate file: ${n_cr} samples"
        echo "  [diag] Call rate first 3 lines:"
        head -3 "${call_rate_file}" | sed 's/^/    [diag]   /'
    else
        echo "  [diag] WARNING: Call rate file is empty"
    fi

    # Compute LRR standard deviation per sample
    # Filter out non-numeric values (nan, -nan, inf, -inf) that gtc2vcf
    # outputs for probes with zero intensities or failed normalization.
    echo "  Computing LRR standard deviation per sample..."
    bcftools query -f '[%SAMPLE\t%LRR\n]' "${vcf_file}" 2>/dev/null | \
        awk -F'\t' '$2 != "." && $2 != "" && $2 ~ /^-?[0-9]/ {
            n[$1]++
            sum[$1] += $2
            sum2[$1] += $2 * $2
        } END {
            for (s in n) {
                if (n[s] > 1) {
                    mean = sum[s] / n[s]
                    var = (sum2[s] / n[s]) - (mean * mean)
                    if (var < 0) var = 0
                    sd = sqrt(var)
                    print s "\t" sd
                } else {
                    print s "\tNA"
                }
            }
        }' | sort -k1,1 > "${lrr_sd_file}" 2>/dev/null || true

    # Diagnostic: log LRR SD file stats
    if [[ -s "${lrr_sd_file}" ]]; then
        n_sd=$(wc -l < "${lrr_sd_file}")
        echo "  [diag] LRR SD file: ${n_sd} samples"
        echo "  [diag] LRR SD first 3 lines:"
        head -3 "${lrr_sd_file}" | sed 's/^/    [diag]   /'
    else
        echo "  [diag] WARNING: LRR SD file is empty (LRR FORMAT field may not exist in VCF)"
    fi

    # Extract computed gender from metadata if available
    # The gtc2vcf --extra TSV uses "gtc" (filename) as the sample identifier
    # column, not "sample_id". Strip the .gtc extension to match VCF sample names.
    if [[ -f "${metadata_tsv}" ]] && head -1 "${metadata_tsv}" | grep -q "computed_gender"; then
        awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) {if($i=="sample_id"||$i=="gtc") si=i; if($i=="computed_gender") gi=i}}
                     NR>1 && si && gi {
                         id = $si
                         sub(/\.gtc$/, "", id)
                         sub(/.*\//, "", id)
                         print id "\t" $gi
                     }' "${metadata_tsv}" | sort -k1,1 > "${gender_file}"
    else
        # Placeholder
        while read -r s; do echo -e "${s}\tNA"; done < "${samples_file}" > "${gender_file}"
    fi

    # Merge all metrics into a single output file
    echo "  Merging metrics..."
    {
        echo -e "sample_id\tcall_rate\tlrr_sd\tcomputed_gender"
        join -t$'\t' -a1 -e 'NA' -o '0,1.2,2.2' "${call_rate_file}" "${lrr_sd_file}" | \
            join -t$'\t' -a1 -e 'NA' -o '0,1.2,1.3,2.2' - "${gender_file}"
    } > "${output_file}" 2>/dev/null || {
        # Fallback: simpler merge approach
        echo -e "sample_id\tcall_rate\tlrr_sd\tcomputed_gender"
        paste "${call_rate_file}" "${lrr_sd_file}" "${gender_file}" | \
            awk -F'\t' '{print $1 "\t" $2 "\t" $4 "\t" $6}'
    } > "${output_file}"

    n_samples=$(( $(wc -l < "${output_file}") - 1 ))
    echo "  QC metrics computed for ${n_samples} samples: ${output_file}"

    # Diagnostic: show first few rows and value ranges
    echo "  [diag] QC output first 5 lines:"
    head -5 "${output_file}" | sed 's/^/    [diag]   /'
    if [[ "${n_samples}" -gt 0 ]]; then
        awk -F'\t' 'NR>1 && $2 != "NA" {
            n++; sum+=$2
            if (n==1 || $2+0 < min) min=$2+0
            if (n==1 || $2+0 > max) max=$2+0
        } END {
            if (n>0) printf "  [diag] Call rate range: %.4f - %.4f (mean %.4f, n=%d)\n", min, max, sum/n, n
            else print "  [diag] WARNING: No valid call rate values found"
        }' "${output_file}"
        awk -F'\t' 'NR>1 && $3 != "NA" {
            n++; sum+=$3
            if (n==1 || $3+0 < min) min=$3+0
            if (n==1 || $3+0 > max) max=$3+0
        } END {
            if (n>0) printf "  [diag] LRR SD range: %.4f - %.4f (mean %.4f, n=%d)\n", min, max, sum/n, n
            else print "  [diag] WARNING: No valid LRR SD values found"
        }' "${output_file}"
    fi

    rm -f "${samples_file}"
}
