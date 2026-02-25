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

    # Compute call rate per sample: fraction of non-missing genotypes
    echo "  Computing call rate per sample..."
    bcftools +fill-tags "${vcf_file}" -Ou -- -t F_MISSING | \
        bcftools query -f '[%SAMPLE\t%F_MISSING\n]' | \
        awk -F'\t' '{
            missing[$1] += $2
            count[$1]++
        } END {
            for (s in missing) {
                avg_missing = missing[s] / count[s]
                call_rate = 1 - avg_missing
                print s "\t" call_rate
            }
        }' | sort -k1,1 > "${call_rate_file}" 2>/dev/null || true

    # If fill-tags approach fails, compute directly
    if [[ ! -s "${call_rate_file}" ]]; then
        echo "  Computing call rate from genotypes directly..."
        bcftools query -f '%CHROM\t%POS[\t%GT]\n' "${vcf_file}" | \
            awk -F'\t' 'BEGIN {OFS="\t"}
            NR==1 {n_samples = NF - 2}
            {
                for (i=3; i<=NF; i++) {
                    total[i]++
                    if ($i != "./." && $i != "." && $i != "./." ) called[i]++
                }
            }
            END {
                for (i=3; i<=n_samples+2; i++) {
                    print i-2, (total[i] > 0 ? called[i]/total[i] : 0)
                }
            }' > "${call_rate_file}"
    fi

    # Compute LRR standard deviation per sample
    echo "  Computing LRR standard deviation per sample..."
    bcftools query -f '[%SAMPLE\t%LRR\n]' "${vcf_file}" 2>/dev/null | \
        awk -F'\t' '$2 != "." && $2 != "" {
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

    # Extract computed gender from metadata if available
    if [[ -f "${metadata_tsv}" ]] && head -1 "${metadata_tsv}" | grep -q "computed_gender"; then
        awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="sample_id") si=i; if($i=="computed_gender") gi=i}
                     NR>1 {print $si "\t" $gi}' "${metadata_tsv}" | sort -k1,1 > "${gender_file}"
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

    rm -f "${samples_file}"
}
