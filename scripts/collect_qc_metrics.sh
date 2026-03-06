#!/usr/bin/env bash
#
# collect_qc_metrics.sh
#
# Helper functions to extract per-sample QC metrics from genotyping output.
# Computes call rate, LRR SD, LRR mean, LRR median, BAF SD, and
# heterozygosity rate per sample in a SINGLE BCF pass.
#
# Performance: Previously required two separate BCF decompression passes
# (one for call rate/LRR SD via awk, one for LRR median via Python).
# Now uses a single pass with compute_sample_qc.py, halving I/O.
#
# Sourced by stage1 and stage2 scripts; not intended to be run directly.
#

collect_qc_metrics() {
    local vcf_file="$1"
    local metadata_tsv="$2"
    local output_file="$3"
    local threads="${4:-1}"

    local script_dir
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    local qc_start qc_elapsed
    qc_start=${SECONDS}
    echo "Extracting per-sample QC metrics..."

    local gender_file qc_metrics_file samples_file
    gender_file=$(mktemp)
    qc_metrics_file=$(mktemp)
    samples_file=$(mktemp)
    trap 'rm -f "${gender_file}" "${qc_metrics_file}" "${samples_file}"' RETURN

    # Get sample names from VCF
    bcftools query -l "${vcf_file}" > "${samples_file}"
    local n_samples_vcf
    n_samples_vcf=$(wc -l < "${samples_file}" | tr -d ' ')

    # -----------------------------------------------------------------
    # Diagnostic: Variant counts (use index when possible)
    # -----------------------------------------------------------------
    echo "  [diag] === VCF variant diagnostics ==="
    local n_total_vars n_intensity_only n_genotypeable
    # Use bcftools index -n for instant count from CSI index
    n_total_vars=$(bcftools index -n "${vcf_file}" 2>/dev/null || \
        bcftools view -H --threads "${threads}" "${vcf_file}" 2>/dev/null | wc -l | tr -d ' ')
    # Count intensity-only probes via bcftools view -c (avoids full BCF scan
    # when the index supports region queries).  Falls back to streaming count.
    n_intensity_only=$(bcftools view -i 'INFO/INTENSITY_ONLY=1' -H --threads "${threads}" "${vcf_file}" 2>/dev/null | wc -l | tr -d ' ')
    n_genotypeable=$(( n_total_vars - n_intensity_only ))
    echo "  [diag] Total variants in VCF:          ${n_total_vars}"
    echo "  [diag] Intensity-only probes:           ${n_intensity_only}"
    echo "  [diag] Genotypeable variants:            ${n_genotypeable}"
    echo "  [diag] Samples in VCF:                   ${n_samples_vcf}"
    echo "  [diag] ==================================="

    # -----------------------------------------------------------------
    # Compute ALL per-sample metrics in a SINGLE BCF pass.
    #
    # Extracts GT:LRR:BAF triples per sample per variant, piped to
    # compute_sample_qc.py which computes call_rate, lrr_sd, lrr_mean,
    # lrr_median, baf_sd, and het_rate simultaneously.
    #
    # This replaces the previous two-pass approach (awk + Python median),
    # halving BCF decompression and I/O cost.
    #
    # Restricted to autosomes only to avoid inflation from sex-chromosome
    # hemizygosity in males.
    # -----------------------------------------------------------------
    echo "  [diag] Starting autosomal GT/LRR/BAF stream extraction..."
    local metric_start metric_elapsed
    metric_start=${SECONDS}
    local autosome_targets
    autosome_targets="$(printf 'chr%s,' {1..22}; printf '%s,' {1..22})"
    autosome_targets="${autosome_targets%,}"
    echo "  Computing QC metrics per sample (${n_samples_vcf} samples, autosomes only, single pass)..."
    bcftools view -e 'INFO/INTENSITY_ONLY=1' -t "${autosome_targets}" --threads "${threads}" "${vcf_file}" 2>/dev/null | \
    bcftools query -f '[\t%GT:%LRR:%BAF]\n' 2>/dev/null | \
        python3 "${script_dir}/compute_sample_qc.py" \
            --num-samples "${n_samples_vcf}" \
            --output "${qc_metrics_file}.raw" 2>&1 | sed 's/^/  /'
    metric_elapsed=$(( SECONDS - metric_start ))
    echo "  [diag] Sample metric extraction completed in ${metric_elapsed}s"

    # Map column indices back to sample names
    if [[ -s "${qc_metrics_file}.raw" ]]; then
        # Skip header from compute_sample_qc.py, map idx to sample names
        awk -F'\t' 'NR==FNR { names[NR-1] = $0; next }
             FNR > 1 { print names[$1+0] "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 }' \
            "${samples_file}" "${qc_metrics_file}.raw" | \
            sort -k1,1 > "${qc_metrics_file}"
        rm -f "${qc_metrics_file}.raw"
    else
        echo "  [diag] WARNING: QC metrics computation returned empty output"
        touch "${qc_metrics_file}"
    fi

    if [[ -s "${qc_metrics_file}" ]]; then
        n_qm=$(wc -l < "${qc_metrics_file}" | tr -d ' ')
        echo "  [diag] QC metrics file: ${n_qm} samples"
        echo "  [diag] QC metrics first 3 lines (sample, call_rate, lrr_sd, lrr_mean, lrr_median, baf_sd, het_rate):"
        head -3 "${qc_metrics_file}" | sed 's/^/    [diag]   /'
    else
        echo "  [diag] WARNING: QC metrics file is empty"
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

    # Merge QC metrics with gender into a single output file
    local merge_start merge_elapsed
    merge_start=${SECONDS}
    echo "  Merging metrics..."
    {
        echo -e "sample_id\tcall_rate\tlrr_sd\tlrr_mean\tlrr_median\tbaf_sd\thet_rate\tcomputed_gender"
        join -t$'\t' -a1 -e 'NA' -o '0,1.2,1.3,1.4,1.5,1.6,1.7,2.2' "${qc_metrics_file}" "${gender_file}"
    } > "${output_file}" 2>/dev/null || {
        # Fallback: simpler merge approach
        echo -e "sample_id\tcall_rate\tlrr_sd\tlrr_mean\tlrr_median\tbaf_sd\thet_rate\tcomputed_gender"
        paste "${qc_metrics_file}" "${gender_file}" | \
            awk -F'\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $9}'
    } > "${output_file}"

    merge_elapsed=$(( SECONDS - merge_start ))
    n_samples=$(( $(wc -l < "${output_file}") - 1 ))
    echo "  QC metrics computed for ${n_samples} samples: ${output_file}"
    echo "  [diag] Merge + finalization completed in ${merge_elapsed}s"

    # Diagnostic: show first few rows and value ranges
    echo "  [diag] QC output first 5 lines:"
    head -5 "${output_file}" | sed 's/^/    [diag]   /'
    if [[ "${n_samples}" -gt 0 ]]; then
        awk -F'\t' '
            NR==1 { for(i=1;i<=NF;i++) if($i=="call_rate") c=i; next }
            c && $c != "NA" {
                n++; sum+=$c
                if (n==1 || $c+0 < min) min=$c+0
                if (n==1 || $c+0 > max) max=$c+0
            } END {
                if (n>0) printf "  [diag] Call rate range: %.4f - %.4f (mean %.4f, n=%d)\n", min, max, sum/n, n
                else print "  [diag] WARNING: No valid call rate values found"
            }' "${output_file}"
        awk -F'\t' '
            NR==1 { for(i=1;i<=NF;i++) if($i=="lrr_sd") c=i; next }
            c && $c != "NA" {
                n++; sum+=$c
                if (n==1 || $c+0 < min) min=$c+0
                if (n==1 || $c+0 > max) max=$c+0
            } END {
                if (n>0) printf "  [diag] LRR SD range: %.4f - %.4f (mean %.4f, n=%d)\n", min, max, sum/n, n
                else print "  [diag] WARNING: No valid LRR SD values found"
            }' "${output_file}"
        awk -F'\t' '
            NR==1 { for(i=1;i<=NF;i++) if($i=="baf_sd") c=i; next }
            c && $c != "NA" {
                n++; sum+=$c
                if (n==1 || $c+0 < min) min=$c+0
                if (n==1 || $c+0 > max) max=$c+0
            } END {
                if (n>0) printf "  [diag] BAF SD range: %.4f - %.4f (mean %.4f, n=%d)\n", min, max, sum/n, n
                else print "  [diag] WARNING: No valid BAF SD values found"
            }' "${output_file}"
        awk -F'\t' '
            NR==1 { for(i=1;i<=NF;i++) if($i=="het_rate") c=i; next }
            c && $c != "NA" {
                n++; sum+=$c
                if (n==1 || $c+0 < min) min=$c+0
                if (n==1 || $c+0 > max) max=$c+0
            } END {
                if (n>0) printf "  [diag] Het rate range: %.4f - %.4f (mean %.4f, n=%d)\n", min, max, sum/n, n
                else print "  [diag] WARNING: No valid het rate values found"
            }' "${output_file}"
    fi

    qc_elapsed=$(( SECONDS - qc_start ))
    echo "  [diag] Total sample QC metric stage time: ${qc_elapsed}s"

    rm -f "${samples_file}"
}
