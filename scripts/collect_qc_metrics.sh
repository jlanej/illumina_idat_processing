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
    local threads="${4:-1}"

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
        bcftools view -H "${vcf_file}" 2>/dev/null | wc -l | tr -d ' ')
    n_intensity_only=$(bcftools view -i 'INFO/INTENSITY_ONLY=1' -H "${vcf_file}" 2>/dev/null | wc -l | tr -d ' ')
    n_genotypeable=$(( n_total_vars - n_intensity_only ))
    echo "  [diag] Total variants in VCF:          ${n_total_vars}"
    echo "  [diag] Intensity-only probes:           ${n_intensity_only}"
    echo "  [diag] Genotypeable variants:            ${n_genotypeable}"
    echo "  [diag] Samples in VCF:                   ${n_samples_vcf}"
    echo "  [diag] ==================================="

    # -----------------------------------------------------------------
    # Compute call rate AND LRR SD per sample in a single BCF pass.
    #
    # Instead of making two separate passes over the (potentially huge)
    # BCF file — one for GT and one for LRR — we extract both fields at
    # once with bcftools query -f '[\t%GT:%LRR]\n'.  The awk script
    # splits each sample's GT:LRR pair and accumulates both call-rate
    # counts and running LRR sum/sum-of-squares simultaneously.
    #
    # This halves the BCF decompression and I/O cost of QC computation.
    #
    # Restricted to autosomes only to avoid inflation from sex-chromosome
    # hemizygosity in males.  Filters non-numeric LRR values (nan, inf)
    # that gtc2vcf can produce for zero-intensity probes.
    # -----------------------------------------------------------------
    echo "  Computing call rate and LRR SD per sample (${n_samples_vcf} samples, autosomes only)..."
    bcftools view -e 'INFO/INTENSITY_ONLY=1' -t ^chrX,chrY,chrM,X,Y,MT --threads "${threads}" "${vcf_file}" 2>/dev/null | \
    bcftools query -f '[\t%GT:%LRR]\n' 2>/dev/null | \
        awk -F'\t' '
        {
            for (i = 2; i <= NF; i++) {
                split($i, a, ":")
                gt = a[1]
                lrr = a[2]
                # Call rate: count called vs total genotypes
                total[i]++
                if (gt != "./." && gt != "." && gt != ".|.") called[i]++
                # LRR SD: accumulate sum and sum-of-squares
                if (lrr != "." && lrr != "" && lrr ~ /^-?[0-9]/) {
                    n_lrr[i]++
                    sum_lrr[i] += lrr
                    sum2_lrr[i] += lrr * lrr
                }
            }
        }
        END {
            for (i = 2; i <= NF; i++) {
                idx = i - 2
                cr = (total[i] > 0) ? called[i] / total[i] : 0
                if (n_lrr[i] > 1) {
                    mean = sum_lrr[i] / n_lrr[i]
                    var = (sum2_lrr[i] / n_lrr[i]) - (mean * mean)
                    if (var < 0) var = 0
                    sd = sqrt(var)
                } else {
                    sd = "NA"
                }
                print idx "\t" cr "\t" sd
            }
        }' > "${qc_metrics_file}.idx" 2>/dev/null || true

    # Map column indices back to sample names (combined call_rate + lrr_sd)
    awk 'NR==FNR { names[NR-1] = $0; next }
         { print names[$1+0] "\t" $2 "\t" $3 }' \
        "${samples_file}" "${qc_metrics_file}.idx" | \
        sort -k1,1 > "${qc_metrics_file}"
    rm -f "${qc_metrics_file}.idx"

    if [[ -s "${qc_metrics_file}" ]]; then
        n_qm=$(wc -l < "${qc_metrics_file}" | tr -d ' ')
        echo "  [diag] QC metrics file: ${n_qm} samples"
        echo "  [diag] QC metrics first 3 lines (sample, call_rate, lrr_sd):"
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
    echo "  Merging metrics..."
    {
        echo -e "sample_id\tcall_rate\tlrr_sd\tcomputed_gender"
        join -t$'\t' -a1 -e 'NA' -o '0,1.2,1.3,2.2' "${qc_metrics_file}" "${gender_file}"
    } > "${output_file}" 2>/dev/null || {
        # Fallback: simpler merge approach
        echo -e "sample_id\tcall_rate\tlrr_sd\tcomputed_gender"
        paste "${qc_metrics_file}" "${gender_file}" | \
            awk -F'\t' '{print $1 "\t" $2 "\t" $3 "\t" $5}'
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
