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

    local gender_file call_rate_file lrr_sd_file samples_file
    gender_file=$(mktemp)
    call_rate_file=$(mktemp)
    lrr_sd_file=$(mktemp)
    samples_file=$(mktemp)
    trap 'rm -f "${gender_file}" "${call_rate_file}" "${lrr_sd_file}" "${samples_file}"' RETURN

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
    # Compute call rate per sample using matrix output (single BCF pass).
    #
    # Instead of emitting one line per sample per variant (O(S×V) lines),
    # output one tab-delimited row per variant with all samples' GT values.
    # awk processes each row and accumulates per-column (per-sample) counts.
    # For 5000 samples × 2.5M variants this reduces pipe volume ~5000×.
    # -----------------------------------------------------------------
    echo "  Computing call rate per sample (${n_samples_vcf} samples, ${n_genotypeable} variants)..."
    bcftools view -e 'INFO/INTENSITY_ONLY=1' "${vcf_file}" 2>/dev/null | \
    bcftools query -f '[\t%GT]\n' 2>/dev/null | \
        awk -F'\t' '
        {
            for (i = 2; i <= NF; i++) {
                total[i]++
                if ($i != "./." && $i != "." && $i != ".|.") called[i]++
            }
        }
        END {
            for (i = 2; i <= NF; i++) {
                idx = i - 2
                cr = (total[i] > 0) ? called[i] / total[i] : 0
                print idx "\t" cr
            }
        }' > "${call_rate_file}.idx" 2>/dev/null || true

    # Map column indices back to sample names
    awk 'NR==FNR { names[NR-1] = $0; next }
         { print names[$1+0] "\t" $2 }' \
        "${samples_file}" "${call_rate_file}.idx" | \
        sort -k1,1 > "${call_rate_file}"
    rm -f "${call_rate_file}.idx"

    if [[ -s "${call_rate_file}" ]]; then
        n_cr=$(wc -l < "${call_rate_file}" | tr -d ' ')
        echo "  [diag] Call rate file: ${n_cr} samples"
        echo "  [diag] Call rate first 3 lines:"
        head -3 "${call_rate_file}" | sed 's/^/    [diag]   /'
    else
        echo "  [diag] WARNING: Call rate file is empty"
    fi

    # -----------------------------------------------------------------
    # Compute LRR standard deviation per sample (autosomes only).
    #
    # Same matrix approach: one row per variant, all samples' LRR values
    # tab-delimited. awk accumulates running sum and sum-of-squares per
    # column for online variance computation.
    #
    # Excludes sex chromosomes and MT to avoid inflation from hemizygosity.
    # Filters non-numeric LRR values (nan, inf) that gtc2vcf can produce.
    # -----------------------------------------------------------------
    echo "  Computing LRR standard deviation per sample..."
    bcftools view -e 'INFO/INTENSITY_ONLY=1' -t ^chrX,chrY,chrM,X,Y,MT "${vcf_file}" 2>/dev/null | \
    bcftools query -f '[\t%LRR]\n' 2>/dev/null | \
        awk -F'\t' '
        {
            for (i = 2; i <= NF; i++) {
                v = $i
                if (v != "." && v != "" && v ~ /^-?[0-9]/) {
                    n[i]++
                    sum[i] += v
                    sum2[i] += v * v
                }
            }
        }
        END {
            for (i = 2; i <= NF; i++) {
                idx = i - 2
                if (n[i] > 1) {
                    mean = sum[i] / n[i]
                    var = (sum2[i] / n[i]) - (mean * mean)
                    if (var < 0) var = 0
                    print idx "\t" sqrt(var)
                } else {
                    print idx "\tNA"
                }
            }
        }' > "${lrr_sd_file}.idx" 2>/dev/null || true

    # Map column indices back to sample names
    awk 'NR==FNR { names[NR-1] = $0; next }
         { print names[$1+0] "\t" $2 }' \
        "${samples_file}" "${lrr_sd_file}.idx" | \
        sort -k1,1 > "${lrr_sd_file}"
    rm -f "${lrr_sd_file}.idx"

    if [[ -s "${lrr_sd_file}" ]]; then
        n_sd=$(wc -l < "${lrr_sd_file}" | tr -d ' ')
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
