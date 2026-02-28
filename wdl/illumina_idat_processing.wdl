version 1.0

## Illumina IDAT Processing Pipeline
##
## Two-stage processing of raw Illumina array IDAT data:
##   Stage 1: Initial genotyping (IDAT → GTC → VCF) with per-sample QC
##   Stage 2: Recluster on high-quality samples, reprocess all samples
##
## Inspired by mocha (https://github.com/freeseek/mocha) and mochawdl
## (https://github.com/freeseek/mochawdl).

workflow illumina_idat_processing {
  input {
    String sample_set_id
    Array[File] green_idat_files
    Array[File] red_idat_files
    Array[String] sample_ids

    # Manifest files
    File bpm_file
    File egt_file
    File csv_file

    # Reference genome
    File ref_fasta
    File ref_fasta_fai

    # QC thresholds for Stage 2
    Float min_call_rate = 0.97
    Float max_lrr_sd = 0.35

    # Processing options
    Int threads = 4
    Boolean skip_stage2 = false
    String docker = "us.gcr.io/mccarroll-mocha/bcftools:1.21-20250101"
  }

  # Stage 1: Initial genotyping
  call idat_to_gtc as stage1_idat2gtc {
    input:
      green_idat_files = green_idat_files,
      red_idat_files = red_idat_files,
      bpm_file = bpm_file,
      egt_file = egt_file,
      filebase = sample_set_id + ".stage1",
      docker = docker,
  }

  call gtc_to_vcf as stage1_gtc2vcf {
    input:
      gtc_files = stage1_idat2gtc.gtc_files,
      bpm_file = bpm_file,
      csv_file = csv_file,
      egt_file = egt_file,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      sample_ids = sample_ids,
      adjust_clusters = false,
      threads = threads,
      filebase = sample_set_id + ".stage1",
      docker = docker,
  }

  call compute_qc_metrics as stage1_qc {
    input:
      vcf_file = stage1_gtc2vcf.vcf_file,
      vcf_idx = stage1_gtc2vcf.vcf_idx,
      filebase = sample_set_id + ".stage1",
      docker = docker,
  }

  # Stage 2: Recluster on high-quality samples and reprocess
  if (!skip_stage2) {
    call select_hq_samples {
      input:
        qc_tsv = stage1_qc.qc_tsv,
        min_call_rate = min_call_rate,
        max_lrr_sd = max_lrr_sd,
        filebase = sample_set_id,
        docker = docker,
    }

    call recluster_egt {
      input:
        original_egt = egt_file,
        bpm_file = bpm_file,
        gtc_files = stage1_idat2gtc.gtc_files,
        hq_samples = select_hq_samples.hq_samples,
        filebase = sample_set_id,
        docker = docker,
    }

    call idat_to_gtc as stage2_idat2gtc {
      input:
        green_idat_files = green_idat_files,
        red_idat_files = red_idat_files,
        bpm_file = bpm_file,
        egt_file = recluster_egt.reclustered_egt,
        filebase = sample_set_id + ".stage2",
        docker = docker,
    }

    call gtc_to_vcf as stage2_gtc2vcf {
      input:
        gtc_files = stage2_idat2gtc.gtc_files,
        bpm_file = bpm_file,
        csv_file = csv_file,
        egt_file = recluster_egt.reclustered_egt,
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        sample_ids = sample_ids,
        adjust_clusters = true,
        threads = threads,
        filebase = sample_set_id + ".stage2",
        docker = docker,
    }

    call compute_qc_metrics as stage2_qc {
      input:
        vcf_file = stage2_gtc2vcf.vcf_file,
        vcf_idx = stage2_gtc2vcf.vcf_idx,
        filebase = sample_set_id + ".stage2",
        docker = docker,
    }
  }

  output {
    File final_vcf = select_first([stage2_gtc2vcf.vcf_file, stage1_gtc2vcf.vcf_file])
    File final_vcf_idx = select_first([stage2_gtc2vcf.vcf_idx, stage1_gtc2vcf.vcf_idx])
    File final_qc_tsv = select_first([stage2_qc.qc_tsv, stage1_qc.qc_tsv])
    File stage1_qc_tsv = stage1_qc.qc_tsv
    Array[File] gtc_files = select_first([stage2_idat2gtc.gtc_files, stage1_idat2gtc.gtc_files])
    File? reclustered_egt_file = recluster_egt.reclustered_egt
    File? hq_samples_file = select_hq_samples.hq_samples
    File? excluded_samples_file = select_hq_samples.excluded_samples
  }

  meta {
    description: "Two-stage Illumina IDAT processing pipeline: initial genotyping with QC, followed by reclustering on high-quality samples."
  }
}

# ============================================================
# Task: Convert IDAT files to GTC files
# ============================================================
task idat_to_gtc {
  input {
    Array[File] green_idat_files
    Array[File] red_idat_files
    File bpm_file
    File egt_file
    String filebase

    String docker
    Int cpu = 1
    Int disk_size = 50
    Float memory = 8.0
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mkdir -p gtc_output

    green_files=(~{sep=' ' green_idat_files})
    red_files=(~{sep=' ' red_idat_files})

    # Build interleaved green/red IDAT pairs.
    # idat2gtc expects pairs: <green.idat> <red.idat> [<green.idat> <red.idat> ...]
    idat_args=()
    for (( i=0; i<${#green_files[@]}; i++ )); do
      idat_args+=("${green_files[i]}" "${red_files[i]}")
    done

    bcftools +idat2gtc \
      --bpm ~{bpm_file} \
      --egt ~{egt_file} \
      --output gtc_output \
      "${idat_args[@]}"

    ls gtc_output/*.gtc > gtc_list.txt
  >>>

  output {
    Array[File] gtc_files = read_lines("gtc_list.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GiB"
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

# ============================================================
# Task: Convert GTC files to VCF
# ============================================================
task gtc_to_vcf {
  input {
    Array[File] gtc_files
    File bpm_file
    File csv_file
    File egt_file
    File ref_fasta
    File ref_fasta_fai
    Array[String] sample_ids
    Boolean adjust_clusters
    Int threads
    String filebase

    String docker
    Int cpu = 4
    Int disk_size = 100
    Float memory = 16.0
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mkdir -p gtc_dir

    # Link GTC files
    for f in ~{sep=' ' gtc_files}; do
      ln -s "$f" gtc_dir/
    done

    bcftools +gtc2vcf \
      --no-version -Ou \
      --bpm ~{bpm_file} \
      --csv ~{csv_file} \
      --egt ~{egt_file} \
      --fasta-ref ~{ref_fasta} \
      --gtcs gtc_dir \
      --extra ~{filebase}.metadata.tsv \
      ~{if adjust_clusters then "--adjust-clusters" else ""} \
      --threads ~{threads} | \
    bcftools sort -Ob -T ./bcftools. \
      -o ~{filebase}.pre_norm.bcf --write-index

    # Diagnostic: count variants before normalization
    N_PRE_NORM=$(bcftools index -n ~{filebase}.pre_norm.bcf 2>/dev/null || \
        bcftools view -H ~{filebase}.pre_norm.bcf | wc -l)
    echo "Variants before normalization: ${N_PRE_NORM}"

    # Normalize with -c ws (warn and swap REF/ALT to match reference).
    # CRITICAL: Do NOT use -c x. When processing cross-build data,
    # -c x silently sets genotypes to missing for REF-mismatched probes,
    # deflating call rates. -c ws swaps REF/ALT to match the reference.
    bcftools norm --no-version -Ob -c ws -f ~{ref_fasta} \
      ~{filebase}.pre_norm.bcf \
      -o ~{filebase}.bcf --write-index 2>~{filebase}.norm_warnings.log

    N_POST_NORM=$(bcftools index -n ~{filebase}.bcf 2>/dev/null || \
        bcftools view -H ~{filebase}.bcf | wc -l)
    N_REF_SWAPS=$(grep -c "REF_MISMATCH" ~{filebase}.norm_warnings.log || true)
    echo "Variants after normalization: ${N_POST_NORM}"
    echo "REF/ALT swaps: ${N_REF_SWAPS}"

    rm -f ~{filebase}.pre_norm.bcf ~{filebase}.pre_norm.bcf.csi
  >>>

  output {
    File vcf_file = "~{filebase}.bcf"
    File vcf_idx = "~{filebase}.bcf.csi"
    File metadata_tsv = "~{filebase}.metadata.tsv"
    File norm_warnings_log = "~{filebase}.norm_warnings.log"
  }

  runtime {
    docker: docker
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GiB"
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

# ============================================================
# Task: Compute per-sample QC metrics (call rate, LRR SD)
# ============================================================
task compute_qc_metrics {
  input {
    File vcf_file
    File vcf_idx
    String filebase

    String docker
    Int cpu = 1
    Int disk_size = 50
    Float memory = 8.0
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail

    # Compute call rate and LRR standard deviation per sample.
    # Uses single-pass matrix output (one row per variant, all samples'
    # values tab-delimited) instead of per-sample BCF scans.

    # Get sample list
    bcftools query -l ~{vcf_file} > samples.txt
    n_samples=$(wc -l < samples.txt)

    # Call rate: one row per variant, tab-delimited GT for all samples.
    # awk accumulates per-column (per-sample) called/total counts.
    bcftools view -e 'INFO/INTENSITY_ONLY=1' ~{vcf_file} | \
    bcftools query -f '[\t%GT]\n' | \
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
            printf "%d\t%.6f\n", idx, cr
        }
    }' > call_rates.idx.tsv

    # LRR SD on autosomes: same matrix approach, online variance computation.
    bcftools view -e 'INFO/INTENSITY_ONLY=1' -t ^chrX,chrY,chrM,X,Y,MT ~{vcf_file} | \
    bcftools query -f '[\t%LRR]\n' | \
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
                printf "%d\t%.6f\n", idx, sqrt(var)
            } else {
                printf "%d\tNA\n", idx
            }
        }
    }' > lrr_sd.idx.tsv

    # Merge: map column indices to sample names and combine
    echo -e "sample_id\tcall_rate\tlrr_sd" > ~{filebase}.qc.tsv
    paste call_rates.idx.tsv lrr_sd.idx.tsv | \
    awk -F'\t' 'NR==FNR { names[NR-1] = $0; next }
                { print names[$1+0] "\t" $2 "\t" $4 }' \
        samples.txt - >> ~{filebase}.qc.tsv
  >>>

  output {
    File qc_tsv = "~{filebase}.qc.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GiB"
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

# ============================================================
# Task: Select high-quality samples based on QC thresholds
# ============================================================
task select_hq_samples {
  input {
    File qc_tsv
    Float min_call_rate
    Float max_lrr_sd
    String filebase

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    awk -F'\t' -v min_cr=~{min_call_rate} -v max_sd=~{max_lrr_sd} '
      NR == 1 { next }
      $2+0 >= min_cr+0 && ($3 == "NA" || $3+0 <= max_sd+0) { print $1 }
    ' ~{qc_tsv} > ~{filebase}.hq_samples.txt

    awk -F'\t' -v min_cr=~{min_call_rate} -v max_sd=~{max_lrr_sd} '
      NR == 1 { next }
      !($2+0 >= min_cr+0 && ($3 == "NA" || $3+0 <= max_sd+0)) { print $1 }
    ' ~{qc_tsv} > ~{filebase}.excluded_samples.txt

    n_hq=$(wc -l < ~{filebase}.hq_samples.txt)
    n_excl=$(wc -l < ~{filebase}.excluded_samples.txt)
    echo "High-quality samples: ${n_hq}" >&2
    echo "Excluded samples: ${n_excl}" >&2
  >>>

  output {
    File hq_samples = "~{filebase}.hq_samples.txt"
    File excluded_samples = "~{filebase}.excluded_samples.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GiB"
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

# ============================================================
# Task: Recompute EGT cluster file from high-quality samples
# ============================================================
task recluster_egt {
  input {
    File original_egt
    File bpm_file
    Array[File] gtc_files
    File hq_samples
    String filebase

    Int min_cluster_samples = 5

    String docker
    Int cpu = 2
    Int disk_size = 50
    Float memory = 16.0
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mkdir -p gtc_dir

    # Link GTC files
    for f in ~{sep=' ' gtc_files}; do
      ln -s "$f" gtc_dir/
    done

    python3 /opt/recluster_egt.py \
      --egt ~{original_egt} \
      --bpm ~{bpm_file} \
      --gtc-dir gtc_dir \
      --hq-samples ~{hq_samples} \
      --output-egt ~{filebase}.reclustered.egt \
      --min-cluster-samples ~{min_cluster_samples}
  >>>

  output {
    File reclustered_egt = "~{filebase}.reclustered.egt"
  }

  runtime {
    docker: docker
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GiB"
    preemptible: preemptible
    maxRetries: maxRetries
  }
}
