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

    call gtc_to_vcf as stage2_gtc2vcf {
      input:
        gtc_files = stage1_idat2gtc.gtc_files,
        bpm_file = bpm_file,
        csv_file = csv_file,
        egt_file = egt_file,
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
    Array[File] gtc_files = stage1_idat2gtc.gtc_files
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

    # Link files to working directory
    for f in "${green_files[@]}" "${red_files[@]}"; do
      ln -s "$f" .
    done

    bcftools +idat2gtc \
      --bpm ~{bpm_file} \
      --egt ~{egt_file} \
      --output gtc_output \
      "${green_files[@]}"

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
    bcftools sort -Ou -T ./bcftools. | \
    bcftools norm --no-version -Ob -c x -f ~{ref_fasta} \
      -o ~{filebase}.bcf --write-index
  >>>

  output {
    File vcf_file = "~{filebase}.bcf"
    File vcf_idx = "~{filebase}.bcf.csi"
    File metadata_tsv = "~{filebase}.metadata.tsv"
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

    # Compute call rate and LRR standard deviation per sample
    echo -e "sample_id\tcall_rate\tlrr_sd" > ~{filebase}.qc.tsv

    # Get sample list
    bcftools query -l ~{vcf_file} > samples.txt

    # Total variants
    n_total=$(bcftools view -H ~{vcf_file} | wc -l)

    # Per-sample call rate: count non-missing genotypes
    while read -r sample; do
      n_called=$(bcftools view -s "${sample}" -H ~{vcf_file} | \
        bcftools query -f '[%GT]\n' | grep -cv '^\./\.' || true)
      if [[ "${n_total}" -gt 0 ]]; then
        call_rate=$(awk "BEGIN {printf \"%.6f\", ${n_called}/${n_total}}")
      else
        call_rate="NA"
      fi

      # LRR standard deviation
      lrr_sd=$(bcftools view -s "${sample}" ~{vcf_file} | \
        bcftools query -f '[%LRR]\n' | \
        awk '$1 != "." && $1 != "" {
          n++; sum += $1; sum2 += $1*$1
        } END {
          if (n > 1) {
            mean = sum/n; var = sum2/n - mean*mean
            if (var < 0) var = 0
            printf "%.6f", sqrt(var)
          } else print "NA"
        }')

      echo -e "${sample}\t${call_rate}\t${lrr_sd}"
    done < samples.txt >> ~{filebase}.qc.tsv
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
