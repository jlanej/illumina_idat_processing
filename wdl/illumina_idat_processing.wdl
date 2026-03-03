version 1.0

## Illumina IDAT Processing Pipeline
##
## Runs the two-stage Illumina IDAT processing pipeline by calling
## run_pipeline.sh inside the pre-built container image:
##   ghcr.io/jlanej/illumina_idat_processing
##
## For a quickstart with Cromwell and Apptainer, run:
##   bash wdl/quickstart.sh

workflow illumina_idat_processing {
  input {
    # All IDAT files (green and red, all in one array)
    Array[File] idat_files

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
    Int     threads     = 4
    Boolean skip_stage2 = false

    String docker = "ghcr.io/jlanej/illumina_idat_processing:main"
  }

  call run_pipeline {
    input:
      idat_files    = idat_files,
      bpm_file      = bpm_file,
      egt_file      = egt_file,
      csv_file      = csv_file,
      ref_fasta     = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      min_call_rate = min_call_rate,
      max_lrr_sd    = max_lrr_sd,
      threads       = threads,
      skip_stage2   = skip_stage2,
      docker        = docker,
  }

  output {
    File  final_vcf            = run_pipeline.final_vcf
    File  final_vcf_idx        = run_pipeline.final_vcf_idx
    File  final_qc_tsv         = run_pipeline.final_qc_tsv
    File  stage1_qc_tsv        = run_pipeline.stage1_qc_tsv
    File? compiled_sample_sheet = run_pipeline.compiled_sample_sheet
    File? qc_diagnostic_report = run_pipeline.qc_diagnostic_report
  }

  meta {
    description: "Two-stage Illumina IDAT processing pipeline. Calls run_pipeline.sh inside the pre-built container image ghcr.io/jlanej/illumina_idat_processing."
  }
}

# ============================================================
# Task: Run the complete pipeline via run_pipeline.sh
# ============================================================
task run_pipeline {
  input {
    Array[File] idat_files
    File bpm_file
    File egt_file
    File csv_file
    File ref_fasta
    File ref_fasta_fai
    Float   min_call_rate
    Float   max_lrr_sd
    Int     threads
    Boolean skip_stage2

    String docker
    Int    cpu         = 8
    Int    disk_size   = 500
    Float  memory      = 32.0
    Int    preemptible = 1
    Int    maxRetries  = 0
  }

  command <<<
    set -euo pipefail

    # Collect all IDAT files into a flat directory for run_pipeline.sh
    mkdir -p idat_dir
    for f in ~{sep=' ' idat_files}; do
      ln -sf "$f" idat_dir/
    done

    mkdir -p pipeline_output

    bash /opt/scripts/run_pipeline.sh \
      --idat-dir      idat_dir \
      --bpm           ~{bpm_file} \
      --egt           ~{egt_file} \
      --csv           ~{csv_file} \
      --ref-fasta     ~{ref_fasta} \
      --output-dir    pipeline_output \
      --threads       ~{threads} \
      --min-call-rate ~{min_call_rate} \
      --max-lrr-sd    ~{max_lrr_sd} \
      --skip-download \
      ~{if skip_stage2 then "--skip-stage2" else ""}
  >>>

  output {
    File final_vcf = if skip_stage2
      then "pipeline_output/stage1/vcf/stage1_initial.bcf"
      else "pipeline_output/stage2/vcf/stage2_reclustered.bcf"
    File final_vcf_idx = if skip_stage2
      then "pipeline_output/stage1/vcf/stage1_initial.bcf.csi"
      else "pipeline_output/stage2/vcf/stage2_reclustered.bcf.csi"
    File final_qc_tsv = if skip_stage2
      then "pipeline_output/stage1/qc/stage1_sample_qc.tsv"
      else "pipeline_output/stage2/qc/stage2_sample_qc.tsv"
    File  stage1_qc_tsv        = "pipeline_output/stage1/qc/stage1_sample_qc.tsv"
    File? compiled_sample_sheet = "pipeline_output/compiled_sample_sheet.tsv"
    File? qc_diagnostic_report = "pipeline_output/qc_diagnostic_report.txt"
  }

  runtime {
    docker:      docker
    cpu:         cpu
    disks:       "local-disk " + disk_size + " HDD"
    memory:      memory + " GiB"
    preemptible: preemptible
    maxRetries:  maxRetries
  }
}
