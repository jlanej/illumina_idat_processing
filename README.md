# Illumina IDAT Processing Pipeline

Process raw Illumina array IDAT data on HPC or command line **without Genome Studio**.

Inspired by [MoChA](https://github.com/freeseek/mocha) and [MoChA WDL](https://github.com/freeseek/mochawdl), this pipeline uses the open-source [bcftools](https://github.com/samtools/bcftools) plugins [idat2gtc](https://github.com/freeseek/gtc2vcf) and [gtc2vcf](https://github.com/freeseek/gtc2vcf) to convert raw Illumina intensity data directly to VCF format with BAF (B Allele Frequency) and LRR (Log R Ratio) intensities.

## Overview

The pipeline runs in three phases:

1. **Manifest Realignment**: Probe flank sequences from the CSV manifest are realigned against the reference genome using BWA, following the [MoChA/gtc2vcf](https://github.com/freeseek/gtc2vcf) best practice. This validates probe positions and enables correct coordinate mapping regardless of the original manifest build — **even if the CSV claims the same build as the reference**, this step catches discrepancies and ensures probes are placed correctly. A detailed mapping summary is output showing mapped, unmapped, multi-mapped, and coordinate-changed probes.

2. **Stage 1 — Initial Genotyping**: Convert IDAT files to GTC (genotype calls) using the Illumina GenCall algorithm (via `idat2gtc`), then to VCF format. Compute per-sample QC metrics including **call rate** and **LRR standard deviation**.

3. **Stage 2 — Reclustering**: Identify high-quality samples from Stage 1 QC metrics, recompute the EGT genotype cluster file from those samples' intensity data, then re-call genotypes for all samples using the new study-specific clusters. The `--adjust-clusters` option is also applied during VCF conversion for additional BAF/LRR correction.

The output VCF (with `GT`, `BAF`, and `LRR` FORMAT fields) is ready for downstream analysis such as phasing with [SHAPEIT5](https://odelaneau.github.io/shapeit5/), mosaic chromosomal alteration detection with [MoChA](https://github.com/freeseek/mocha), or imputation with [IMPUTE5](https://jmarchini.org/software/#impute-5).

## Quick Start

### Using Docker / Apptainer (recommended for HPC)

A pre-built container image is available from GitHub Container Registry:

```bash
# Docker
docker pull ghcr.io/jlanej/illumina_idat_processing:main
docker run -v /path/to/data:/data ghcr.io/jlanej/illumina_idat_processing:main \
    /opt/scripts/run_pipeline.sh \
    --idat-dir /data/idats --array-name GSA-24v3-0_A1 --output-dir /data/output

# Apptainer / Singularity (HPC)
apptainer pull docker://ghcr.io/jlanej/illumina_idat_processing:main
apptainer exec --bind /path/to/data:/data illumina_idat_processing_main.sif \
    bash /opt/scripts/run_pipeline.sh \
    --idat-dir /data/idats --array-name GSA-24v3-0_A1 --output-dir /data/output
```

### 1000 Genomes Example

Process 1000 Genomes Omni2.5 IDAT data (downloads automatically):

```bash
# Process 200 samples (quick test)
bash scripts/process_1000g.sh --output-dir ./1000g_output --num-samples 200

# Process all ~2141 samples on HPC with Apptainer
apptainer exec --bind $PWD illumina_idat_processing_main.sif \
    bash /opt/scripts/process_1000g.sh \
    --output-dir $PWD/1000g_output --num-samples all --threads 8
```

### From Source

```bash
# 1. Install dependencies (bcftools + plugins)
bash scripts/install_dependencies.sh --prefix $HOME/bin
export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"

# 2. Run the pipeline (minimal arguments - auto-downloads manifests and reference)
bash scripts/run_pipeline.sh \
    --idat-dir /path/to/idat/files \
    --array-name GSA-24v3-0_A1 \
    --output-dir /path/to/output
```

The pipeline will:
- Download the correct BPM, EGT, and CSV manifest files for the specified array
- Download the GRCh38 reference genome
- Run Stage 1 (initial genotyping + QC)
- Run Stage 2 (recluster on high-quality samples + reprocess)

## Requirements

- Linux (tested on Ubuntu/Debian, should work on any HPC)
- GCC 5+ (for compiling bcftools and plugins)
- Python 3.6+ with NumPy (for EGT reclustering in Stage 2)
- Standard tools: `wget`, `samtools`, `make`
- ~30 GB disk for GRCh38 reference genome

All other dependencies (bcftools, gtc2vcf, idat2gtc) are installed by `scripts/install_dependencies.sh`.

## Usage

### Full Pipeline

```bash
bash scripts/run_pipeline.sh [OPTIONS]
```

**Required:**
| Option | Description |
|--------|-------------|
| `--idat-dir DIR` | Directory containing IDAT files |
| `--output-dir DIR` | Output directory |

**Manifest options (one of):**
| Option | Description |
|--------|-------------|
| `--array-name NAME` | Known array name (auto-downloads manifests) |
| `--manifest-dir DIR` | Directory containing BPM, EGT, and CSV files |
| `--bpm FILE` / `--egt FILE` / `--csv FILE` | Individual manifest files |

**Reference options:**
| Option | Description |
|--------|-------------|
| `--ref-fasta FILE` | Reference FASTA (overrides auto-download) |
| `--ref-dir DIR` | Directory with reference genome |
| `--genome NAME` | Genome build: `GRCh37` or `GRCh38` (default: `GRCh38`) |

**Processing options:**
| Option | Description |
|--------|-------------|
| `--threads INT` | Number of threads (default: 1) |
| `--min-call-rate FLOAT` | Min call rate for HQ samples in Stage 2 (default: 0.97) |
| `--max-lrr-sd FLOAT` | Max LRR SD for HQ samples in Stage 2 (default: 0.35) |
| `--skip-stage2` | Skip Stage 2 reclustering |
| `--skip-download` | Do not auto-download manifests or reference |

### Individual Scripts

Each stage can also be run independently:

```bash
# Install dependencies
bash scripts/install_dependencies.sh --prefix /path/to/install

# Download manifests for a known array
bash scripts/download_manifests.sh --array-name GSA-24v3-0_A1 --output-dir manifests/

# Download reference genome
bash scripts/download_reference.sh --genome GRCh38 --output-dir reference/

# Stage 1: Initial genotyping
bash scripts/stage1_initial_genotyping.sh \
    --idat-dir /path/to/idats \
    --bpm manifests/GSA-24v3-0_A1.bpm \
    --egt manifests/GSA-24v3-0_A1_ClusterFile.egt \
    --csv manifests/GSA-24v3-0_A1.csv \
    --ref-fasta reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --output-dir output/stage1

# Stage 2: Recluster and reprocess
bash scripts/stage2_recluster.sh \
    --idat-dir /path/to/idats \
    --bpm manifests/GSA-24v3-0_A1.bpm \
    --egt manifests/GSA-24v3-0_A1_ClusterFile.egt \
    --csv manifests/GSA-24v3-0_A1.csv \
    --ref-fasta reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --stage1-dir output/stage1 \
    --output-dir output/stage2
```

### WDL Workflow

An optional [WDL](https://openwdl.org/) workflow is provided at `wdl/illumina_idat_processing.wdl` for use with [Cromwell](https://cromwell.readthedocs.io/) or [Terra](https://terra.bio/):

```json
{
  "illumina_idat_processing.sample_set_id": "my_cohort",
  "illumina_idat_processing.green_idat_files": ["gs://bucket/sample1_Grn.idat"],
  "illumina_idat_processing.red_idat_files": ["gs://bucket/sample1_Red.idat"],
  "illumina_idat_processing.sample_ids": ["sample1"],
  "illumina_idat_processing.bpm_file": "gs://bucket/manifests/array.bpm",
  "illumina_idat_processing.egt_file": "gs://bucket/manifests/array.egt",
  "illumina_idat_processing.csv_file": "gs://bucket/manifests/array.csv",
  "illumina_idat_processing.ref_fasta": "gs://bucket/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
  "illumina_idat_processing.ref_fasta_fai": "gs://bucket/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"
}
```

## Supported Arrays

The following arrays have pre-configured manifest download URLs in `config/manifest_urls.tsv`:

| Array | Description |
|-------|-------------|
| `GSA-24v3-0_A1` | Global Screening Array v3.0 (GRCh37) |
| `GSA-24v3-0_A2` | Global Screening Array v3.0 (GRCh38) |
| `GSA-24v2-0_A1` | Global Screening Array v2.0 (GRCh37) |
| `GSA-24v1-0_A1` | Global Screening Array v1.0 (GRCh37) |
| `MEGA_Consortium_v2-0_A1` | Multi-Ethnic Global Array v2.0 |
| `InfiniumOmni2-5-8v1-5_A1` | Infinium Omni2.5 v1.5 |
| `HumanCNV370v1_C` | HumanCNV370-Duo v1.0 |
| `HumanOmni2-5-4v1-Multi_B` | HumanOmni2.5-4 v1 Multi (1000G) |

Additional arrays can be processed by providing manifest files directly with `--bpm`, `--egt`, and `--csv`, or by adding entries to `config/manifest_urls.tsv`.

## Output Files

```
output/
├── stage1/
│   ├── gtc/                    # GTC genotype files (original clusters)
│   ├── vcf/
│   │   └── stage1_initial.bcf  # Initial VCF with BAF + LRR
│   └── qc/
│       ├── stage1_sample_qc.tsv  # Per-sample call rate and LRR SD
│       └── gtc_metadata.tsv      # GTC file metadata
├── stage2/
│   ├── gtc/                    # Re-called GTC files (new clusters)
│   ├── vcf/
│   │   └── stage2_reclustered.bcf  # Reclustered VCF (final output)
│   ├── qc/
│   │   ├── stage2_sample_qc.tsv    # Updated QC metrics
│   │   ├── high_quality_samples.txt # Samples used for reclustering
│   │   └── excluded_samples.txt     # Low-quality samples
│   └── clusters/
│       └── reclustered.egt          # Study-specific EGT cluster file
├── manifests/                  # Downloaded manifest files (if auto-downloaded)
├── realigned_manifests/        # Realigned manifest and mapping summary
│   ├── *.realigned.csv         # CSV with validated/remapped coordinates
│   └── realign_summary.txt     # Detailed mapping statistics
└── reference/                  # Downloaded reference (if auto-downloaded)
```

### QC Metrics File Format

The `*_sample_qc.tsv` files contain:

| Column | Description |
|--------|-------------|
| `sample_id` | Sample identifier |
| `call_rate` | Fraction of variants with non-missing genotype calls |
| `lrr_sd` | Standard deviation of Log R Ratio (lower is better) |
| `computed_gender` | Inferred gender (1=male, 2=female) |

## Processing Rationale

### Manifest Realignment

Before genotyping, the pipeline always realigns probe flank sequences from the CSV manifest against the reference genome using BWA (`scripts/realign_manifest.sh`). This follows the [MoChA WDL](https://github.com/freeseek/mochawdl) best practice and is critical for two reasons:

1. **Cross-build remapping**: If the manifest was designed for GRCh37 but processing targets GRCh38 (or vice versa), `gtc2vcf` does **not** automatically remap coordinates. The flank alignment produces a remapped CSV with correct positions for the target reference.

2. **Position validation**: Even when the manifest claims the same build, probes can have incorrect positions due to assembly updates, manufacturer errors, or custom array designs. Realignment validates every probe against the actual reference sequence.

The pipeline outputs a detailed mapping summary (`realign_summary.txt`) reporting:
- Total probes, mapped/unmapped counts, mapping quality distribution
- Coordinate changes: positions confirmed, changed, newly placed, or lost

### Two-Stage Genotyping

Standard Illumina genotyping uses pre-built cluster files (EGT) derived from a training set. However, batch effects, population differences, and varying sample quality can degrade genotype calls when the default clusters don't match the study population.

**Stage 1** applies the default EGT clusters to produce initial genotype calls and computes per-sample QC metrics. Samples with low call rate (< 0.97) or high LRR standard deviation (> 0.35) are flagged as low quality.

**Stage 2** performs true reclustering in three steps:

1. **Recompute EGT**: Using `scripts/recluster_egt.py`, extracts normalized intensity data (theta, R coordinates) from Stage 1 GTC files of high-quality samples, grouped by their genotype calls (AA/AB/BB). For each probe, recomputes the cluster means and standard deviations in (theta, R) space. Probes with insufficient data fall back to the original cluster definitions. The result is a new EGT file tailored to the study population.

2. **Re-call genotypes**: Runs `idat2gtc` with the new EGT file, producing improved genotype calls for ALL samples. Because the GenCall algorithm uses the EGT cluster definitions to assign each sample's intensities to the correct genotype, study-specific clusters yield more accurate genotype calls than the default Illumina training set.

3. **Convert to VCF**: Runs `gtc2vcf` with `--adjust-clusters` for additional median adjustment of BAF and LRR values.

This approach is fundamentally different from just using `--adjust-clusters` in `gtc2vcf` (which only adjusts BAF/LRR post-hoc without changing the underlying genotype calls). By recomputing the actual EGT cluster definitions, we improve the genotype calls themselves — analogous to the iterative clustering performed in Genome Studio, but fully automated for HPC.

## Downstream Analysis

After running this pipeline, the output VCF can be used directly with:

- **[MoChA](https://github.com/freeseek/mocha)** — Detection of mosaic chromosomal alterations
- **[SHAPEIT5](https://odelaneau.github.io/shapeit5/)** — Genotype phasing
- **[IMPUTE5](https://jmarchini.org/software/)** / **[Beagle5](https://faculty.washington.edu/browning/beagle/)** — Genotype imputation
- **[PLINK2](https://www.cog-genomics.org/plink/2.0/)** — GWAS and population genetics

### Copy Number Variation (CNV) Calling

The pipeline's BAF and LRR output is directly suitable for CNV detection. See **[docs/cnv_calling_methods.md](docs/cnv_calling_methods.md)** for a detailed survey of compatible methods, including PennCNV, bcftools +mocha (already installed), QuantiSNP, and EnsembleCNV.

## References

- Loh P., Genovese G., McCarroll S., Price A. et al. *Insights about clonal expansions from 8,342 mosaic chromosomal alterations.* Nature 559, 350–355 (2018). [DOI: 10.1038/s41586-018-0321-x](https://doi.org/10.1038/s41586-018-0321-x)
- Liu A. et al. *Genetic drivers and cellular selection of female mosaic X chromosome loss.* Nature (2024). [DOI: 10.1038/s41586-024-07533-7](https://doi.org/10.1038/s41586-024-07533-7)
- [gtc2vcf](https://github.com/freeseek/gtc2vcf) — Convert Illumina/Affymetrix intensity data to VCF without Windows
- [MoChA WDL](https://github.com/freeseek/mochawdl) — WDL pipelines for the full MoChA workflow

## License

MIT License. See [LICENSE](LICENSE).
