# Alternative Run Methods

This document describes alternative ways to run the Illumina IDAT Processing Pipeline beyond the recommended Apptainer approach documented in the main [README](../README.md).

## Docker

If you have Docker installed (e.g., on a workstation or cloud VM), you can run the pipeline directly with the pre-built container image.

```bash
# Pull the image
docker pull ghcr.io/jlanej/illumina_idat_processing:main

# Run the full pipeline
docker run --rm \
    -v /path/to/data:/data \
    ghcr.io/jlanej/illumina_idat_processing:main \
    bash /opt/scripts/run_pipeline.sh \
    --idat-dir /data/idats \
    --array-name GSA-24v3-0_A1 \
    --output-dir /data/output

# With pre-downloaded manifests and reference
docker run --rm \
    -v /path/to/data:/data \
    ghcr.io/jlanej/illumina_idat_processing:main \
    bash /opt/scripts/run_pipeline.sh \
    --idat-dir /data/idats \
    --bpm /data/manifests/array.bpm \
    --egt /data/manifests/array.egt \
    --csv /data/manifests/array.csv \
    --ref-fasta /data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --output-dir /data/output
```

### 1000 Genomes Example (Docker)

```bash
docker run --rm \
    -v $PWD:/data \
    ghcr.io/jlanej/illumina_idat_processing:main \
    bash /opt/scripts/process_1000g.sh \
    --output-dir /data/1000g_output \
    --num-samples 200
```

## From Source

If you prefer to run the pipeline scripts directly without a container, you can install all dependencies from source.

### 1. Install dependencies

```bash
# Install bcftools with gtc2vcf/idat2gtc/mocha plugins
bash scripts/install_dependencies.sh --prefix $HOME/bin
export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"
```

### 2. Install Python packages

```bash
pip install numpy matplotlib
```

### 3. Run the pipeline

```bash
bash scripts/run_pipeline.sh \
    --idat-dir /path/to/idat/files \
    --array-name GSA-24v3-0_A1 \
    --output-dir /path/to/output
```

The pipeline will auto-download the correct BPM, EGT, and CSV manifest files for the specified array, as well as the GRCh38 reference genome.

### From-source requirements

- Linux (tested on Ubuntu/Debian; should work on any HPC)
- GCC 5+ (for compiling bcftools and plugins)
- Python 3.6+ with NumPy and Matplotlib
- Standard tools: `wget`, `samtools`, `bwa`, `make`
- ~30 GB disk for GRCh38 reference genome

All other dependencies (bcftools, gtc2vcf, idat2gtc) are installed by `scripts/install_dependencies.sh`.
