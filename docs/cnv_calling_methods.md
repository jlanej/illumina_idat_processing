# Copy Number Variation (CNV) Calling Methods

This document surveys CNV calling approaches that can operate on the BAF (B Allele Frequency) and LRR (Log R Ratio) intensities produced by this pipeline's output VCF. The methods below are ranked by their compatibility with our data product and ease of integration.

## Pipeline Data Product

The pipeline outputs a BCF/VCF with per-sample `GT`, `BAF`, and `LRR` FORMAT fields. These are the standard inputs for SNP-array-based CNV detection:

- **LRR** — Log R Ratio: normalized total intensity relative to the expected cluster position. Deviations indicate copy number changes (positive → gain, negative → loss).
- **BAF** — B Allele Frequency: relative proportion of the B allele signal. Heterozygous deletions shift BAF away from 0.5; LOH events collapse BAF to 0 or 1.

## Recommended Methods

### 1. PennCNV

- **Website**: <https://penncnv.openbioinformatics.org/>
- **Reference**: Wang K. et al. *PennCNV: An integrated hidden Markov model designed for high-resolution copy number variation detection in whole-genome SNP genotyping data.* Genome Research 17:1665–1674 (2007). [DOI: 10.1101/gr.6861907](https://doi.org/10.1101/gr.6861907)
- **Method**: Hidden Markov Model (HMM) that jointly models LRR and BAF along each chromosome to segment copy number states (0, 1, 2, 3, 4 copies) and detect LOH.
- **Input**: Per-SNP LRR and BAF values, plus a population frequency of B allele (PFB) file and GC content model for wave correction.
- **Integration**: **High**. PennCNV is the most widely used array-based CNV caller. The pipeline's per-probe LRR and BAF can be exported to PennCNV's tabular input format with a straightforward `bcftools query` extraction:
  ```bash
  bcftools query -f '%CHROM\t%POS\t%ID\t[%LRR]\t[%BAF]\n' stage2_reclustered.bcf
  ```
  PFB files can be computed from cohort BAF values. GC model files for GRCh38 are available from the PennCNV distribution.
- **Strengths**: Mature, well-validated, supports family-based calling (trio/joint), GC wave adjustment, and has extensive documentation.
- **Limitations**: Single-sample HMM; does not share information across samples within a cohort. Sensitive to LRR noise (benefits from Stage 2 reclustering).

### 2. bcftools +mocha (MoChA)

- **Website**: <https://github.com/freeseek/mocha>
- **Reference**: Loh P. et al. *Insights about clonal expansions from 8,342 mosaic chromosomal alterations.* Nature 559:350–355 (2018). [DOI: 10.1038/s41586-018-0321-x](https://doi.org/10.1038/s41586-018-0321-x)
- **Method**: Phased BAF analysis using an HMM to detect mosaic chromosomal alterations (mCAs), including mosaic gains, losses, and copy-neutral LOH. Operates on phased genotypes with BAF/LRR.
- **Input**: Phased VCF/BCF with `BAF` and `LRR` FORMAT fields — exactly what this pipeline produces.
- **Integration**: **Very high**. MoChA is already a bcftools plugin in our toolchain (installed by `install_dependencies.sh`). After phasing with SHAPEIT5 or Eagle, MoChA runs directly on the pipeline output:
  ```bash
  bcftools +mocha --input-stats --LRR-weight 0.2 \
      --call-rate 0.97 --bdev-LRR-BAF 6.0 \
      phased.bcf --output calls.tsv --mosaic-calls mosaic.tsv
  ```
- **Strengths**: Designed for the exact data type we produce. Detects mosaic events at low cell fractions. Already part of our dependency stack.
- **Limitations**: Primarily targets mosaic (sub-clonal) events rather than constitutional germline CNVs. Requires phased genotypes as input. Best suited for large cohorts.

### 3. QuantiSNP

- **Website**: <https://sites.google.com/site/quantisnp/>
- **Reference**: Colella S. et al. *QuantiSNP: an Objective Bayes Hidden-Markov Model to detect and accurately map copy number variation using SNP genotyping data.* Nucleic Acids Research 35:2013–2025 (2007). [DOI: 10.1093/nar/gkm076](https://doi.org/10.1093/nar/gkm076)
- **Method**: Objective Bayes HMM for CNV detection using LRR and BAF, with automatic parameter estimation (no need for manual training).
- **Input**: Per-SNP LRR and BAF values in tabular format.
- **Integration**: **Moderate**. Same extraction approach as PennCNV. QuantiSNP uses MATLAB Runtime, which may require additional setup on HPC systems.
- **Strengths**: Objective Bayes framework avoids manual tuning. Good detection of small CNVs.
- **Limitations**: Requires MATLAB Compiler Runtime. Less actively maintained than PennCNV. Single-sample calling.

### 4. EnsembleCNV

- **Website**: <https://github.com/ShenLab/EnsembleCNV>
- **Reference**: Zhang Z. et al. *EnsembleCNV: an ensemble machine learning algorithm to identify and genotype copy number variation using SNP array data.* Nucleic Acids Research 47:e39 (2019). [DOI: 10.1093/nar/gkz068](https://doi.org/10.1093/nar/gkz068)
- **Method**: Ensemble approach that combines calls from multiple individual CNV callers (PennCNV, QuantiSNP, iPattern) and uses a random forest classifier to produce consensus CNV calls with genotype assignments.
- **Input**: Raw intensity data (LRR/BAF) processed through multiple callers.
- **Integration**: **Moderate**. Requires running multiple callers first, then merging. Adds complexity but improves specificity.
- **Strengths**: Higher accuracy than any single caller. Provides CNV genotypes (copy number per sample per region) suitable for association testing.
- **Limitations**: Heavier computational requirement (runs multiple callers). More complex pipeline setup.

### 5. iPattern

- **Website**: Available via Illumina GenomeStudio or as standalone tool
- **Reference**: Pinto D. et al. *Comprehensive assessment of array-based platforms and calling algorithms for detection of copy number variants.* Nature Biotechnology 29:512–520 (2011). [DOI: 10.1038/nbt.1852](https://doi.org/10.1038/nbt.1852)
- **Method**: Pattern-recognition algorithm for CNV detection using LRR and BAF, with population-level frequency filtering.
- **Input**: LRR and BAF values from Illumina arrays.
- **Integration**: **Low-Moderate**. Typically integrated with GenomeStudio. Standalone usage requires Illumina-specific file formats.
- **Strengths**: Good for rare CNV detection. Accounts for population frequency to reduce false positives.
- **Limitations**: Tightly coupled to Illumina ecosystem. Less portable than open-source alternatives.

## Integration Summary

| Method | Input Match | Already Installed | Cohort-Aware | Ease of Integration |
|--------|------------|-------------------|-------------|-------------------|
| PennCNV | LRR + BAF | No | No (single-sample) | High |
| bcftools +mocha | BAF + LRR (phased VCF) | **Yes** | Yes | Very High |
| QuantiSNP | LRR + BAF | No | No (single-sample) | Moderate |
| EnsembleCNV | LRR + BAF (multi-caller) | No | Yes (ensemble) | Moderate |
| iPattern | LRR + BAF | No | No | Low-Moderate |

## Recommended Integration Path

1. **Immediate**: Use `bcftools +mocha` for mosaic chromosomal alteration detection — it is already installed and operates directly on the pipeline output after phasing.

2. **Short-term**: Add PennCNV as a complementary caller for constitutional (germline) CNV detection. A wrapper script could extract LRR/BAF from the pipeline VCF, generate PFB and GC model files, and run PennCNV:
   ```bash
   # Example extraction for PennCNV input
   bcftools query -f '%CHROM\t%POS\t%ID\t[%LRR]\t[%BAF]\n' \
       stage2_reclustered.bcf > sample.penncnv
   ```

3. **Long-term**: Consider EnsembleCNV to combine PennCNV and other callers for improved accuracy in large cohort studies.
