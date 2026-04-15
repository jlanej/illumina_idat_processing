# Minor Allele Frequency (MAF) in the Pipeline

This document maps every point in the pipeline where minor allele frequency is
used — either as a filter applied to input variants or as a column in an output
file.  It is intended as a quick reference for analysts and developers.

---

## Quick Reference

| Context | Script / file | Default threshold | Role |
|---------|--------------|-------------------|------|
| PCA variant filtering | `ancestry_pca.sh --min-maf` | ≥ 0.05 | **Input filter** (hard-exclude) |
| Cross-ancestry QC flag | `collate_variant_qc.py` `MAF_THRESHOLD` | ≥ 0.01 | **Computed flag** |
| Autosomal variant QC output | `collated_variant_qc.tsv` `all_maf` / `{ANC}_maf` | — | **Output column** |
| Cross-ancestry pass flag | `collated_variant_qc.tsv` `all_ancestries_maf_pass` | ≥ 0.01 | **Output column** |
| chrX variant QC output | `collated_variant_qc.tsv` `chrX_maf` | — | **Output column** |
| chrY variant QC output | `collated_variant_qc.tsv` `chrY_male_maf` | — | **Output column** |
| Per-ancestry sex-chr output | `collated_variant_qc.tsv` `{ANC}_chrX_maf` / `{ANC}_chrY_male_maf` | — | **Output column** |
| QC comparison plot | `plot_qc_comparison.py` | — | Visualization only |
| HTML report histograms | `generate_report.py` | — | Visualization only |

---

## MAF as an Input Filter

### `ancestry_pca.sh` — PCA variant set

**Flag**: `--min-maf FLOAT` (default **0.05**)

Before running PCA, `ancestry_pca.sh` calls `plink2 --maf` to remove variants
below the specified frequency.  The stricter threshold (0.05 vs the 0.01 used
elsewhere) is intentional: very low-frequency variants contribute little to
population structure and can destabilize PC loadings.

```bash
plink2 \
  --bcf input.bcf \
  --autosome \
  --geno 0.02 \
  --hwe 1e-6 \
  --maf 0.05 \   # <-- hard filter; variants below this are dropped
  --mind 0.02 \
  --make-bed \
  --out pca_qc
```

The threshold is configurable:
```bash
ancestry_pca.sh --vcf input.bcf --output-dir pca/ --min-maf 0.01
```

After PCA, all samples are *projected* onto the computed PCs — including those
that were excluded from the training set — so the MAF filter applies only to
the variant set used for PC computation, not to downstream analyses.

---

## MAF as Computed Pass / Fail Flags

### `collate_variant_qc.py` — cross-ancestry MAF pass flag

**Constant**: `MAF_THRESHOLD = 0.01` (hard-coded, line 41)

When ancestry-stratified variant QC results are collated, the script evaluates
whether each variant clears MAF ≥ 0.01 in **every** ancestry group that was
analyzed.  The result is written as the `all_ancestries_maf_pass` column (1 /
0) in `collated_variant_qc.tsv`.

```
all_ancestries_maf_pass = 1  iff  {ANC}_maf >= 0.01  for every ancestry ANC
```

If no ancestry-stratified data are available the column is written as `NA`.

---

## MAF as Output Columns

All MAF columns are computed from `plink2 --freq`, which writes `.afreq` files
(one row per variant, `ALT_FREQS` field).  MAF is then derived as
`min(ALT_FREQS, 1 − ALT_FREQS)` and stored as a 6-decimal-place string.

### `collated_variant_qc.tsv` — autosomal columns

Produced by `collate_variant_qc.py` (fed from `compute_variant_qc.sh` output).

| Column | Description |
|--------|-------------|
| `all_maf` | Full-cohort autosomal MAF |
| `{ANC}_maf` | Per-ancestry autosomal MAF (e.g. `EUR_maf`, `AFR_maf`) |
| `all_ancestries_maf_pass` | 1 if `{ANC}_maf ≥ 0.01` for all ancestries |

### `collated_variant_qc.tsv` — sex-chromosome columns

Produced by `collate_variant_qc.py` (fed from `compute_sex_chr_variant_qc.sh`
output).  Available only when sex-chromosome QC data exist.

| Column | Description |
|--------|-------------|
| `chrX_maf` | chrX MAF, all samples (non-PAR/XTR excluded) |
| `chrY_male_maf` | chrY MAF, males only (PAR excluded) |
| `{ANC}_chrX_maf` | chrX MAF within ancestry group |
| `{ANC}_chrY_male_maf` | chrY male MAF within ancestry group |

`chrX_female_maf` is **not** a column; MAF on chrX is reported for all samples
together because the ALT allele frequency on chrX is most informative at the
population level.  HWE (diploid assumption) is restricted to females.

### Raw `.afreq` files

`compute_variant_qc.sh` and `compute_sex_chr_variant_qc.sh` write intermediate
`plink2` `.afreq` files that feed into the steps above:

| File | Scope |
|------|-------|
| `<output-dir>/<prefix>.afreq` | Autosomal, full cohort (per stage) |
| `sex_chr_qc/chrX/chrX_all.afreq` | chrX, all samples |
| `sex_chr_qc/chrY/chrY_male.afreq` | chrY, males only |

These files are not meant for direct downstream consumption but are readable
with any tab-separated reader (`ID`, `REF`, `ALT`, `ALT_FREQS`, `OBS_CT`
columns).

---

## MAF in Diagnostics and Visualization

### `compute_variant_qc.sh` — per-run summary

After computing `.afreq`, the script prints a MAF distribution summary to
stdout (not saved to a file):

```
--- Allele Frequency Distribution ---
Monomorphic (MAF<0.1%):     3,421 (0.14%)
Rare (MAF<1%):             52,183 (2.18%)
Low frequency (1-5%):     137,482 (5.76%)
Common (MAF>=5%):       2,196,177 (91.92%)
Mean MAF: 0.1892
```

### `plot_qc_comparison.py` — before/after reclustering plot

Reads the `.afreq` files from Stage 1 and Stage 2 variant QC directories and
plots overlapping histograms of the MAF distribution.  Saved to
`stage2/qc/comparison/qc_comparison_variants.png` (third panel).

A stable MAF distribution across stages confirms that reclustering did not
introduce systematic allele-frequency distortion.  Systematic shifts can
indicate population stratification artefacts or allele swaps in the cluster
file.

### `generate_report.py` — interactive HTML report

The HTML report (`generate_report.py`) reads `collated_variant_qc.tsv` and
generates:

- **MAF histogram** (pre-binned at 0.01 intervals over [0, 0.5]) — one panel
  per ancestry group plus a full-cohort panel; interactive slider to zoom the
  frequency range.
- **MAF pass-rate bar chart** — percentage of variants with MAF ≥ 0.01, shown
  per ancestry group so that ancestry-specific differences are immediately
  visible.
- **QC thresholds table** — lists MAF ≥ 0.01 (variant QC) and MAF ≥ 0.05
  (PCA) with citations.

---

## Threshold Rationale

| Threshold | Where applied | Rationale |
|-----------|--------------|-----------|
| MAF ≥ 0.01 | `all_ancestries_maf_pass` flag; QC report | Rare variants have low statistical power in GWAS and higher genotyping error sensitivity (Marees et al. 2018) |
| MAF ≥ 0.05 | `ancestry_pca.sh` default | Ensures stable PC loadings by excluding low-frequency variants that contribute noise rather than population signal (Price et al. 2006; Marees et al. 2018) |

Both thresholds are adjustable at the command line (`--min-maf` for the PCA
step); the collation threshold is a constant in `collate_variant_qc.py` and
would require a code change to modify.

---

## References

- Anderson C.A. et al. *Data quality control in genetic case-control
  association studies.* Nat Protoc **5**, 1564–1573 (2010).
  [DOI: 10.1038/nprot.2010.116](https://doi.org/10.1038/nprot.2010.116)
- Marees A.T. et al. *A tutorial on conducting genome-wide association studies:
  Quality control and statistical analysis.* Int J Methods Psychiatr Res
  **27**, e1608 (2018).
  [DOI: 10.1002/mpr.1608](https://doi.org/10.1002/mpr.1608)
- Price A.L. et al. *Principal components analysis corrects for stratification
  in genome-wide association studies.* Nat Genet **38**, 904–909 (2006).
  [DOI: 10.1038/ng1847](https://doi.org/10.1038/ng1847)
