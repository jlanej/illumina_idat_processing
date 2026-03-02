# QC Comparison: Pre vs Post Reclustering

After Stage 2 completes, the pipeline automatically runs `scripts/plot_qc_comparison.py` to generate a comprehensive before/after comparison of QC metrics. The output is written to `stage2/qc/comparison/`.

## Overview

Reclustering builds a study-specific EGT cluster file from the high-quality Stage 1 samples and then re-calls all samples using that new cluster file. The QC comparison step quantifies exactly how much this improves genotype quality.

The comparison covers:

- **Per-sample call rate** — fraction of variants with a genotype call
- **Per-sample LRR SD** — standard deviation of the Log R Ratio (proxy for signal noise, critical for CNV/mosaic detection)
- **Variant missingness** — fraction of samples missing a genotype call per variant
- **Hardy-Weinberg equilibrium (HWE)** — conformance to expected allele frequency ratios
- **Minor allele frequency (MAF)** — allele frequency distribution across variants

## Output Files

| File | Description |
|------|-------------|
| `qc_comparison_report.txt` | Text summary report with matched-sample statistics and threshold counts |
| `qc_comparison_samples.png` | 2×3 panel of per-sample call rate and LRR SD comparisons |
| `qc_comparison_variants.png` | 1×3 panel of variant missingness, HWE, and MAF comparisons |

## Running the Comparison Manually

The comparison script is run automatically by `stage2_recluster.sh`, but can also be run manually:

```bash
python3 scripts/plot_qc_comparison.py \
    --stage1-sample-qc stage1/qc/stage1_sample_qc.tsv \
    --stage2-sample-qc stage2/qc/stage2_sample_qc.tsv \
    --stage1-variant-qc-dir stage1/qc/variant_qc \
    --stage2-variant-qc-dir stage2/qc/variant_qc \
    --output-dir stage2/qc/comparison
```

Variant QC directories are optional; the script produces sample-level comparisons even without them.

## Interpreting the Plots

### Sample QC Plot (`qc_comparison_samples.png`)

The plot contains two rows (call rate and LRR SD), each with three panels:

**Left — Distribution overlay**: Histograms of Stage 1 (blue) and Stage 2 (orange) values. A red dashed line marks the QC threshold (call rate ≥ 0.97, LRR SD ≤ 0.35). After reclustering, you expect:
  - Call rate: slight rightward shift (more high-quality calls)
  - LRR SD: clear leftward shift (lower noise = better signal)

**Center — Scatter plot (Stage 1 vs Stage 2)**: Each dot is one sample. Points above the diagonal (call rate) or below it (LRR SD) represent samples that improved. A tight cluster below the diagonal in the LRR SD panel indicates broad, consistent improvement.

**Right — Change distribution**: Histogram of per-sample metric change (Stage 2 − Stage 1). A positive mean (call rate) or negative mean (LRR SD) indicates improvement. The green line marks the mean change.

### Variant QC Plot (`qc_comparison_variants.png`)

Three panels show distributions before (blue) and after (orange) reclustering:

**Variant Missingness**: Lower values are better. After reclustering, the peak near zero should be taller and more variants should have near-zero missingness.

**Hardy-Weinberg Equilibrium (−log10 p-value)**: Plotted on a −log10 scale; variants well under a red dashed line at p=1e-6 are considered HWE failures. Reclustering should reduce extreme HWE deviations caused by cluster misassignment.

**MAF Distribution**: Should be broadly stable between stages. Systematic shifts could indicate population stratification or allele swaps in the cluster file.

## Example Results: 1000 Genomes (1,000 samples)

The `tests/data/comparison/` directory contains example output generated from 1,000 samples of the 1000 Genomes Project processed with the HumanOmni2.5-4v1 array.

### Text Report (`qc_comparison_report.txt`)

```
==================================================================
  QC Comparison Report: Pre vs Post Reclustering
==================================================================

Matched samples: 1000
Stage 1 only:    0
Stage 2 only:    0

------------------------------------------------------------------
  Sample Call Rate (autosomes)
------------------------------------------------------------------
  Stage 1:  mean=0.9965  median=0.9968  min=0.9811  max=0.9985  n=1000
  Stage 2:  mean=0.9971  median=0.9975  min=0.9816  max=0.9988  n=1000
  Change:  mean=+0.0005  median=+0.0005
  Improved: 938/1000 (93.8%) samples
  Passing CR>=0.95: 1000 -> 1000 (+0)
  Passing CR>=0.97: 1000 -> 1000 (+0)
  Passing CR>=0.98: 1000 -> 1000 (+0)
  Passing CR>=0.99: 993 -> 992 (-1)

------------------------------------------------------------------
  Sample LRR SD (autosomes)
------------------------------------------------------------------
  Stage 1:  mean=0.1997  median=0.1946  min=0.1540  max=0.5088  n=1000
  Stage 2:  mean=0.1461  median=0.1355  min=0.1071  max=0.5288  n=1000
  Change:  mean=-0.0536  median=-0.0573
  Improved: 967/1000 (96.7%) samples
  Passing LRR_SD<=0.2: 631 -> 935 (+304)
  Passing LRR_SD<=0.25: 952 -> 968 (+16)
  Passing LRR_SD<=0.3: 987 -> 985 (-2)
  Passing LRR_SD<=0.35: 996 -> 994 (-2)

------------------------------------------------------------------
  Variant Missingness
------------------------------------------------------------------
  Stage 1: n=2389263  mean=0.003461  median=0.001000
  Stage 2: n=2389263  mean=0.002942  median=0.001000
  Variants <=0.01: 2284343 -> 2302478 (+18135)
  Variants <=0.02: 2342084 -> 2355303 (+13219)
  Variants <=0.05: 2371433 -> 2378433 (+7000)

------------------------------------------------------------------
  HWE p-value
------------------------------------------------------------------
  Stage 1: n=2389327  mean=0.267920  median=0.148299
  Stage 2: n=2389327  mean=0.268628  median=0.149922
  Variants >=1e-06: 2204543 -> 2205649 (+1106)
  Variants >=1e-10: 2308327 -> 2308915 (+588)
  Variants >=1e-20: 2371073 -> 2371239 (+166)

------------------------------------------------------------------
  Minor Allele Frequency
------------------------------------------------------------------
  Stage 1: n=2389231  mean=nan  median=0.181225
  Stage 2: n=2389231  mean=nan  median=0.181500
  Variants >=0.01: 1889522 -> 1890989 (+1467)
  Variants >=0.05: 1387760 -> 1388465 (+705)

==================================================================
```

> **Note on `mean=nan` in the MAF section**: The `mean=nan` reflects a NaN (not-a-number) value output by plink2 for the mean allele frequency — this can occur when the distribution includes monomorphic variants (MAF = 0) or extreme values. The `median` is a more robust summary and is always shown.

### Key Takeaways

**LRR SD is the most dramatic improvement**: The median LRR SD drops from 0.1946 to 0.1355 (a 30% reduction), and samples passing the strict LRR SD ≤ 0.20 threshold jump from 63% (631/1000) to 94% (935/1000). This directly translates to improved sensitivity for mosaic chromosomal alteration detection with MoChA.

**Call rate improves broadly**: 93.8% of samples see improved call rate, with a mean gain of 0.05 percentage points. All 1,000 samples pass the standard 0.97 threshold in both stages, but the distribution tightens.

**Variant missingness decreases**: 18,135 additional variants meet the common 1% missingness cutoff after reclustering, reflecting better overall genotype calling at the variant level.

**HWE and MAF are stable**: Small improvements in HWE conformance (1,106 more variants passing p ≥ 1e-6) confirm that reclustering corrects some genotype misassignments without introducing systematic allele frequency distortion.

### Example Plots

**Sample QC:**

![Sample QC comparison](../tests/data/comparison/qc_comparison_samples.png)

**Variant QC:**

![Variant QC comparison](../tests/data/comparison/qc_comparison_variants.png)

## What to Expect in Your Data

The magnitude of improvement depends on:

- **How different your population is from the Illumina training cohort**: Arrays trained on predominantly European cohorts often show larger LRR SD improvements when applied to admixed or non-European populations.
- **Sample quality**: Low-input or degraded DNA samples may show little improvement if the signal noise is truly biological.
- **Cohort size**: Reclustering requires enough high-quality samples (~100+) to estimate cluster means and standard deviations reliably. Small cohorts (< 50 samples) may fall back to original cluster definitions for many probes.

A typical expectation for a well-powered study with ≥ 500 samples:
- LRR SD reduction: 15–40%
- Samples with improved call rate: > 80%
- New variants passing 1% missingness: > 10,000
