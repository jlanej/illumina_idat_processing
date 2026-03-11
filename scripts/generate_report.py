#!/usr/bin/env python3
"""
generate_report.py

Generate a comprehensive HTML pipeline report with publication-quality
figures, summary statistics, and auto-generated methods text.

This is the final step of the pipeline, producing a single self-contained
HTML report with:
  - Pipeline configuration and runtime summary
  - Sample QC dashboard (call rate vs LRR SD scatter, BAF SD, het rate distributions)
  - Per-sample detail table with QC failure reasons and downloadable pass/fail lists
  - Manifest realignment summary
  - Reclustering comparison (if Stage 2 was run)
  - Sex check visualization with annotation guide for common sex chromosome patterns
  - Ancestry PCA scatter plots (PC1 vs PC2, PC3 vs PC4) with k-means clustering
  - Variant QC summary
  - Peddy pedigree/sex/ancestry QC
  - Auto-generated methods paragraph for publications
  - Summary statistics table in publication-ready format

Usage:
    python3 generate_report.py --output-dir /path/to/pipeline_output

All input files are auto-detected from the standard output directory structure.
"""

import argparse
import base64
import io
import json
import math
import os
import shutil
import subprocess
import sys
import statistics
from datetime import datetime


# ---------------------------------------------------------------------------
# GWAS QC best-practice thresholds
# References:
#   - Anderson et al. (2010) Nat Protoc 5:1564-73. doi:10.1038/nprot.2010.116
#   - Marees et al. (2018) Int J Methods Psychiatr Res 27:e1608. doi:10.1002/mpr.1608
#   - Turner et al. (2011) Curr Protoc Hum Genet Ch.1:Unit1.19.
#     doi:10.1002/0471142905.hg0119s68
# ---------------------------------------------------------------------------
GWAS_THRESHOLDS = {
    'call_rate_min': 0.97,
    'lrr_sd_max': 0.35,
    'baf_sd_max': 0.15,
    'het_rate_sd_multiplier': 3,
    'variant_missingness_max': 0.02,
    'hwe_pvalue_min': 1e-6,
    'maf_min': 0.01,
    'inbreeding_f_max': 0.05,
}


GWAS_METHODS_CITATIONS = [
    {
        'citation': 'Anderson et al. 2010, Nat Protoc 5:1564-1573',
        'doi': '10.1038/nprot.2010.116',
        'cited_for': 'Foundational GWAS sample/variant QC guidance',
        'pipeline_application': (
            'Supports sample call rate and heterozygosity outlier review; '
            'supports variant call rate and HWE filtering guidance.'
        ),
    },
    {
        'citation': 'Marees et al. 2018, Int J Methods Psychiatr Res 27:e1608',
        'doi': '10.1002/mpr.1608',
        'cited_for': 'Modern practical GWAS QC workflow recommendations',
        'pipeline_application': (
            'Supports contamination/noise interpretation (BAF SD, heterozygosity), '
            'variant QC thresholds, and stricter PCA MAF filtering for stability.'
        ),
    },
    {
        'citation': 'Turner et al. 2011, Curr Protoc Hum Genet Unit 1.19',
        'doi': '10.1002/0471142905.hg0119s68',
        'cited_for': 'Array-focused genotype QC and filtering practice',
        'pipeline_application': (
            'Supports LRR/BAF quality interpretation, variant missingness thresholds, '
            'and inbreeding/relatedness flagging context.'
        ),
    },
    {
        'citation': 'Danecek et al. 2021, Gigascience 10(2):giab008',
        'doi': '10.1093/gigascience/giab008',
        'cited_for': 'bcftools framework for variant processing',
        'pipeline_application': (
            'Supports IDAT-to-GTC/VCF processing and downstream variant manipulations.'
        ),
    },
    {
        'citation': 'Li and Durbin 2009, Bioinformatics 25(14):1754-1760',
        'doi': '10.1093/bioinformatics/btp324',
        'cited_for': 'BWA alignment algorithm',
        'pipeline_application': (
            'Supports manifest probe flank realignment to the selected reference genome.'
        ),
    },
    {
        'citation': 'Chang et al. 2015, Gigascience 4:7',
        'doi': '10.1186/s13742-015-0047-8',
        'cited_for': 'PLINK analytical framework for large-scale genotype QC',
        'pipeline_application': (
            'Supports variant-level QC metrics, HWE/missingness/MAF filtering, '
            'and LD pruning in ancestry PCA preparation.'
        ),
    },
    {
        'citation': 'Abraham et al. 2017, Bioinformatics 33(17):2776-2778',
        'doi': '10.1093/bioinformatics/btx299',
        'cited_for': 'flashpca2 randomized PCA implementation',
        'pipeline_application': (
            'Supports computationally efficient ancestry PCA for large genotype matrices.'
        ),
    },
    {
        'citation': 'Pedersen and Quinlan 2017, Am J Hum Genet 100(3):406-413',
        'doi': '10.1016/j.ajhg.2017.01.017',
        'cited_for': 'peddy pedigree/sex/ancestry QC',
        'pipeline_application': (
            'Supports automated sex check, ancestry prediction, and relatedness '
            'validation from VCF genotypes via the peddy tool.'
        ),
    },
    {
        'citation': 'Genovese et al. 2024, Bioinformatics 40(2):btae038',
        'doi': '10.1093/bioinformatics/btae038',
        'cited_for': 'bcftools/liftover assembly-coordinate conversion',
        'pipeline_application': (
            'Supports coordinate conversion from non-GRCh38 builds to GRCh38 '
            'for peddy site matching via bcftools +liftover.'
        ),
    },
    {
        'citation': 'Peterson et al. 2019, Am J Hum Genet 105(5):921-935',
        'doi': '10.1016/j.ajhg.2019.09.022',
        'cited_for': 'Within-ancestry PCA and ancestry-stratified QC',
        'pipeline_application': (
            'Supports ancestry-stratified variant QC and within-ancestry PCA '
            'to resolve finer population structure and avoid Wahlund effect in '
            'HWE testing across mixed-ancestry cohorts.'
        ),
    },
]


def get_tool_versions():
    """Collect versions of key tools."""
    versions = {}
    for tool, cmd in [
        ('bcftools', ['bcftools', '--version']),
        ('plink2', ['plink2', '--version']),
        ('bwa', ['bwa']),  # bwa prints version to stderr
        ('flashpca', ['flashpca', '--version']),
        ('peddy', ['python3', '-m', 'peddy', '--version']),
    ]:
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
            output = (proc.stdout + proc.stderr).strip().split('\n')[0]
            versions[tool] = output
        except Exception:
            versions[tool] = ''
    return versions


def read_tsv(filepath):
    """Read a TSV file, return (header, list of row dicts)."""
    if not os.path.exists(filepath):
        return [], []
    rows = []
    with open(filepath) as f:
        header_line = f.readline().strip()
        if not header_line:
            return [], []
        header = header_line.split('\t')
        for line in f:
            fields = line.strip().split('\t')
            row = {}
            for i, col in enumerate(header):
                row[col] = fields[i] if i < len(fields) else 'NA'
            rows.append(row)
    return header, rows


def read_text(filepath):
    """Read a text file, return contents."""
    if not os.path.exists(filepath):
        return None
    with open(filepath) as f:
        return f.read()


def read_csv(filepath):
    """Read a CSV file, return (header, list of row dicts)."""
    import csv
    if not os.path.exists(filepath):
        return [], []
    rows = []
    with open(filepath) as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames or []
        for row in reader:
            rows.append(row)
    return header, rows


def safe_float(val, default=None):
    """Convert string to float, returning default on failure."""
    if val in ('NA', '', None):
        return default
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def fig_to_base64(fig):
    """Convert matplotlib figure to base64-encoded PNG for HTML embedding."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')


def create_sample_qc_dashboard(qc_rows, stage_label=''):
    """Create a multi-panel sample QC dashboard figure."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
    except (ImportError, AttributeError):
        return None

    # Extract values
    call_rates = [safe_float(r.get('call_rate')) for r in qc_rows]
    lrr_sds = [safe_float(r.get('lrr_sd')) for r in qc_rows]
    baf_sds = [safe_float(r.get('baf_sd')) for r in qc_rows]
    het_rates = [safe_float(r.get('het_rate')) for r in qc_rows]

    call_rates = [v for v in call_rates if v is not None]
    lrr_sds = [v for v in lrr_sds if v is not None]
    baf_sds = [v for v in baf_sds if v is not None]
    het_rates = [v for v in het_rates if v is not None]

    n_panels = 2  # always have call rate vs LRR SD and distributions
    if baf_sds:
        n_panels += 1
    if het_rates:
        n_panels += 1

    fig = plt.figure(figsize=(18, 5 * ((n_panels + 2) // 3)))
    gs = gridspec.GridSpec((n_panels + 2) // 3, 3, figure=fig, hspace=0.35, wspace=0.3)

    title_suffix = f' ({stage_label})' if stage_label else ''
    fig.suptitle(f'Sample QC Dashboard{title_suffix}',
                 fontsize=16, fontweight='bold', y=0.98)

    panel_idx = 0

    # Panel 1: Call Rate vs LRR SD scatter
    if call_rates and lrr_sds:
        ax = fig.add_subplot(gs[panel_idx // 3, panel_idx % 3])
        # Color by pass/fail
        pass_cr = []
        pass_sd = []
        fail_cr = []
        fail_sd = []
        for r in qc_rows:
            cr = safe_float(r.get('call_rate'))
            sd = safe_float(r.get('lrr_sd'))
            if cr is None or sd is None:
                continue
            if cr >= 0.97 and sd <= 0.35:
                pass_cr.append(cr)
                pass_sd.append(sd)
            else:
                fail_cr.append(cr)
                fail_sd.append(sd)

        if pass_cr:
            ax.scatter(pass_cr, pass_sd, alpha=0.4, s=12, c='#2166AC',
                       label=f'Pass (n={len(pass_cr)})', edgecolors='none')
        if fail_cr:
            ax.scatter(fail_cr, fail_sd, alpha=0.6, s=18, c='#B2182B', marker='x',
                       label=f'Fail (n={len(fail_cr)})')

        ax.axhline(y=0.35, color='red', linestyle='--', alpha=0.5, linewidth=1)
        ax.axvline(x=0.97, color='red', linestyle='--', alpha=0.5, linewidth=1)
        ax.set_xlabel('Call Rate', fontsize=11)
        ax.set_ylabel('LRR SD', fontsize=11)
        ax.set_title('Call Rate vs LRR SD', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9, loc='upper left')
        ax.grid(True, alpha=0.2)
        panel_idx += 1

    # Panel 2: Call Rate + LRR SD distributions
    if call_rates:
        ax = fig.add_subplot(gs[panel_idx // 3, panel_idx % 3])
        ax.hist(call_rates, bins=50, color='#2166AC', alpha=0.7, edgecolor='white')
        ax.axvline(x=0.97, color='red', linestyle='--', alpha=0.7, linewidth=1.5,
                   label='Threshold (0.97)')
        mean_cr = statistics.mean(call_rates)
        ax.axvline(x=mean_cr, color='green', linestyle='-', alpha=0.7, linewidth=1.5,
                   label=f'Mean ({mean_cr:.4f})')
        ax.set_xlabel('Call Rate', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title('Call Rate Distribution', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2, axis='y')
        panel_idx += 1

    # Panel 3: LRR SD distribution
    if lrr_sds:
        ax = fig.add_subplot(gs[panel_idx // 3, panel_idx % 3])
        ax.hist(lrr_sds, bins=50, color='#D6604D', alpha=0.7, edgecolor='white')
        ax.axvline(x=0.35, color='red', linestyle='--', alpha=0.7, linewidth=1.5,
                   label='Threshold (0.35)')
        mean_sd = statistics.mean(lrr_sds)
        ax.axvline(x=mean_sd, color='green', linestyle='-', alpha=0.7, linewidth=1.5,
                   label=f'Mean ({mean_sd:.4f})')
        ax.set_xlabel('LRR SD', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title('LRR SD Distribution', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2, axis='y')
        panel_idx += 1

    # Panel 4: BAF SD distribution
    if baf_sds:
        ax = fig.add_subplot(gs[panel_idx // 3, panel_idx % 3])
        ax.hist(baf_sds, bins=50, color='#4393C3', alpha=0.7, edgecolor='white')
        ax.axvline(x=0.15, color='red', linestyle='--', alpha=0.7, linewidth=1.5,
                   label='Contamination flag (0.15)')
        mean_baf = statistics.mean(baf_sds)
        ax.axvline(x=mean_baf, color='green', linestyle='-', alpha=0.7, linewidth=1.5,
                   label=f'Mean ({mean_baf:.4f})')
        ax.set_xlabel('BAF SD (heterozygous sites)', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title('BAF SD Distribution', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2, axis='y')
        panel_idx += 1

    # Panel 5: Heterozygosity rate distribution
    if het_rates:
        ax = fig.add_subplot(gs[panel_idx // 3, panel_idx % 3])
        ax.hist(het_rates, bins=50, color='#5AAE61', alpha=0.7, edgecolor='white')
        mean_hr = statistics.mean(het_rates)
        sd_hr = statistics.stdev(het_rates) if len(het_rates) > 1 else 0
        ax.axvline(x=mean_hr, color='green', linestyle='-', alpha=0.7, linewidth=1.5,
                   label=f'Mean ({mean_hr:.4f})')
        if sd_hr > 0:
            ax.axvline(x=mean_hr - 3 * sd_hr, color='red', linestyle='--', alpha=0.5,
                       linewidth=1, label=f'±3 SD ({3*sd_hr:.4f})')
            ax.axvline(x=mean_hr + 3 * sd_hr, color='red', linestyle='--', alpha=0.5,
                       linewidth=1)
        ax.set_xlabel('Heterozygosity Rate', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title('Heterozygosity Rate Distribution', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2, axis='y')
        panel_idx += 1

    return fig


def _read_pca_projections(pca_file):
    """Read PCA projections as a dict of sample_id -> list(PC values)."""
    if not os.path.exists(pca_file):
        return {}

    pca_data = {}
    with open(pca_file) as f:
        header = f.readline()  # skip header
        for line in f:
            fields = line.strip().split('\t') if '\t' in line else line.strip().split()
            if len(fields) < 4:
                continue
            # IID is second column in PLINK/flashpca output (FID, IID, PC1, PC2, ...)
            sample_id = fields[1] if len(fields) > 1 else fields[0]
            try:
                pcs = [float(x) for x in fields[2:]]
                pca_data[sample_id] = pcs
            except (ValueError, IndexError):
                pass
    return pca_data


def create_pca_plot(pca_file, qc_rows):
    """Create PCA scatter plots (PC1 vs PC2, PC3 vs PC4) colored by predicted sex."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except (ImportError, AttributeError):
        return None

    if not os.path.exists(pca_file):
        return None

    pca_data = _read_pca_projections(pca_file)

    if not pca_data:
        return None

    # Build sex map from QC rows
    sex_map = {}
    for r in qc_rows:
        sid = r.get('sample_id', '')
        sex = r.get('computed_gender', 'NA')
        sex_map[sid] = sex

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Ancestry PCA', fontsize=14, fontweight='bold')

    sex_colors = {
        '1': ('#4393C3', 'Male'), 'M': ('#4393C3', 'Male'),
        '2': ('#D6604D', 'Female'), 'F': ('#D6604D', 'Female'),
    }

    for ax_idx, (pc_x, pc_y, title) in enumerate([
        (0, 1, 'PC1 vs PC2'), (2, 3, 'PC3 vs PC4')
    ]):
        ax = axes[ax_idx]

        # Group samples by sex for legend
        groups = {}  # label -> (x_vals, y_vals, color)
        for sid, pcs in pca_data.items():
            if len(pcs) <= max(pc_x, pc_y):
                continue
            sex = sex_map.get(sid, 'NA')
            if sex in sex_colors:
                color, label = sex_colors[sex]
            else:
                color, label = '#999999', 'Unknown'
            groups.setdefault(label, ([], [], color))
            groups[label][0].append(pcs[pc_x])
            groups[label][1].append(pcs[pc_y])

        for label, (xs, ys, color) in sorted(groups.items()):
            ax.scatter(xs, ys, c=color, label=f'{label} (n={len(xs)})',
                       alpha=0.5, s=10, edgecolors='none')

        ax.set_xlabel(f'PC{pc_x + 1}', fontsize=11)
        ax.set_ylabel(f'PC{pc_y + 1}', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2)

    plt.tight_layout()
    return fig


def compute_summary_stats(qc_rows):
    """Compute summary statistics table for the report."""
    stats = {}
    n = len(qc_rows)
    stats['n_samples'] = n

    for metric, col in [
        ('call_rate', 'call_rate'), ('lrr_sd', 'lrr_sd'),
        ('lrr_mean', 'lrr_mean'), ('lrr_median', 'lrr_median'),
        ('baf_sd', 'baf_sd'), ('het_rate', 'het_rate'),
    ]:
        vals = [safe_float(r.get(col)) for r in qc_rows]
        vals = [v for v in vals if v is not None]
        if vals:
            stats[f'{metric}_mean'] = statistics.mean(vals)
            stats[f'{metric}_median'] = statistics.median(vals)
            stats[f'{metric}_sd'] = statistics.stdev(vals) if len(vals) > 1 else 0
            stats[f'{metric}_min'] = min(vals)
            stats[f'{metric}_max'] = max(vals)
            stats[f'{metric}_n'] = len(vals)

    # Count pass/fail
    n_pass = sum(1 for r in qc_rows
                 if safe_float(r.get('call_rate'), 0) >= GWAS_THRESHOLDS['call_rate_min']
                 and safe_float(r.get('lrr_sd'), 999) <= GWAS_THRESHOLDS['lrr_sd_max'])
    stats['n_pass'] = n_pass
    stats['n_fail'] = n - n_pass
    stats['pass_rate'] = (n_pass / n * 100) if n > 0 else 0

    # BAF SD contamination flags
    stats['n_baf_flag'] = sum(
        1 for r in qc_rows
        if safe_float(r.get('baf_sd'), 0) > GWAS_THRESHOLDS['baf_sd_max'])

    # Het rate outliers (±3 SD from mean)
    het_mean = stats.get('het_rate_mean')
    het_sd = stats.get('het_rate_sd')
    if het_mean is not None and het_sd is not None and het_sd > 0:
        stats['n_het_outlier'] = sum(
            1 for r in qc_rows
            if safe_float(r.get('het_rate')) is not None
            and abs(safe_float(r.get('het_rate')) - het_mean)
            > GWAS_THRESHOLDS['het_rate_sd_multiplier'] * het_sd)
    else:
        stats['n_het_outlier'] = 0

    # Sex counts
    males = sum(1 for r in qc_rows if r.get('computed_gender') in ('1', 'M'))
    females = sum(1 for r in qc_rows if r.get('computed_gender') in ('2', 'F'))
    stats['n_male'] = males
    stats['n_female'] = females

    return stats


def generate_methods_text(stats, tool_versions, genome='CHM13'):
    """Generate a methods paragraph suitable for publication."""
    n = stats.get('n_samples', 0)
    cr = stats.get('call_rate_mean', 0)
    sd = stats.get('lrr_sd_mean', 0)
    n_pass = stats.get('n_pass', 0)

    bcftools_ver = tool_versions.get('bcftools', '')
    plink2_ver = tool_versions.get('plink2', '')

    # Format tool citations cleanly
    bcf_cite = f"bcftools {bcftools_ver}" if bcftools_ver else 'bcftools'
    plink_cite = f"plink2 {plink2_ver}" if plink2_ver else 'plink2'

    text = (
        f"Genotyping intensity data (IDAT files) were processed using the "
        f"Illumina IDAT Processing Pipeline. Raw IDAT files were converted to "
        f"GTC format using the GenCall algorithm (bcftools +idat2gtc), "
        f"then to VCF format with B Allele Frequency (BAF) and Log R Ratio (LRR) "
        f"intensities (bcftools +gtc2vcf). "
        f"Probe coordinates were validated by realigning CSV manifest flank "
        f"sequences against the {genome} reference genome using BWA-MEM. "
        f"\n\n"
        f"A two-stage genotyping approach was employed. Stage 1 genotype calls "
        f"were generated using manufacturer-provided EGT cluster definitions. "
        f"Per-sample QC metrics (call rate, LRR standard deviation, BAF standard "
        f"deviation at heterozygous sites, and heterozygosity rate) were computed "
        f"on autosomal variants. Following common GWAS QC guidance "
        f"(Anderson et al., 2010; Marees et al., 2018; Turner et al., 2011), "
        f"high-quality samples (call rate ≥ 0.97 and "
        f"LRR SD ≤ 0.35; n = {n_pass} of {n}) were used to recompute "
        f"study-specific genotype cluster definitions (EGT file) in Stage 2. "
        f"All samples were then re-genotyped using the study-specific clusters, "
        f"with additional BAF/LRR median adjustment (--adjust-clusters). "
        f"\n\n"
        f"Variant-level QC included per-variant missingness, Hardy-Weinberg "
        f"equilibrium tests (mid-p adjustment), allele frequency estimation, "
        f"per-sample inbreeding coefficients (F statistic), and transition/"
        f"transversion ratio computation ({plink_cite}). "
        f"Ancestry principal components were computed using stringent variant QC "
        f"(missingness < 2%, HWE p ≥ 1e-6, MAF ≥ 5%), LD pruning "
        f"(window = 1000 kb, step = 1, r² < 0.1), and flashpca2. "
        f"This stricter PCA MAF threshold was used to improve loading stability. "
        f"\n\n"
        f"Sample ancestry was predicted using peddy (Pedersen and Quinlan, 2017), "
        f"which projects samples onto 1000 Genomes principal components for "
        f"ancestry classification. For each ancestry group with sufficient "
        f"sample size (≥ 100 by default), ancestry-stratified variant QC "
        f"(missingness, HWE, allele frequency) was performed to avoid inflated "
        f"HWE deviation due to the Wahlund effect in mixed-ancestry samples "
        f"(Anderson et al., 2010). Within-ancestry PCA was also computed to "
        f"resolve finer population structure masked in multi-ancestry PCA "
        f"(Peterson et al., 2019). Ancestry-specific PCs and variant QC "
        f"metrics are collated alongside the full-cohort results. "
        f"After processing, the mean call rate was {cr:.4f} and the mean "
        f"LRR SD was {sd:.4f}."
    )
    return text


def generate_html_report(output_dir, genome=None):
    """Generate the complete HTML report."""

    # Auto-detect files from standard directory structure
    final_qc = None
    stage1_qc = None
    for path in [
        os.path.join(output_dir, 'stage2', 'qc', 'stage2_sample_qc.tsv'),
        os.path.join(output_dir, 'stage1', 'qc', 'stage1_sample_qc.tsv'),
    ]:
        if os.path.exists(path):
            if final_qc is None:
                final_qc = path
            if 'stage1' in path:
                stage1_qc = path

    if final_qc is None:
        print("Error: No sample QC file found.", file=sys.stderr)
        return None

    # Read QC data
    _, qc_rows = read_tsv(final_qc)
    stage_label = 'Stage 2' if 'stage2' in final_qc else 'Stage 1'

    _, stage1_rows = read_tsv(stage1_qc) if stage1_qc else ([], [])

    # Tool versions
    tool_versions = get_tool_versions()

    # Summary statistics
    stats = compute_summary_stats(qc_rows)
    stage1_stats = None
    if stage1_rows:
        stage1_stats = compute_summary_stats(stage1_rows)

    # Realignment summary
    realign_text = read_text(os.path.join(output_dir, 'realigned_manifests', 'realign_summary.txt'))

    # Comparison report
    comparison_text = None
    for comp_dir in [
        os.path.join(output_dir, 'stage2', 'qc', 'comparison'),
    ]:
        comp_path = os.path.join(comp_dir, 'qc_comparison_report.txt')
        if os.path.exists(comp_path):
            comparison_text = read_text(comp_path)
            break

    # Variant QC summary
    variant_qc_text = None
    for vqc_dir in [
        os.path.join(output_dir, 'stage2', 'qc', 'variant_qc'),
        os.path.join(output_dir, 'stage1', 'qc', 'variant_qc'),
    ]:
        vqc_path = os.path.join(vqc_dir, 'variant_qc_summary.txt')
        if os.path.exists(vqc_path):
            variant_qc_text = read_text(vqc_path)
            break

    # Ancestry-stratified variant QC summaries
    ancestry_vqc_texts = {}  # ancestry -> summary text
    ancestry_qc_dir = os.path.join(output_dir, 'ancestry_stratified_qc')
    ancestry_strat_summary = read_text(
        os.path.join(ancestry_qc_dir, 'ancestry_stratified_summary.txt'))
    if os.path.isdir(ancestry_qc_dir):
        for entry in sorted(os.listdir(ancestry_qc_dir)):
            entry_path = os.path.join(ancestry_qc_dir, entry)
            if os.path.isdir(entry_path) and entry not in ('tmp',):
                vqc_path = os.path.join(entry_path, 'variant_qc',
                                        'variant_qc_summary.txt')
                if os.path.exists(vqc_path):
                    ancestry_vqc_texts[entry] = read_text(vqc_path)

    # Collated variant QC (for data-driven plots)
    collated_vqc_file = os.path.join(ancestry_qc_dir, 'collated_variant_qc.tsv')

    # Diagnostic report
    diag_text = read_text(os.path.join(output_dir, 'qc_diagnostic_report.txt'))

    # PCA file
    pca_file = os.path.join(output_dir, 'ancestry_pca', 'pca_projections.tsv')
    sex_check_table = os.path.join(output_dir, 'sex_check', 'sex_check_chrXY_lrr.tsv')
    pca_summary = read_text(os.path.join(output_dir, 'ancestry_pca', 'qc_summary.txt'))
    mock_notice = read_text(os.path.join(output_dir, 'mock_data_notice.txt'))

    # Peddy directory
    peddy_dir = os.path.join(output_dir, 'peddy')

    # Generate methods text — auto-detect genome from reference/ dir if not
    # explicitly provided, so the methods paragraph names the correct build.
    if genome is None:
        ref_dir = os.path.join(output_dir, 'reference')
        if os.path.isdir(ref_dir):
            for name in os.listdir(ref_dir):
                if 'chm13' in name.lower() and name.endswith(('.fa', '.fasta', '.fna')):
                    genome = 'CHM13'
                    break
                if 'grch38' in name.lower() or 'hg38' in name.lower() or 'GCA_000001405' in name:
                    genome = 'GRCh38'
                    break
                if 'g1k_v37' in name.lower() or 'grch37' in name.lower() or 'hg19' in name.lower():
                    genome = 'GRCh37'
                    break
        if genome is None:
            genome = 'CHM13'  # default
    methods_text = generate_methods_text(stats, tool_versions, genome)

    # ---- Create figures ----
    figures = {}

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Sample QC dashboard
        fig = create_sample_qc_dashboard(qc_rows, stage_label)
        if fig:
            figures['qc_dashboard'] = fig_to_base64(fig)
            fig.savefig(os.path.join(output_dir, 'qc_dashboard.png'),
                        dpi=150, bbox_inches='tight', facecolor='white')
            plt.close(fig)

        # PCA plot
        if os.path.exists(pca_file):
            fig = create_pca_plot(pca_file, qc_rows)
            if fig:
                figures['pca'] = fig_to_base64(fig)
                fig.savefig(os.path.join(output_dir, 'pca_scatter.png'),
                            dpi=150, bbox_inches='tight', facecolor='white')
                plt.close(fig)

        # Embed existing comparison plots if they exist
        for name, filename in [
            ('comparison_samples', 'qc_comparison_samples.png'),
            ('comparison_variants', 'qc_comparison_variants.png'),
        ]:
            for comp_dir in [
                os.path.join(output_dir, 'stage2', 'qc', 'comparison'),
            ]:
                img_path = os.path.join(comp_dir, filename)
                if os.path.exists(img_path):
                    with open(img_path, 'rb') as f:
                        figures[name] = base64.b64encode(f.read()).decode('utf-8')
                    break

        # Embed sex check plot
        sex_check_path = os.path.join(output_dir, 'sex_check', 'sex_check_chrXY_lrr.png')
        if os.path.exists(sex_check_path):
            with open(sex_check_path, 'rb') as f:
                figures['sex_check'] = base64.b64encode(f.read()).decode('utf-8')

    except (ImportError, AttributeError):
        print("Warning: matplotlib not available. Report will have no figures.",
              file=sys.stderr)

    # ---- Build HTML ----
    pca_json = _prepare_pca_json(pca_file, qc_rows) if os.path.exists(pca_file) else '[]'
    sex_check_json = _prepare_sex_check_json(sex_check_table)
    peddy_json = _prepare_peddy_json(peddy_dir)
    collated_vqc_json = _prepare_collated_vqc_json(collated_vqc_file)
    html = _build_html(stats, stage1_stats,
                       figures, realign_text, comparison_text,
                       variant_qc_text, diag_text, methods_text,
                       pca_summary, tool_versions, stage_label,
                       qc_rows=qc_rows, pca_json=pca_json,
                       sex_check_json=sex_check_json,
                       peddy_json=peddy_json,
                       ancestry_vqc_texts=ancestry_vqc_texts,
                       ancestry_strat_summary=ancestry_strat_summary,
                       collated_vqc_json=collated_vqc_json,
                       mock_notice=mock_notice or '')

    report_path = os.path.join(output_dir, 'pipeline_report.html')
    with open(report_path, 'w') as f:
        f.write(html)

    # Also write summary stats TSV for publication use
    stats_path = os.path.join(output_dir, 'summary_statistics.tsv')
    _write_summary_tsv(stats, stats_path)

    # Write methods text
    methods_path = os.path.join(output_dir, 'methods_text.txt')
    with open(methods_path, 'w') as f:
        f.write(methods_text)

    citations_path = os.path.join(output_dir, 'citations_summary.tsv')
    _write_citations_summary_tsv(citations_path)

    _consolidate_summary_outputs(output_dir)

    return report_path


def _consolidate_summary_outputs(output_dir):
    """Copy top-level report artifacts into a consolidated summary directory."""
    summary_dir = os.path.join(output_dir, 'summary')
    os.makedirs(summary_dir, exist_ok=True)

    # Top-level outputs commonly used for quick review.
    for filename in [
        'pipeline_report.html',
        'summary_statistics.tsv',
        'methods_text.txt',
        'citations_summary.tsv',
        'compiled_sample_sheet.tsv',
        'qc_diagnostic_report.txt',
        'qc_dashboard.png',
        'pca_scatter.png',
    ]:
        src = os.path.join(output_dir, filename)
        if os.path.isfile(src):
            shutil.copy2(src, os.path.join(summary_dir, filename))

    # Representative nested outputs frequently referenced in reports.
    for relpath, outname in [
        (os.path.join('ancestry_pca', 'pca_projections.tsv'), 'pca_projections.tsv'),
        (os.path.join('ancestry_pca', 'pca_eigenvalues.txt'), 'pca_eigenvalues.txt'),
        (os.path.join('ancestry_pca', 'qc_summary.txt'), 'ancestry_pca_qc_summary.txt'),
        (os.path.join('ancestry_pca', 'pca_qc_samples.txt'), 'ancestry_pca_qc_samples.txt'),
        (os.path.join('sex_check', 'sex_check_chrXY_lrr.png'), 'sex_check_chrXY_lrr.png'),
        (os.path.join('sex_check', 'sex_check_chrXY_lrr.tsv'), 'sex_check_chrXY_lrr.tsv'),
        (os.path.join('realigned_manifests', 'realign_summary.txt'), 'realign_summary.txt'),
        (os.path.join('stage2', 'qc', 'comparison', 'qc_comparison_report.txt'), 'qc_comparison_report.txt'),
        (os.path.join('stage2', 'qc', 'comparison', 'qc_comparison_samples.png'), 'qc_comparison_samples.png'),
        (os.path.join('stage2', 'qc', 'comparison', 'qc_comparison_variants.png'), 'qc_comparison_variants.png'),
        (os.path.join('stage2', 'qc', 'variant_qc', 'variant_qc_summary.txt'), 'variant_qc_summary_stage2.txt'),
        (os.path.join('stage1', 'qc', 'variant_qc', 'variant_qc_summary.txt'), 'variant_qc_summary_stage1.txt'),
        (os.path.join('peddy', 'peddy.het_check.csv'), 'peddy.het_check.csv'),
        (os.path.join('peddy', 'peddy.sex_check.csv'), 'peddy.sex_check.csv'),
        (os.path.join('peddy', 'peddy.ped_check.csv'), 'peddy.ped_check.csv'),
        (os.path.join('peddy', 'peddy.peddy.ped'), 'peddy.peddy.ped'),
        (os.path.join('peddy', 'peddy_final.ped'), 'peddy_final.ped'),
        (os.path.join('ancestry_stratified_qc', 'collated_variant_qc.tsv'),
         'collated_variant_qc.tsv'),
        (os.path.join('ancestry_stratified_qc', 'ancestry_stratified_summary.txt'),
         'ancestry_stratified_summary.txt'),
    ]:
        src = os.path.join(output_dir, relpath)
        if os.path.isfile(src):
            shutil.copy2(src, os.path.join(summary_dir, outname))


def _write_summary_tsv(stats, filepath):
    """Write summary statistics as a publication-ready TSV."""
    with open(filepath, 'w') as f:
        f.write("Metric\tValue\n")
        f.write(f"Total samples\t{stats.get('n_samples', 'NA')}\n")
        f.write(f"Samples passing QC\t{stats.get('n_pass', 'NA')}\n")
        f.write(f"Samples failing QC\t{stats.get('n_fail', 'NA')}\n")
        f.write(f"Males\t{stats.get('n_male', 'NA')}\n")
        f.write(f"Females\t{stats.get('n_female', 'NA')}\n")
        for metric, label in [
            ('call_rate', 'Call rate'),
            ('lrr_sd', 'LRR SD'),
            ('lrr_mean', 'LRR mean'),
            ('lrr_median', 'LRR median'),
            ('baf_sd', 'BAF SD (het sites)'),
            ('het_rate', 'Heterozygosity rate'),
        ]:
            mean = stats.get(f'{metric}_mean')
            median = stats.get(f'{metric}_median')
            sd = stats.get(f'{metric}_sd')
            mn = stats.get(f'{metric}_min')
            mx = stats.get(f'{metric}_max')
            if mean is not None:
                f.write(f"{label} (mean ± SD)\t{mean:.4f} ± {sd:.4f}\n")
                f.write(f"{label} (median)\t{median:.4f}\n")
                f.write(f"{label} (range)\t{mn:.4f} – {mx:.4f}\n")


def _write_citations_summary_tsv(filepath):
    """Write a compiled best-practice citations summary TSV."""
    with open(filepath, 'w') as f:
        f.write("Citation\tDOI\tCited for\tPipeline application\n")
        for row in GWAS_METHODS_CITATIONS:
            f.write(
                f"{row['citation']}\t{row['doi']}\t"
                f"{row['cited_for']}\t{row['pipeline_application']}\n"
            )


def _prepare_sample_json(qc_rows, stats):
    """Serialize sample QC data to JSON for interactive Plotly charts."""
    def _clean(v):
        if v is None:
            return None
        if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
            return None
        return round(v, 6)

    het_mean = stats.get('het_rate_mean', 0) or 0
    het_sd = stats.get('het_rate_sd', 0) or 0

    samples = []
    for r in qc_rows:
        cr = safe_float(r.get('call_rate'))
        lsd = safe_float(r.get('lrr_sd'))
        bsd = safe_float(r.get('baf_sd'))
        hr = safe_float(r.get('het_rate'))
        cr_pass = cr is not None and cr >= GWAS_THRESHOLDS['call_rate_min']
        lsd_pass = lsd is not None and lsd <= GWAS_THRESHOLDS['lrr_sd_max']
        passes = cr_pass and lsd_pass
        het_outlier = (hr is not None and het_sd > 0
                       and abs(hr - het_mean)
                       > GWAS_THRESHOLDS['het_rate_sd_multiplier'] * het_sd)
        # Build a human-readable failure reason string for the report table
        fail_reasons = []
        if cr is None:
            fail_reasons.append('CR missing')
        elif not cr_pass:
            fail_reasons.append(f'CR < {GWAS_THRESHOLDS["call_rate_min"]}')
        if lsd is None:
            fail_reasons.append('LRR SD missing')
        elif not lsd_pass:
            fail_reasons.append(f'LRR SD > {GWAS_THRESHOLDS["lrr_sd_max"]}')
        samples.append({
            'id': r.get('sample_id', ''),
            'call_rate': _clean(cr),
            'lrr_sd': _clean(lsd),
            'lrr_mean': _clean(safe_float(r.get('lrr_mean'))),
            'lrr_median': _clean(safe_float(r.get('lrr_median'))),
            'baf_sd': _clean(bsd),
            'het_rate': _clean(hr),
            'gender': r.get('computed_gender', 'NA'),
            'qc_pass': passes,
            'fail_reason': ', '.join(fail_reasons) if fail_reasons else '',
            'baf_flag': bsd is not None and bsd > GWAS_THRESHOLDS['baf_sd_max'],
            'het_outlier': het_outlier,
        })
    return json.dumps(samples)


def _prepare_pca_json(pca_file, qc_rows):
    """Serialize PCA projection data to JSON for interactive Plotly charts."""
    def _clean(v):
        if v is None:
            return None
        if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
            return None
        return round(v, 6)

    pca_data = _read_pca_projections(pca_file)
    if not pca_data:
        return '[]'

    sex_map = {
        r.get('sample_id', ''): r.get('computed_gender', 'NA')
        for r in (qc_rows or [])
    }
    points = []
    for sid, pcs in sorted(pca_data.items()):
        if len(pcs) < 2:
            continue
        points.append({
            'id': sid,
            'gender': sex_map.get(sid, 'NA'),
            'pc1': _clean(pcs[0]),
            'pc2': _clean(pcs[1]),
            'pc3': _clean(pcs[2]) if len(pcs) > 2 else None,
            'pc4': _clean(pcs[3]) if len(pcs) > 3 else None,
            # Preserve the full PC vector for interactive, user-configurable clustering in the HTML report.
            'all_pcs': [_clean(v) for v in pcs],
        })
    return json.dumps(points)


def _prepare_sex_check_json(sex_check_file):
    """Serialize sex check chrX/chrY median LRR data for interactive plotting."""
    if not os.path.exists(sex_check_file):
        return '[]'

    _, rows = read_tsv(sex_check_file)
    points = []
    for row in rows:
        x = safe_float(row.get('chrx_lrr_median'))
        y = safe_float(row.get('chry_lrr_median'))
        if x is None or y is None:
            continue
        points.append({
            'id': row.get('sample_id', ''),
            'gender': row.get('computed_gender', 'NA'),
            'chrx_lrr_median': round(x, 6),
            'chry_lrr_median': round(y, 6),
        })
    return json.dumps(points)


def _prepare_peddy_json(peddy_dir):
    """Serialize peddy het_check data to JSON for ancestry PCA overlay.

    Returns JSON array of objects with sample_id, ancestry prediction,
    peddy PCs, and het ratio.
    """
    het_file = os.path.join(peddy_dir, 'peddy.het_check.csv')
    if not os.path.exists(het_file):
        return '[]'

    _, rows = read_csv(het_file)
    points = []
    for row in rows:
        pc1 = safe_float(row.get('PC1'))
        pc2 = safe_float(row.get('PC2'))
        if pc1 is None or pc2 is None:
            continue
        points.append({
            'id': row.get('sample_id', ''),
            'ancestry': row.get('ancestry-prediction', 'UNKNOWN'),
            'ancestry_prob': round(safe_float(row.get('ancestry-prob', 0)) or 0, 4),
            'pc1': round(pc1, 6),
            'pc2': round(pc2, 6),
            'pc3': round(safe_float(row.get('PC3')) or 0, 6),
            'pc4': round(safe_float(row.get('PC4')) or 0, 6),
            'het_ratio': round(safe_float(row.get('het_ratio')) or 0, 6),
        })
    return json.dumps(points)


def _prepare_collated_vqc_json(collated_vqc_file):
    """Prepare collated variant QC data as JSON for interactive plots.

    Returns JSON with per-ancestry summary statistics (count of variants
    at various thresholds) rather than raw per-variant data, which would
    be too large for embedding in HTML.
    """
    if not os.path.exists(collated_vqc_file):
        return '{}'

    _, rows = read_tsv(collated_vqc_file)
    if not rows:
        return '{}'

    # Identify ancestry prefixes from headers
    header = list(rows[0].keys()) if rows else []
    prefixes = set()
    for col in header:
        if col.endswith('_call_rate') and col != 'all_call_rate':
            prefixes.add(col.replace('_call_rate', ''))
    groups = ['all'] + sorted(prefixes)

    result = {}
    for group in groups:
        cr_col = f'{group}_call_rate'
        hwe_col = f'{group}_hwe_p'
        maf_col = f'{group}_maf'

        cr_vals = [safe_float(r.get(cr_col)) for r in rows
                   if safe_float(r.get(cr_col)) is not None]
        hwe_vals = [safe_float(r.get(hwe_col)) for r in rows
                    if safe_float(r.get(hwe_col)) is not None]
        maf_vals = [safe_float(r.get(maf_col)) for r in rows
                    if safe_float(r.get(maf_col)) is not None]

        n_variants = len(cr_vals) or len(hwe_vals) or len(maf_vals)
        if n_variants == 0:
            continue

        # Missingness thresholds
        n_miss_5 = sum(1 for v in cr_vals if v < 0.95)
        n_miss_2 = sum(1 for v in cr_vals if v < 0.98)
        n_miss_1 = sum(1 for v in cr_vals if v < 0.99)

        # HWE thresholds
        n_hwe_6 = sum(1 for v in hwe_vals if v < 1e-6)
        n_hwe_10 = sum(1 for v in hwe_vals if v < 1e-10)
        n_hwe_20 = sum(1 for v in hwe_vals if v < 1e-20)

        # MAF distribution
        n_mono = sum(1 for v in maf_vals if v < 0.001)
        n_rare = sum(1 for v in maf_vals if v < 0.01)
        n_low = sum(1 for v in maf_vals if 0.01 <= v < 0.05)
        n_common = sum(1 for v in maf_vals if v >= 0.05)

        result[group] = {
            'n_variants': n_variants,
            'n_cr': len(cr_vals),
            'n_hwe': len(hwe_vals),
            'n_maf': len(maf_vals),
            'mean_cr': round(sum(cr_vals) / len(cr_vals), 6) if cr_vals else None,
            'n_miss_gt5pct': n_miss_5,
            'n_miss_gt2pct': n_miss_2,
            'n_miss_gt1pct': n_miss_1,
            'n_hwe_lt1e6': n_hwe_6,
            'n_hwe_lt1e10': n_hwe_10,
            'n_hwe_lt1e20': n_hwe_20,
            'n_mono': n_mono,
            'n_rare': n_rare,
            'n_low_freq': n_low,
            'n_common': n_common,
            'mean_maf': round(sum(maf_vals) / len(maf_vals), 6) if maf_vals else None,
            'cr_values': [round(v, 6) for v in cr_vals[:5000]],
            'hwe_values': [round(v, 10) for v in hwe_vals[:5000]],
            'maf_values': [round(v, 6) for v in maf_vals[:5000]],
        }

    return json.dumps(result)


def _get_report_css():
    """Return CSS for the HTML report."""
    return """
:root {
    --primary: #2563eb;
    --primary-dark: #1e40af;
    --primary-light: #dbeafe;
    --success: #059669;
    --success-light: #d1fae5;
    --warning: #d97706;
    --warning-light: #fef3c7;
    --danger: #dc2626;
    --danger-light: #fecaca;
    --bg: #f8fafc;
    --card-bg: #ffffff;
    --text: #1e293b;
    --text-muted: #64748b;
    --border: #e2e8f0;
}
* { box-sizing: border-box; }
body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto,
                 'Helvetica Neue', Arial, sans-serif;
    background: var(--bg); color: var(--text);
    line-height: 1.6; margin: 0; padding: 0;
}
.header {
    background: linear-gradient(135deg, var(--primary-dark), var(--primary));
    color: white; padding: 2rem 1.5rem; text-align: center;
}
.header h1 { font-size: 1.7rem; font-weight: 700; margin: 0 0 0.3rem; }
.header .subtitle { opacity: 0.85; font-size: 0.85rem; margin: 0; }

.mock-banner {
    background: #fffbeb; color: #92400e; border: 1px solid #fde68a;
    border-left: 4px solid #d97706; border-radius: 10px;
    padding: 0.85rem 1rem; margin-bottom: 1rem; font-size: 0.86rem;
}

.nav-bar {
    background: var(--card-bg); border-bottom: 1px solid var(--border);
    padding: 0.6rem 1rem; position: sticky; top: 0; z-index: 100;
    display: flex; gap: 0.4rem; flex-wrap: wrap; justify-content: center;
    box-shadow: 0 1px 3px rgba(0,0,0,0.08);
}
.nav-bar a {
    color: var(--primary); text-decoration: none; font-size: 0.8rem;
    font-weight: 500; padding: 0.2rem 0.5rem; border-radius: 4px;
    transition: background 0.2s;
}
.nav-bar a:hover { background: var(--primary-light); }

.container { max-width: 1400px; margin: 0 auto; padding: 1.5rem; }

section { margin-bottom: 2.5rem; }
section h2 {
    font-size: 1.3rem; font-weight: 700; color: var(--text);
    border-bottom: 2px solid var(--primary); padding-bottom: 0.4rem;
    margin: 0 0 1rem;
}
.card {
    background: var(--card-bg); border: 1px solid var(--border);
    border-radius: 10px; padding: 1.25rem; margin-bottom: 1rem;
    box-shadow: 0 1px 3px rgba(0,0,0,0.04);
}
.metric-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
    gap: 0.75rem; margin-bottom: 1.25rem;
}
.metric-card {
    background: var(--card-bg); border: 1px solid var(--border);
    border-radius: 10px; padding: 1rem; text-align: center;
    box-shadow: 0 1px 3px rgba(0,0,0,0.04);
    transition: transform 0.15s, box-shadow 0.15s;
}
.metric-card:hover {
    transform: translateY(-2px); box-shadow: 0 4px 12px rgba(0,0,0,0.08);
}
.metric-card .value { font-size: 1.6rem; font-weight: 700; color: var(--text); }
.metric-card .label {
    font-size: 0.7rem; color: var(--text-muted); margin-top: 0.2rem;
    text-transform: uppercase; letter-spacing: 0.04em;
}
.metric-card.pass { border-left: 4px solid var(--success); }
.metric-card.fail { border-left: 4px solid var(--danger); }
.metric-card.warn { border-left: 4px solid var(--warning); }

.pass-bar-track {
    background: #f1f5f9; border-radius: 999px; height: 10px;
    overflow: hidden; margin: 0.4rem 0;
}
.pass-bar-fill {
    height: 100%; border-radius: 999px;
    background: linear-gradient(90deg, var(--success), #34d399);
    transition: width 0.6s ease;
}

/* Tables */
.stats-table { border-collapse: collapse; width: 100%; font-size: 0.85rem; }
.stats-table th, .stats-table td {
    border: 1px solid var(--border); padding: 0.5rem 0.8rem; text-align: right;
}
.stats-table th {
    background: var(--primary); color: white; font-weight: 600; font-size: 0.8rem;
}
.stats-table td:first-child { text-align: left; font-weight: 500; }
.stats-table tr:nth-child(even) { background: #f8fafc; }
.stats-table tr:hover { background: #eff6ff; }

.threshold-table { border-collapse: collapse; width: 100%; font-size: 0.82rem; }
.threshold-table th, .threshold-table td {
    border: 1px solid var(--border); padding: 0.5rem 0.8rem;
}
.threshold-table th {
    background: #f1f5f9; font-weight: 600; text-align: left;
}
.threshold-table td:nth-child(2) { font-family: monospace; font-weight: 600; }
.badge {
    display: inline-block; padding: 0.1rem 0.5rem; border-radius: 999px;
    font-size: 0.7rem; font-weight: 600;
}
.badge-pass { background: var(--success-light); color: #065f46; }
.badge-warn { background: var(--warning-light); color: #92400e; }
.badge-info { background: var(--primary-light); color: #1e40af; }

/* Sample table */
.sample-controls { display: flex; gap: 0.75rem; align-items: center;
                   flex-wrap: wrap; margin-bottom: 0.75rem; }
.search-input {
    flex: 1; min-width: 200px; max-width: 320px;
    padding: 0.45rem 0.9rem; border: 1px solid var(--border);
    border-radius: 8px; font-size: 0.82rem; outline: none;
    transition: border-color 0.2s;
}
.search-input:focus { border-color: var(--primary); }
.sample-count { font-size: 0.8rem; color: var(--text-muted); }

.sample-table-wrap {
    max-height: 480px; overflow-y: auto;
    border: 1px solid var(--border); border-radius: 8px;
}
.sample-table { border-collapse: collapse; width: 100%; font-size: 0.78rem; }
.sample-table th {
    position: sticky; top: 0; background: var(--primary); color: white;
    padding: 0.45rem 0.65rem; cursor: pointer; user-select: none;
    font-size: 0.75rem; white-space: nowrap; z-index: 1;
}
.sample-table th:hover { background: var(--primary-dark); }
.sample-table th::after { content: ' ⇅'; opacity: 0.5; font-size: 0.65rem; }
.sample-table th.sort-asc::after { content: ' ↑'; opacity: 1; }
.sample-table th.sort-desc::after { content: ' ↓'; opacity: 1; }
.sample-table td {
    padding: 0.35rem 0.65rem; border-bottom: 1px solid var(--border);
    white-space: nowrap;
}
.sample-table tr:hover { background: #eff6ff; }
.cell-pass { color: var(--success); }
.cell-fail { color: var(--danger); font-weight: 600; }
.cell-warn { color: var(--warning); font-weight: 500; }
.status-badge {
    display: inline-block; padding: 0.1rem 0.45rem; border-radius: 999px;
    font-size: 0.68rem; font-weight: 700;
}
.status-pass { background: var(--success-light); color: #065f46; }
.status-fail { background: var(--danger-light); color: #991b1b; }
.status-flag { background: #fef3c7; color: #92400e; }

/* Plots */
.plot-grid {
    display: grid; grid-template-columns: 1fr 1fr; gap: 1rem;
}
@media (max-width: 900px) { .plot-grid { grid-template-columns: 1fr; } }
.plot-box {
    background: var(--card-bg); border: 1px solid var(--border);
    border-radius: 10px; padding: 0.5rem; min-height: 340px;
}

/* Pre-formatted and methods */
.report-pre {
    background: #1e293b; color: #e2e8f0; padding: 1rem; border-radius: 8px;
    overflow-x: auto; font-size: 0.78rem; line-height: 1.5;
    font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
}
.methods-text {
    background: var(--card-bg); border-left: 4px solid var(--primary);
    padding: 1.1rem 1.4rem; margin: 0.75rem 0; font-style: italic;
    border-radius: 0 8px 8px 0; line-height: 1.7;
}
.figure-block img {
    max-width: 100%; border: 1px solid var(--border); border-radius: 8px;
    margin: 0.5rem 0;
}
.comparison-change-good { color: var(--success); font-weight: 600; }
.comparison-change-bad { color: var(--danger); font-weight: 600; }
.footer {
    text-align: center; padding: 1.5rem; color: var(--text-muted);
    font-size: 0.75rem; border-top: 1px solid var(--border); margin-top: 2rem;
}
.tools-list { list-style: none; padding: 0; }
.tools-list li {
    padding: 0.35rem 0; border-bottom: 1px solid var(--border); font-size: 0.85rem;
}
.tools-list code {
    background: #f1f5f9; padding: 0.15rem 0.4rem; border-radius: 4px;
    font-size: 0.82rem;
}
.plot-controls-note {
    font-size: 0.78rem; color: var(--text-muted); margin-top: 0.5rem;
}
.plot-controls {
    display: flex; flex-wrap: wrap; gap: 0.6rem; align-items: center;
    margin-bottom: 0.75rem;
}
.plot-controls label {
    font-size: 0.78rem; color: var(--text-muted);
}
.plot-slider {
    width: 180px;
}
.plot-input {
    padding: 0.3rem 0.45rem; border: 1px solid var(--border); border-radius: 6px;
    font-size: 0.78rem; background: var(--card-bg);
}
.plot-btn {
    border: 1px solid var(--border); border-radius: 6px; background: #f8fafc;
    padding: 0.32rem 0.65rem; font-size: 0.76rem; cursor: pointer;
    transition: background 0.15s, color 0.15s, border-color 0.15s;
    white-space: nowrap;
}
.plot-btn:hover { border-color: var(--primary); color: var(--primary); background: var(--primary-light); }
.plot-btn.active { background: var(--primary); color: white; border-color: var(--primary); }
.plot-btn.active:hover { background: var(--primary-dark); border-color: var(--primary-dark); }
.plot-btn:disabled { opacity: 0.45; cursor: not-allowed; pointer-events: none; }
.btn-group {
    display: flex; border-radius: 6px; overflow: hidden;
    border: 1px solid var(--border);
}
.btn-group .plot-btn {
    border-radius: 0; border: none; border-right: 1px solid var(--border);
}
.btn-group .plot-btn:last-child { border-right: none; }
.btn-group .plot-btn.active { border-color: transparent; }
.control-sep {
    width: 1px; background: var(--border); height: 22px;
    align-self: center; flex-shrink: 0; margin: 0 0.15rem;
}
.plot-status { font-size: 0.75rem; color: var(--text-muted); }
.vqc-tabs {
    display: flex; gap: 0.25rem; flex-wrap: wrap; margin-bottom: 0.75rem;
}
.vqc-tab {
    border: 1px solid var(--border); border-radius: 6px 6px 0 0;
    background: #f1f5f9; padding: 0.4rem 0.9rem; font-size: 0.82rem;
    cursor: pointer; font-weight: 500; color: var(--text-muted);
    transition: background 0.15s, color 0.15s;
}
.vqc-tab:hover { background: var(--primary-light); color: var(--primary); }
.vqc-tab.active {
    background: var(--primary); color: white; border-color: var(--primary);
    font-weight: 600;
}
.vqc-panel { margin-bottom: 0.5rem; }
"""


def _get_report_js():
    """Return JavaScript for interactive Plotly charts, sample table, and sorting."""
    return r"""
document.addEventListener('DOMContentLoaded', function() {
    var el = document.getElementById('qc-data');
    if (!el) return;
    var rawData = JSON.parse(el.textContent);
    if (!rawData || rawData.length === 0) return;

    var C = {pass:'#059669', fail:'#dc2626', warn:'#d97706',
             blue:'#3b82f6', thresh:'#dc2626', mean:'#059669'};
    var baseLayout = {
        font: {family: '-apple-system, BlinkMacSystemFont, Segoe UI, Roboto, sans-serif', size: 12},
        paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
        margin: {t: 45, r: 25, b: 55, l: 60}, hovermode: 'closest'
    };
    var cfg = {responsive: true, displaylogo: false,
               modeBarButtonsToRemove: ['lasso2d','select2d']};

    /* ---- Scatter: Call Rate vs LRR SD ---- */
    var pX=[], pY=[], pI=[], fX=[], fY=[], fI=[];
    rawData.forEach(function(s) {
        if (s.call_rate === null || s.lrr_sd === null) return;
        if (s.qc_pass) { pX.push(s.call_rate); pY.push(s.lrr_sd); pI.push(s.id); }
        else { fX.push(s.call_rate); fY.push(s.lrr_sd); fI.push(s.id); }
    });
    var scTraces = [];
    if (pX.length) scTraces.push({x:pX,y:pY,text:pI,mode:'markers',type:'scatter',
        name:'Pass ('+pX.length+')',marker:{color:C.pass,size:6,opacity:0.6},
        hovertemplate:'<b>%{text}</b><br>Call Rate: %{x:.4f}<br>LRR SD: %{y:.4f}<extra>Pass</extra>'});
    if (fX.length) scTraces.push({x:fX,y:fY,text:fI,mode:'markers',type:'scatter',
        name:'Fail ('+fX.length+')',marker:{color:C.fail,size:8,symbol:'x',opacity:0.8},
        hovertemplate:'<b>%{text}</b><br>Call Rate: %{x:.4f}<br>LRR SD: %{y:.4f}<extra>Fail</extra>'});
    var scLay = Object.assign({}, baseLayout, {
        title:{text:'Call Rate vs LRR SD',font:{size:14}},
        xaxis:{title:'Call Rate',gridcolor:'#f1f5f9',zeroline:false},
        yaxis:{title:'LRR SD',gridcolor:'#f1f5f9',zeroline:false},
        shapes:[
            {type:'line',x0:0.97,x1:0.97,y0:0,y1:1,yref:'paper',line:{color:C.thresh,width:1.5,dash:'dash'}},
            {type:'line',x0:0,x1:1,y0:0.35,y1:0.35,xref:'paper',line:{color:C.thresh,width:1.5,dash:'dash'}}
        ],
        annotations:[
            {x:0.97,y:1.04,yref:'paper',text:'CR \u2265 0.97',showarrow:false,font:{size:10,color:C.thresh}},
            {x:1.04,xref:'paper',y:0.35,text:'LRR SD \u2264 0.35',showarrow:false,font:{size:10,color:C.thresh}}
        ]
    });
    var d1 = document.getElementById('plot-scatter');
    if (d1) Plotly.newPlot(d1, scTraces, scLay, cfg);

    /* ---- Histogram helper ---- */
    function hist(divId, vals, title, xLbl, thresh, thLbl, dir) {
        var d = document.getElementById(divId);
        if (!d || !vals.length) return;
        var mean = vals.reduce(function(a,b){return a+b;},0)/vals.length;
        var shapes = [
            {type:'line',x0:thresh,x1:thresh,y0:0,y1:1,yref:'paper',
             line:{color:C.thresh,width:2,dash:'dash'}},
            {type:'line',x0:mean,x1:mean,y0:0,y1:1,yref:'paper',
             line:{color:C.mean,width:2}}
        ];
        var ann = [
            {x:thresh,y:1.04,yref:'paper',text:thLbl,showarrow:false,font:{size:10,color:C.thresh}},
            {x:mean,y:1.09,yref:'paper',text:'Mean: '+mean.toFixed(4),showarrow:false,font:{size:10,color:C.mean}}
        ];
        Plotly.newPlot(d, [{x:vals,type:'histogram',
            marker:{color:C.blue,opacity:0.7,line:{color:'white',width:0.5}},
            hovertemplate:xLbl+': %{x}<br>Count: %{y}<extra></extra>'}],
            Object.assign({}, baseLayout, {
                title:{text:title,font:{size:14}},
                xaxis:{title:xLbl,gridcolor:'#f1f5f9'},
                yaxis:{title:'Count',gridcolor:'#f1f5f9'},
                shapes:shapes, annotations:ann
            }), cfg);
    }

    var cr=[], lr=[], bf=[], hr=[];
    rawData.forEach(function(s){
        if(s.call_rate!==null)cr.push(s.call_rate);
        if(s.lrr_sd!==null)lr.push(s.lrr_sd);
        if(s.baf_sd!==null)bf.push(s.baf_sd);
        if(s.het_rate!==null)hr.push(s.het_rate);
    });
    hist('plot-cr',cr,'Call Rate Distribution','Call Rate',0.97,'\u2265 0.97','min');
    hist('plot-lrr',lr,'LRR SD Distribution','LRR SD',0.35,'\u2264 0.35','max');

    if (bf.length) {
        hist('plot-baf',bf,'BAF SD Distribution (Het Sites)','BAF SD',0.15,'Flag: 0.15','max');
        var bc=document.getElementById('baf-container'); if(bc) bc.style.display='';
    }

    /* ---- Het Rate (special: ±3 SD lines) ---- */
    if (hr.length) {
        var hm=hr.reduce(function(a,b){return a+b;},0)/hr.length;
        var hv=hr.reduce(function(a,b){return a+(b-hm)*(b-hm);},0)/hr.length;
        var hs=Math.sqrt(hv);
        var hShp=[{type:'line',x0:hm,x1:hm,y0:0,y1:1,yref:'paper',line:{color:C.mean,width:2}}];
        var hAnn=[{x:hm,y:1.09,yref:'paper',text:'Mean: '+hm.toFixed(4),showarrow:false,font:{size:10,color:C.mean}}];
        if(hs>0){
            hShp.push({type:'line',x0:hm-3*hs,x1:hm-3*hs,y0:0,y1:1,yref:'paper',line:{color:C.thresh,width:1.5,dash:'dash'}});
            hShp.push({type:'line',x0:hm+3*hs,x1:hm+3*hs,y0:0,y1:1,yref:'paper',line:{color:C.thresh,width:1.5,dash:'dash'}});
            hAnn.push({x:hm-3*hs,y:1.04,yref:'paper',text:'-3 SD',showarrow:false,font:{size:10,color:C.thresh}});
            hAnn.push({x:hm+3*hs,y:1.04,yref:'paper',text:'+3 SD',showarrow:false,font:{size:10,color:C.thresh}});
        }
        var hd=document.getElementById('plot-het');
        if(hd) Plotly.newPlot(hd,[{x:hr,type:'histogram',
            marker:{color:C.blue,opacity:0.7,line:{color:'white',width:0.5}},
            hovertemplate:'Het Rate: %{x}<br>Count: %{y}<extra></extra>'}],
            Object.assign({},baseLayout,{
                title:{text:'Heterozygosity Rate Distribution',font:{size:14}},
                xaxis:{title:'Heterozygosity Rate',gridcolor:'#f1f5f9'},
                yaxis:{title:'Count',gridcolor:'#f1f5f9'},
                shapes:hShp, annotations:hAnn
            }),cfg);
        var hc=document.getElementById('het-container'); if(hc) hc.style.display='';
    }

    function normSex(g) {
        var sx = (g || 'NA').toString();
        if (sx === 'M' || sx === '1') return 'M';
        if (sx === 'F' || sx === '2') return 'F';
        return 'NA';
    }

    function renderDensityToggleButtons() {
        return []; // Scatter/Density controls are handled by external HTML buttons
    }

    function kmeans(points, k, maxIter) {
        if (!points.length || k < 1) return null;
        var dim = points[0].length;
        var n = points.length;
        k = Math.max(1, Math.min(k, n));

        var centroids = [];
        // k-means++ seeding: choose one random observed point, then sample subsequent
        // centroids with probability proportional to squared distance from nearest centroid.
        var first = Math.floor(Math.random() * n);
        centroids.push(points[first].slice());
        while (centroids.length < k) {
            var distSq = new Array(n).fill(0);
            var total = 0;
            for (var i = 0; i < n; i++) {
                var minD = Infinity;
                for (var c = 0; c < centroids.length; c++) {
                    var d = 0;
                    for (var j = 0; j < dim; j++) {
                        var delta = points[i][j] - centroids[c][j];
                        d += delta * delta;
                    }
                    if (d < minD) minD = d;
                }
                distSq[i] = minD;
                total += minD;
            }
            if (total <= 0) {
                centroids.push(points[Math.floor(Math.random() * n)].slice());
                continue;
            }
            var target = Math.random() * total;
            var acc = 0;
            var chosen = n - 1;
            for (var idx = 0; idx < n; idx++) {
                acc += distSq[idx];
                if (acc >= target) {
                    chosen = idx;
                    break;
                }
            }
            centroids.push(points[chosen].slice());
        }

        var labels = new Array(n).fill(0);
        for (var iter = 0; iter < maxIter; iter++) {
            var changed = false;
            for (var p = 0; p < n; p++) {
                var best = 0;
                var bestDist = Infinity;
                for (var cc = 0; cc < k; cc++) {
                    var dist = 0;
                    for (var dd = 0; dd < dim; dd++) {
                        var diff = points[p][dd] - centroids[cc][dd];
                        dist += diff * diff;
                    }
                    if (dist < bestDist) {
                        bestDist = dist;
                        best = cc;
                    }
                }
                if (labels[p] !== best) {
                    labels[p] = best;
                    changed = true;
                }
            }

            var sums = Array.from({length: k}, function() { return new Array(dim).fill(0); });
            var counts = new Array(k).fill(0);
            for (var pp = 0; pp < n; pp++) {
                var lab = labels[pp];
                counts[lab] += 1;
                for (var d2 = 0; d2 < dim; d2++) {
                    sums[lab][d2] += points[pp][d2];
                }
            }
            for (var c2 = 0; c2 < k; c2++) {
                if (!counts[c2]) {
                    centroids[c2] = points[Math.floor(Math.random() * n)].slice();
                    continue;
                }
                for (var d3 = 0; d3 < dim; d3++) {
                    centroids[c2][d3] = sums[c2][d3] / counts[c2];
                }
            }
            if (!changed) break;
        }
        return {labels: labels, k: k};
    }

    /* ---- Interactive PCA ---- */
    var pcaEl = document.getElementById('pca-data');
    if (pcaEl) {
        var pcaData = JSON.parse(pcaEl.textContent || '[]');
        if (pcaData && pcaData.length) {
            // Color palette for cluster labels in interactive ancestry views.
            var CLUSTER_COLORS = ['#2563eb', '#dc2626', '#059669', '#d97706', '#7c3aed', '#0f766e', '#be123c', '#4d7c0f', '#334155', '#c2410c'];
            var ANCESTRY_COLORS = {
                'EUR': '#2563eb', 'AFR': '#dc2626', 'EAS': '#059669',
                'AMR': '#d97706', 'SAS': '#7c3aed', 'UNKNOWN': '#6b7280'
            };
            var KMEANS_MAX_ITER = 60;
            /* Build a map from sample_id -> peddy ancestry prediction for overlay recoloring.
               This allows the overlay button to recolor the pipeline's own samples by their
               peddy-predicted ancestry without mixing incompatible PC coordinate spaces. */
            var peddyAncestryMap = {};
            (function() {
                var el = document.getElementById('peddy-data');
                if (!el) return;
                try {
                    JSON.parse(el.textContent || '[]').forEach(function(p) {
                        if (p.id && p.ancestry) peddyAncestryMap[p.id] = p.ancestry;
                    });
                } catch(e) {}
            })();
            var pcaState = {
                mode: 'sex',
                assignments: null,
                k: null,
                pcsUsed: null,
                binSize: 35,
                viewMode: 'scatter',
                peddyOverlayActive: false,
                plots2d: [],
                hasPc3: pcaData.some(function(p){ return p.pc3 !== null && p.pc3 !== undefined; }),
                hasPc4: pcaData.some(function(p){ return p.pc4 !== null && p.pc4 !== undefined; })
            };

            function pointLabel(point, idx) {
                if (pcaState.peddyOverlayActive && peddyAncestryMap[point.id]) {
                    return 'peddy:' + peddyAncestryMap[point.id];
                }
                if (pcaState.mode === 'cluster' && pcaState.assignments && pcaState.assignments[idx] !== undefined) {
                    return 'Cluster ' + (pcaState.assignments[idx] + 1);
                }
                var sx = normSex(point.gender);
                if (sx === 'M') return 'Male';
                if (sx === 'F') return 'Female';
                return 'Unknown';
            }

            function pointColor(point, idx) {
                if (pcaState.peddyOverlayActive && peddyAncestryMap[point.id]) {
                    return ANCESTRY_COLORS[peddyAncestryMap[point.id]] || '#6b7280';
                }
                if (pcaState.mode === 'cluster' && pcaState.assignments && pcaState.assignments[idx] !== undefined) {
                    return CLUSTER_COLORS[pcaState.assignments[idx] % CLUSTER_COLORS.length];
                }
                var sx = normSex(point.gender);
                if (sx === 'M') return '#4393C3';
                if (sx === 'F') return '#D6604D';
                return '#6b7280';
            }

            function pcaPlot(divId, xKey, yKey, title) {
                var d = document.getElementById(divId);
                if (!d) return null;
                var grouped = {};
                var allX = [];
                var allY = [];
                pcaData.forEach(function(p, idx) {
                    if (p[xKey] === null || p[yKey] === null || p[xKey] === undefined || p[yKey] === undefined) return;
                    var label = pointLabel(p, idx);
                    if (!grouped[label]) grouped[label] = {x:[], y:[], text:[], color:pointColor(p, idx)};
                    grouped[label].x.push(p[xKey]);
                    grouped[label].y.push(p[yKey]);
                    grouped[label].text.push(p.id);
                    allX.push(p[xKey]);
                    allY.push(p[yKey]);
                });
                var traces = Object.keys(grouped).sort().map(function(label) {
                    var g = grouped[label];
                    return {
                        x:g.x, y:g.y, text:g.text, mode:'markers', type:'scatter',
                        name:label+' (n='+g.x.length+')',
                        marker:{color:g.color, size:5, opacity:0.65},
                        hovertemplate:'<b>%{text}</b><br>'+xKey.toUpperCase()+': %{x:.4f}<br>'+yKey.toUpperCase()+': %{y:.4f}<extra>'+label+'</extra>'
                    };
                });
                if (!traces.length) return null;
                var densityIdx = traces.length;
                traces.push({
                    x:allX, y:allY, type:'histogram2d', name:'Density',
                    nbinsx:pcaState.binSize, nbinsy:pcaState.binSize,
                    colorscale:[
                        [0.0, 'rgba(37,99,235,0)'],
                        [0.05, 'rgba(37,99,235,0.2)'],
                        [0.4, 'rgba(37,99,235,0.5)'],
                        [1.0, 'rgba(30,64,175,0.85)']
                    ],
                    zmin:0, showscale:true, colorbar:{title:'Count'},
                    hovertemplate:xKey.toUpperCase()+': %{x:.4f}<br>'+yKey.toUpperCase()+': %{y:.4f}<br>Count: %{z}<extra>Density</extra>',
                    visible:false
                });
                Plotly.newPlot(d, traces, Object.assign({}, baseLayout, {
                    title:{text:title,font:{size:14}},
                    xaxis:{title:xKey.toUpperCase(),gridcolor:'#f1f5f9',zeroline:false},
                    yaxis:{title:yKey.toUpperCase(),gridcolor:'#f1f5f9',zeroline:false},
                    updatemenus: renderDensityToggleButtons(traces.length, densityIdx)
                }), cfg);
                return {div:d, densityIndex:densityIdx, traceCount:traces.length, n:allX.length};
            }

            function pcaPlot3d() {
                var d3 = document.getElementById('plot-pca3d');
                if (!d3 || !pcaState.hasPc3) return;
                var grouped = {};
                pcaData.forEach(function(p, idx) {
                    if (p.pc1 === null || p.pc2 === null || p.pc3 === null) return;
                    var label = pointLabel(p, idx);
                    if (!grouped[label]) grouped[label] = {x:[], y:[], z:[], text:[], color:pointColor(p, idx)};
                    grouped[label].x.push(p.pc1);
                    grouped[label].y.push(p.pc2);
                    grouped[label].z.push(p.pc3);
                    grouped[label].text.push(p.id);
                });
                var traces = Object.keys(grouped).sort().map(function(label) {
                    var g = grouped[label];
                    return {
                        type:'scatter3d', mode:'markers', name:label+' (n='+g.x.length+')',
                        x:g.x, y:g.y, z:g.z, text:g.text,
                        marker:{size:4, opacity:0.75, color:g.color},
                        hovertemplate:'<b>%{text}</b><br>PC1: %{x:.4f}<br>PC2: %{y:.4f}<br>PC3: %{z:.4f}<extra>'+label+'</extra>'
                    };
                });
                Plotly.newPlot(d3, traces, Object.assign({}, baseLayout, {
                    title:{text:'PC1 vs PC2 vs PC3 (3D)',font:{size:14}},
                    scene: {
                        xaxis:{title:'PC1',gridcolor:'#f1f5f9'},
                        yaxis:{title:'PC2',gridcolor:'#f1f5f9'},
                        zaxis:{title:'PC3',gridcolor:'#f1f5f9'}
                    },
                    margin:{t:45, r:5, b:10, l:5}
                }), cfg);
            }

            function refreshPcaPlots() {
                pcaState.plots2d = [];
                var p12 = pcaPlot('plot-pca12', 'pc1', 'pc2', 'PC1 vs PC2');
                if (p12) pcaState.plots2d.push(p12);
                var c34 = document.getElementById('pca34-container');
                if (c34) c34.style.display = '';
                if (pcaState.hasPc3 && pcaState.hasPc4) {
                    var p34 = pcaPlot('plot-pca34', 'pc3', 'pc4', 'PC3 vs PC4');
                    if (p34) pcaState.plots2d.push(p34);
                } else if (c34) {
                    c34.style.display = 'none';
                }
                pcaPlot3d();
                applyPcaViewMode(pcaState.viewMode);
            }

            function applyPcaViewMode(mode) {
                pcaState.viewMode = mode;
                var scatterBtn = document.getElementById('pca-scatter-btn');
                var densityBtn = document.getElementById('pca-density-btn');
                if (scatterBtn) scatterBtn.classList.toggle('active', mode === 'scatter');
                if (densityBtn) densityBtn.classList.toggle('active', mode === 'density');
                pcaState.plots2d.forEach(function(meta) {
                    var visible = Array.from({length: meta.traceCount}, function(_, i) {
                        return mode === 'density' ? i === meta.densityIndex : i !== meta.densityIndex;
                    });
                    Plotly.restyle(meta.div, {visible: visible});
                });
            }

            function applyPcaDensityBins(value) {
                pcaState.binSize = value;
                pcaState.plots2d.forEach(function(meta) {
                    Plotly.restyle(meta.div, {'nbinsx': value, 'nbinsy': value}, [meta.densityIndex]);
                });
            }

            function exportClusters() {
                if (!pcaState.assignments) return;
                var lines = ['sample_id\tcluster\tk\tpcs_used'];
                pcaData.forEach(function(p, idx) {
                    var label = pcaState.assignments[idx];
                    if (label === undefined) return;
                    lines.push([p.id, String(label + 1), String(pcaState.k), String(pcaState.pcsUsed)].join('\t'));
                });
                var blob = new Blob([lines.join('\n') + '\n'], {type:'text/tab-separated-values;charset=utf-8;'});
                var url = URL.createObjectURL(blob);
                var a = document.createElement('a');
                a.href = url;
                a.download = 'ancestry_cluster_assignments.tsv';
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                URL.revokeObjectURL(url);
            }

            function runClustering() {
                var kInput = document.getElementById('pca-cluster-k');
                var pcsInput = document.getElementById('pca-cluster-pcs');
                var statusEl = document.getElementById('pca-cluster-status');
                var exportBtn = document.getElementById('pca-cluster-export');
                if (!kInput || !pcsInput || !statusEl || !exportBtn) return;

                var k = parseInt(kInput.value, 10);
                var pcsUsed = parseInt(pcsInput.value, 10);
                if (!isFinite(k) || !isFinite(pcsUsed)) {
                    statusEl.textContent = 'Invalid clustering settings.';
                    return;
                }
                pcsUsed = Math.max(2, pcsUsed);
                var vectors = [];
                var idxMap = [];
                var availablePcCount = pcaData.reduce(function(acc, p) {
                    return Math.max(acc, (p.all_pcs || []).length);
                }, 0);
                if (availablePcCount && pcsUsed > availablePcCount) {
                    pcsUsed = availablePcCount;
                    pcsInput.value = String(pcsUsed);
                }
                pcaData.forEach(function(p, idx) {
                    var pcs = p.all_pcs || [];
                    if (!pcs || pcs.length < pcsUsed) return;
                    var vec = pcs.slice(0, pcsUsed);
                    if (vec.some(function(v){ return v === null || v === undefined; })) return;
                    vectors.push(vec);
                    idxMap.push(idx);
                });
                if (vectors.length < 2) {
                    statusEl.textContent = 'Not enough PCA points for clustering.';
                    return;
                }
                var result = kmeans(vectors, k, KMEANS_MAX_ITER);
                if (!result) {
                    statusEl.textContent = 'Clustering failed.';
                    return;
                }
                pcaState.assignments = new Array(pcaData.length);
                for (var i = 0; i < idxMap.length; i++) {
                    pcaState.assignments[idxMap[i]] = result.labels[i];
                }
                pcaState.mode = 'cluster';
                pcaState.k = result.k;
                pcaState.pcsUsed = pcsUsed;
                statusEl.textContent = 'Computed ' + result.k + ' clusters using first ' + pcsUsed + ' PCs.';
                exportBtn.disabled = false;
                var modeSel = document.getElementById('pca-color-mode');
                if (modeSel) modeSel.value = 'cluster';
                refreshPcaPlots();
            }

            refreshPcaPlots();

            var pcaSlider = document.getElementById('pca-density-bins');
            var pcaSliderVal = document.getElementById('pca-density-bins-value');
            if (pcaSlider) {
                pcaSlider.value = String(pcaState.binSize);
                if (pcaSliderVal) pcaSliderVal.textContent = String(pcaState.binSize);
                pcaSlider.addEventListener('input', function() {
                    var v = parseInt(this.value, 10);
                    if (!isFinite(v)) return;
                    if (pcaSliderVal) pcaSliderVal.textContent = String(v);
                    applyPcaDensityBins(v);
                });
            }

            var pcaScatterBtn = document.getElementById('pca-scatter-btn');
            var pcaDensityBtn = document.getElementById('pca-density-btn');
            if (pcaScatterBtn) pcaScatterBtn.addEventListener('click', function() { applyPcaViewMode('scatter'); });
            if (pcaDensityBtn) pcaDensityBtn.addEventListener('click', function() { applyPcaViewMode('density'); });

            var toggleBtn = document.getElementById('pca-view-toggle');
            var grid2d = document.getElementById('pca-2d-grid');
            var view3d = document.getElementById('pca-3d-container');
            if (toggleBtn && grid2d && view3d) {
                if (!pcaState.hasPc3) {
                    toggleBtn.disabled = true;
                    toggleBtn.textContent = '3D view unavailable';
                } else {
                    toggleBtn.addEventListener('click', function() {
                        var show3d = view3d.style.display === 'none';
                        view3d.style.display = show3d ? '' : 'none';
                        grid2d.style.display = show3d ? 'none' : '';
                        this.textContent = show3d ? 'Show 2D view' : 'Show 3D view';
                        this.classList.toggle('active', show3d);
                    });
                }
            }

            var colorSel = document.getElementById('pca-color-mode');
            if (colorSel) {
                colorSel.addEventListener('change', function() {
                    pcaState.mode = this.value === 'cluster' ? 'cluster' : 'sex';
                    refreshPcaPlots();
                });
            }
            var clusterBtn = document.getElementById('pca-cluster-run');
            if (clusterBtn) clusterBtn.addEventListener('click', runClustering);
            var exportBtn2 = document.getElementById('pca-cluster-export');
            if (exportBtn2) exportBtn2.addEventListener('click', exportClusters);

            /* Peddy ancestry color overlay: recolors pipeline samples by predicted ancestry.
               Uses peddy's ancestry predictions from het_check data matched by sample ID.
               The pipeline's own PC coordinates are preserved — only the coloring changes. */
            var overlayBtn = document.getElementById('peddy-overlay-toggle');
            if (overlayBtn) {
                overlayBtn.addEventListener('click', function() {
                    var isActive = this.classList.toggle('active');
                    pcaState.peddyOverlayActive = isActive;
                    this.textContent = isActive ? 'Hide peddy ancestry colors' : 'Color by peddy ancestry';
                    refreshPcaPlots();
                });
            }

            /* Peddy PCs section toggle: shows or hides the separate Peddy PCA section.
               Peddy PCs live in peddy's own coordinate space (1000G reference panel projection)
               and are NOT comparable to the pipeline PCA coordinate space. */
            var peddyPcsBtn = document.getElementById('peddy-pcs-toggle');
            if (peddyPcsBtn) {
                var peddySectionEl = document.getElementById('peddy');
                if (peddySectionEl) peddySectionEl.style.display = 'none';
                peddyPcsBtn.addEventListener('click', function() {
                    var hidden = !peddySectionEl || peddySectionEl.style.display === 'none';
                    if (peddySectionEl) peddySectionEl.style.display = hidden ? '' : 'none';
                    this.textContent = hidden ? 'Hide peddy PCs' : 'Show peddy PCs';
                    this.classList.toggle('active', hidden);
                });
            }
        }
    }

    /* ---- Interactive sex check ---- */
    var sexEl = document.getElementById('sex-data');
    if (sexEl) {
        var sexData = JSON.parse(sexEl.textContent || '[]');
        if (sexData && sexData.length) {
            var grouped2 = {};
            var allX2 = [];
            var allY2 = [];
            sexData.forEach(function(p) {
                if (p.chrx_lrr_median === null || p.chry_lrr_median === null) return;
                var sx = normSex(p.gender);
                var label = sx === 'M' ? 'Male' : (sx === 'F' ? 'Female' : 'Unknown');
                var color = sx === 'M' ? '#4393C3' : (sx === 'F' ? '#D6604D' : '#6b7280');
                if (!grouped2[label]) {
                    grouped2[label] = {x:[], y:[], text:[], color:color};
                }
                grouped2[label].x.push(p.chrx_lrr_median);
                grouped2[label].y.push(p.chry_lrr_median);
                grouped2[label].text.push(p.id);
                allX2.push(p.chrx_lrr_median);
                allY2.push(p.chry_lrr_median);
            });

            var traces2 = Object.keys(grouped2).sort().map(function(label) {
                var g = grouped2[label];
                return {
                    x:g.x, y:g.y, text:g.text, mode:'markers', type:'scatter',
                    name:label+' (n='+g.x.length+')',
                    marker:{color:g.color, size:5, opacity:0.65},
                    hovertemplate:'<b>%{text}</b><br>Median chrX LRR: %{x:.4f}<br>Median chrY LRR: %{y:.4f}<extra>'+label+'</extra>'
                };
            });
            if (traces2.length) {
                var densityIdx2 = traces2.length;
                traces2.push({
                    x:allX2, y:allY2, type:'histogram2d', name:'Density',
                    nbinsx:35, nbinsy:35,
                    colorscale:[
                        [0.0, 'rgba(37,99,235,0)'],
                        [0.05, 'rgba(37,99,235,0.2)'],
                        [0.4, 'rgba(37,99,235,0.5)'],
                        [1.0, 'rgba(30,64,175,0.85)']
                    ],
                    zmin:0, showscale:true, colorbar:{title:'Count'},
                    hovertemplate:'Median chrX LRR: %{x:.4f}<br>Median chrY LRR: %{y:.4f}<br>Count: %{z}<extra>Density</extra>',
                    visible:false
                });
                var sexDiv = document.getElementById('plot-sex-check');
                if (sexDiv) {
                    var sexTraceCount = traces2.length;
                    Plotly.newPlot(sexDiv, traces2, Object.assign({}, baseLayout, {
                        title:{text:'Sex Check: Median chrX vs chrY LRR',font:{size:14}},
                        xaxis:{title:'Median chrX LRR',gridcolor:'#f1f5f9',zeroline:false},
                        yaxis:{title:'Median chrY LRR',gridcolor:'#f1f5f9',zeroline:false}
                    }), cfg);

                    var sexSlider = document.getElementById('sex-density-bins');
                    var sexSliderVal = document.getElementById('sex-density-bins-value');
                    if (sexSlider) {
                        sexSlider.value = '35';
                        if (sexSliderVal) sexSliderVal.textContent = '35';
                        sexSlider.addEventListener('input', function() {
                            var bins = parseInt(this.value, 10);
                            if (!isFinite(bins)) return;
                            if (sexSliderVal) sexSliderVal.textContent = String(bins);
                            Plotly.restyle(sexDiv, {'nbinsx': bins, 'nbinsy': bins}, [densityIdx2]);
                        });
                    }

                    var sexScatterBtn = document.getElementById('sex-scatter-btn');
                    var sexDensityBtn = document.getElementById('sex-density-btn');
                    function setSexViewMode(mode) {
                        var visible = Array.from({length: sexTraceCount}, function(_, i) {
                            return mode === 'density' ? i === densityIdx2 : i !== densityIdx2;
                        });
                        Plotly.restyle(sexDiv, {visible: visible});
                        if (sexScatterBtn) sexScatterBtn.classList.toggle('active', mode === 'scatter');
                        if (sexDensityBtn) sexDensityBtn.classList.toggle('active', mode === 'density');
                    }
                    if (sexScatterBtn) sexScatterBtn.addEventListener('click', function() { setSexViewMode('scatter'); });
                    if (sexDensityBtn) sexDensityBtn.addEventListener('click', function() { setSexViewMode('density'); });
                }
            }
        }
    }

    /* ---- Per-sample table ---- */
    var tbody = document.getElementById('sample-tbody');
    if (tbody) {
        function cc(v,t,dir){
            if(v===null)return '';
            if(dir==='min')return v>=t?'cell-pass':'cell-fail';
            return v<=t?'cell-pass':'cell-fail';
        }
        function f4(v){return v!==null?v.toFixed(4):'NA';}

        function renderTable(filterMode) {
            tbody.innerHTML = '';
            var visible = 0;
            rawData.forEach(function(s) {
                if (filterMode === 'pass' && !s.qc_pass) return;
                if (filterMode === 'fail' && s.qc_pass) return;
                if (filterMode === 'flag' && !s.baf_flag && !s.het_outlier) return;
                var searchVal = (document.getElementById('sample-search') || {}).value || '';
                if (searchVal && s.id.toLowerCase().indexOf(searchVal.toLowerCase()) === -1) return;
                visible++;
                var tr = document.createElement('tr');
                var crC=cc(s.call_rate,0.97,'min'), lrC=cc(s.lrr_sd,0.35,'max');
                var bfC=s.baf_flag?'cell-warn':'';
                var htC=s.het_outlier?'cell-warn':'';
                var st=s.qc_pass?'<span class="status-badge status-pass">PASS</span>'
                                 :'<span class="status-badge status-fail">FAIL</span>';
                // Flags column: short badges for BAF and het-outlier warnings
                var flags = '';
                if (s.baf_flag) flags += '<span class="status-badge status-flag" style="margin-right:2px">BAF\u2191</span>';
                if (s.het_outlier) flags += '<span class="status-badge status-flag">Het</span>';
                tr.innerHTML=
                    '<td>'+s.id+'</td>'+
                    '<td class="'+crC+'">'+f4(s.call_rate)+'</td>'+
                    '<td class="'+lrC+'">'+f4(s.lrr_sd)+'</td>'+
                    '<td>'+f4(s.lrr_mean)+'</td>'+
                    '<td>'+f4(s.lrr_median)+'</td>'+
                    '<td class="'+bfC+'">'+f4(s.baf_sd)+'</td>'+
                    '<td class="'+htC+'">'+f4(s.het_rate)+'</td>'+
                    '<td>'+s.gender+'</td>'+
                    '<td>'+st+'</td>'+
                    '<td>'+flags+'</td>'+
                    '<td style="font-size:0.75rem;color:#dc2626">'+s.fail_reason+'</td>';
                tbody.appendChild(tr);
            });
            var sc = document.getElementById('shown-count');
            if (sc) sc.textContent = visible + ' of ' + rawData.length + ' samples';
        }

        var currentFilter = 'all';
        renderTable(currentFilter);

        /* ---- Filter dropdown ---- */
        var filterSel = document.getElementById('sample-filter');
        if (filterSel) filterSel.addEventListener('change', function() {
            currentFilter = this.value;
            renderTable(currentFilter);
        });
    }

    /* ---- Search ---- */
    var si=document.getElementById('sample-search'), sc=document.getElementById('shown-count');
    if(si&&tbody){si.addEventListener('input',function(){
        renderTable(currentFilter);
    });}

    /* ---- Download helpers ---- */
    function downloadIds(ids, filename) {
        var blob = new Blob([ids.join('\n') + '\n'], {type: 'text/plain'});
        var a = document.createElement('a');
        a.href = URL.createObjectURL(blob);
        a.download = filename;
        a.click();
        URL.revokeObjectURL(a.href);
    }
    var dlPass = document.getElementById('dl-pass');
    if (dlPass) dlPass.addEventListener('click', function() {
        downloadIds(rawData.filter(function(s){return s.qc_pass;}).map(function(s){return s.id;}), 'qc_pass_samples.txt');
    });
    var dlFail = document.getElementById('dl-fail');
    if (dlFail) dlFail.addEventListener('click', function() {
        downloadIds(rawData.filter(function(s){return !s.qc_pass;}).map(function(s){return s.id;}), 'qc_fail_samples.txt');
    });
    var dlFlagged = document.getElementById('dl-flagged');
    if (dlFlagged) dlFlagged.addEventListener('click', function() {
        downloadIds(rawData.filter(function(s){return s.baf_flag||s.het_outlier;}).map(function(s){return s.id;}), 'qc_flagged_samples.txt');
    });

    /* ---- Sort ---- */
    var ths=document.querySelectorAll('.sample-table th[data-col]');
    ths.forEach(function(th){th.addEventListener('click',function(){
        var c=parseInt(this.dataset.col), asc=this.dataset.asc!=='true';
        this.dataset.asc=asc;
        var rows=Array.from(tbody.getElementsByTagName('tr'));
        rows.sort(function(a,b){
            var va=a.cells[c].textContent,vb=b.cells[c].textContent;
            // strip badge html for status column
            if(va.indexOf('PASS')>-1)va='1'; if(va.indexOf('FAIL')>-1)va='0';
            if(vb.indexOf('PASS')>-1)vb='1'; if(vb.indexOf('FAIL')>-1)vb='0';
            var na=parseFloat(va),nb=parseFloat(vb);
            if(!isNaN(na)&&!isNaN(nb)) return asc?na-nb:nb-na;
            return asc?va.localeCompare(vb):vb.localeCompare(va);
        });
        rows.forEach(function(r){tbody.appendChild(r);});
        ths.forEach(function(h){h.classList.remove('sort-asc','sort-desc');});
        this.classList.add(asc?'sort-asc':'sort-desc');
    });});

    /* ---- Peddy Ancestry PCA ---- */
    var peddyEl = document.getElementById('peddy-data');
    if (peddyEl) {
        var peddyData = JSON.parse(peddyEl.textContent || '[]');
        if (peddyData && peddyData.length) {
            var PEDDY_ANCESTRY_COLORS = {
                'EUR': '#2563eb', 'AFR': '#dc2626', 'EAS': '#059669',
                'AMR': '#d97706', 'SAS': '#7c3aed', 'UNKNOWN': '#6b7280'
            };

            function peddyPcaPlot(divId, xKey, yKey, title) {
                var d = document.getElementById(divId);
                if (!d) return;
                var grouped = {};
                peddyData.forEach(function(p) {
                    if (p[xKey] === null || p[yKey] === null) return;
                    var anc = p.ancestry || 'UNKNOWN';
                    if (!grouped[anc]) grouped[anc] = {x:[], y:[], text:[], color: PEDDY_ANCESTRY_COLORS[anc] || '#6b7280'};
                    grouped[anc].x.push(p[xKey]);
                    grouped[anc].y.push(p[yKey]);
                    grouped[anc].text.push(p.id + ' (' + anc + ', p=' + (p.ancestry_prob || 0).toFixed(2) + ')');
                });
                var traces = Object.keys(grouped).sort().map(function(anc) {
                    var g = grouped[anc];
                    return {
                        x:g.x, y:g.y, text:g.text, mode:'markers', type:'scatter',
                        name:anc+' (n='+g.x.length+')',
                        marker:{color:g.color, size:6, opacity:0.7},
                        hovertemplate:'<b>%{text}</b><br>'+xKey.toUpperCase()+': %{x:.4f}<br>'+yKey.toUpperCase()+': %{y:.4f}<extra>'+anc+'</extra>'
                    };
                });
                if (!traces.length) return;
                Plotly.newPlot(d, traces, Object.assign({}, baseLayout, {
                    title:{text:title,font:{size:14}},
                    xaxis:{title:xKey.toUpperCase()+' (peddy)',gridcolor:'#f1f5f9',zeroline:false},
                    yaxis:{title:yKey.toUpperCase()+' (peddy)',gridcolor:'#f1f5f9',zeroline:false}
                }), cfg);
            }

            peddyPcaPlot('plot-peddy-pca12', 'pc1', 'pc2', 'Peddy Ancestry PCA: PC1 vs PC2');
            peddyPcaPlot('plot-peddy-pca34', 'pc3', 'pc4', 'Peddy Ancestry PCA: PC3 vs PC4');
        }
    }

    /* ---- Variant QC tabs ---- */
    var vqcTabs = document.querySelectorAll('.vqc-tab');
    if (vqcTabs.length) {
        vqcTabs.forEach(function(tab) {
            tab.addEventListener('click', function() {
                vqcTabs.forEach(function(t) { t.classList.remove('active'); });
                this.classList.add('active');
                var panels = document.querySelectorAll('.vqc-panel');
                panels.forEach(function(p) { p.style.display = 'none'; });
                var target = document.getElementById(this.dataset.tab);
                if (target) target.style.display = 'block';
            });
        });
    }

    /* ---- Variant QC interactive plots ---- */
    var vqcEl = document.getElementById('collated-vqc-data');
    if (vqcEl) {
        var vqcData = {};
        try { vqcData = JSON.parse(vqcEl.textContent || '{}'); } catch(e) {}
        Object.keys(vqcData).forEach(function(group) {
            var gd = vqcData[group];
            var label = group === 'all' ? 'All Samples' : group;

            /* Call rate histogram */
            var crDiv = document.getElementById('plot-vqc-cr-' + group);
            if (crDiv && gd.cr_values && gd.cr_values.length) {
                Plotly.newPlot(crDiv, [{
                    x: gd.cr_values, type: 'histogram',
                    marker: {color: C.blue, opacity: 0.7, line: {color: 'white', width: 0.5}},
                    hovertemplate: 'Call Rate: %{x}<br>Count: %{y}<extra></extra>'
                }], Object.assign({}, baseLayout, {
                    title: {text: 'Variant Call Rate (' + label + ')', font: {size: 13}},
                    xaxis: {title: 'Call Rate', gridcolor: '#f1f5f9'},
                    yaxis: {title: 'Count', gridcolor: '#f1f5f9'},
                    shapes: [{type:'line', x0:0.98, x1:0.98, y0:0, y1:1, yref:'paper',
                              line:{color:C.thresh, width:1.5, dash:'dash'}}],
                    annotations: [{x:0.98, y:1.04, yref:'paper', text:'\u2265 0.98',
                                   showarrow:false, font:{size:10, color:C.thresh}}]
                }), cfg);
            }

            /* MAF histogram */
            var mafDiv = document.getElementById('plot-vqc-maf-' + group);
            if (mafDiv && gd.maf_values && gd.maf_values.length) {
                Plotly.newPlot(mafDiv, [{
                    x: gd.maf_values, type: 'histogram',
                    marker: {color: '#059669', opacity: 0.7, line: {color: 'white', width: 0.5}},
                    hovertemplate: 'MAF: %{x}<br>Count: %{y}<extra></extra>'
                }], Object.assign({}, baseLayout, {
                    title: {text: 'MAF Distribution (' + label + ')', font: {size: 13}},
                    xaxis: {title: 'Minor Allele Frequency', gridcolor: '#f1f5f9'},
                    yaxis: {title: 'Count', gridcolor: '#f1f5f9'},
                    shapes: [{type:'line', x0:0.01, x1:0.01, y0:0, y1:1, yref:'paper',
                              line:{color:C.thresh, width:1.5, dash:'dash'}}],
                    annotations: [{x:0.01, y:1.04, yref:'paper', text:'MAF \u2265 0.01',
                                   showarrow:false, font:{size:10, color:C.thresh}}]
                }), cfg);
            }
        });
    }
});
"""


def _build_html(stats, stage1_stats, figures, realign_text,
                comparison_text, variant_qc_text, diag_text,
                methods_text, pca_summary, tool_versions, stage_label,
                qc_rows=None, pca_json='[]', sex_check_json='[]',
                peddy_json='[]', ancestry_vqc_texts=None,
                ancestry_strat_summary=None, collated_vqc_json='{}',
                mock_notice=''):
    """Build the HTML report string with interactive Plotly charts."""

    if ancestry_vqc_texts is None:
        ancestry_vqc_texts = {}

    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # --- Prepare interactive data ---
    sample_json = _prepare_sample_json(qc_rows, stats) if qc_rows else '[]'
    mock_notice_html = ''
    if mock_notice:
        escaped_notice = (mock_notice.strip().replace('&', '&amp;')
                          .replace('<', '&lt;')
                          .replace('>', '&gt;')
                          .replace('"', '&quot;')
                          .replace("'", '&#39;'))
        mock_notice_html = f'<div class="mock-banner"><strong>Example Output (Mock Data):</strong> {escaped_notice}</div>'

    # --- Static content helpers ---
    def _stat_row(label, metric, s):
        mean = s.get(f'{metric}_mean')
        median = s.get(f'{metric}_median')
        sd = s.get(f'{metric}_sd')
        if mean is None:
            return ''
        return (f'<tr><td>{label}</td>'
                f'<td>{mean:.4f}</td><td>{median:.4f}</td>'
                f'<td>{sd:.4f}</td>'
                f'<td>{s.get(f"{metric}_min", 0):.4f}</td>'
                f'<td>{s.get(f"{metric}_max", 0):.4f}</td></tr>')

    def _pre_block(text, title=''):
        if not text:
            return ''
        escaped = text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
        return f'<div class="card"><h3>{title}</h3><pre class="report-pre">{escaped}</pre></div>'

    def _fig_block(key, title=''):
        if key not in figures:
            return ''
        return (f'<div class="card figure-block"><h3>{title}</h3>'
                f'<img src="data:image/png;base64,{figures[key]}" '
                f'alt="{title}"></div>')

    # Build stats rows
    stats_rows = ''
    for label, metric in [
        ('Call Rate', 'call_rate'), ('LRR SD', 'lrr_sd'),
        ('LRR Mean', 'lrr_mean'), ('LRR Median', 'lrr_median'),
        ('BAF SD (het sites)', 'baf_sd'), ('Heterozygosity Rate', 'het_rate'),
    ]:
        stats_rows += _stat_row(label, metric, stats)

    # Stage comparison table
    comparison_table = ''
    if stage1_stats:
        comparison_table = '<div class="card"><h3>Stage 1 → Stage 2 Comparison</h3><table class="stats-table">'
        comparison_table += '<tr><th>Metric</th><th>Stage 1 Mean</th><th>Stage 2 Mean</th><th>Change</th></tr>'
        for label, metric in [('Call Rate', 'call_rate'), ('LRR SD', 'lrr_sd'),
                              ('BAF SD', 'baf_sd'), ('Het Rate', 'het_rate')]:
            s1 = stage1_stats.get(f'{metric}_mean')
            s2 = stats.get(f'{metric}_mean')
            if s1 is not None and s2 is not None:
                diff = s2 - s1
                good = ((diff > 0 and metric == 'call_rate')
                        or (diff < 0 and metric in ('lrr_sd', 'baf_sd')))
                cls = 'comparison-change-good' if good else 'comparison-change-bad'
                comparison_table += (
                    f'<tr><td>{label}</td><td>{s1:.4f}</td><td>{s2:.4f}</td>'
                    f'<td class="{cls}">{diff:+.4f}</td></tr>')
        comparison_table += '</table></div>'

    # Tool versions
    tools_html = '<ul class="tools-list">'
    for tool, ver in tool_versions.items():
        display_ver = ver if ver else 'not installed'
        tools_html += f'<li><code>{tool}</code>: {display_ver}</li>'
    tools_html += '</ul>'

    def _esc(text):
        return (text.replace('&', '&amp;')
                .replace('<', '&lt;')
                .replace('>', '&gt;'))

    citation_rows = ''
    for row in GWAS_METHODS_CITATIONS:
        doi = row['doi']
        citation_rows += (
            f"<tr><td>{_esc(row['citation'])}</td>"
            f"<td><a href=\"https://doi.org/{doi}\" target=\"_blank\" rel=\"noopener\">doi:{doi}</a></td>"
            f"<td>{_esc(row['cited_for'])}</td>"
            f"<td>{_esc(row['pipeline_application'])}</td></tr>"
        )

    # Pass rate
    pass_rate = stats.get('pass_rate', 0)
    pass_color = '#059669' if pass_rate >= 90 else '#d97706' if pass_rate >= 70 else '#dc2626'

    # --- CSS and JS (regular strings, no f-string escaping needed) ---
    css = _get_report_css()
    js = _get_report_js()

    # --- Assemble HTML ---
    # The HTML body uses f-strings; CSS/JS are inserted as pre-built strings.
    html_body = f"""
    <div class="header">
        <h1>🧬 Illumina IDAT Processing Pipeline Report</h1>
        <p class="subtitle">Generated: {now}</p>
    </div>

    <nav class="nav-bar">
        <a href="#overview">Overview</a>
        <a href="#thresholds">QC Thresholds</a>
        <a href="#qc-dashboard">QC Dashboard</a>
        <a href="#sample-table">Sample Details</a>
        <a href="#comparison">Stage Comparison</a>
        <a href="#variant-qc">Variant QC</a>
        <a href="#sex-check">Sex Check</a>
        <a href="#pca">Ancestry PCA</a>
        <a href="#peddy">Peddy QC</a>
        <a href="#realignment">Realignment</a>
        <a href="#diagnostics">Diagnostics</a>
        <a href="#citations">Citations</a>
        <a href="#methods">Methods</a>
        <a href="#tools">Tools</a>
    </nav>

    <div class="container">
    {mock_notice_html}

    <section id="overview">
        <h2>📊 Overview</h2>
        <div class="metric-grid">
            <div class="metric-card">
                <div class="value">{stats.get('n_samples', 0)}</div>
                <div class="label">Total Samples</div>
            </div>
            <div class="metric-card pass">
                <div class="value">{stats.get('n_pass', 0)}</div>
                <div class="label">Passing QC</div>
            </div>
            <div class="metric-card fail">
                <div class="value">{stats.get('n_fail', 0)}</div>
                <div class="label">Failing QC</div>
            </div>
            <div class="metric-card">
                <div class="value">{stats.get('n_male', 0)} / {stats.get('n_female', 0)}</div>
                <div class="label">Male / Female</div>
            </div>
            <div class="metric-card">
                <div class="value">{stats.get('call_rate_mean', 0):.4f}</div>
                <div class="label">Mean Call Rate</div>
            </div>
            <div class="metric-card">
                <div class="value">{stats.get('lrr_sd_mean', 0):.4f}</div>
                <div class="label">Mean LRR SD</div>
            </div>
            <div class="metric-card warn">
                <div class="value">{stats.get('n_baf_flag', 0)}</div>
                <div class="label">BAF SD Flags</div>
            </div>
            <div class="metric-card warn">
                <div class="value">{stats.get('n_het_outlier', 0)}</div>
                <div class="label">Het Rate Outliers</div>
            </div>
        </div>

        <div class="card" style="max-width:400px">
            <strong>QC Pass Rate</strong>
            <div class="pass-bar-track">
                <div class="pass-bar-fill" style="width:{pass_rate:.1f}%;background:linear-gradient(90deg,{pass_color},{pass_color}dd)"></div>
            </div>
            <span style="font-size:0.85rem;font-weight:600;color:{pass_color}">{pass_rate:.1f}%</span>
            <span style="font-size:0.75rem;color:#64748b"> ({stats.get('n_pass',0)} of {stats.get('n_samples',0)})</span>
        </div>

        <div class="card">
            <h3>Summary Statistics ({stage_label})</h3>
            <table class="stats-table">
                <tr><th>Metric</th><th>Mean</th><th>Median</th><th>SD</th><th>Min</th><th>Max</th></tr>
                {stats_rows}
            </table>
        </div>

        {comparison_table}
    </section>

    <section id="thresholds">
        <h2>🎯 GWAS QC Best Practice Thresholds</h2>
        <div class="card">
            <p style="font-size:0.82rem;color:#64748b;margin-bottom:0.75rem">
                Reference thresholds based on Anderson et al. (2010), Marees et al. (2018),
                and Turner et al. (2011). These defaults are used to highlight samples and
                variants requiring attention, and should be interpreted with study design,
                ancestry composition, and downstream analysis goals in mind.
            </p>
            <table class="threshold-table">
                <tr><th>Metric</th><th>Threshold</th><th>Type</th><th>Rationale / Evidence</th></tr>
                <tr>
                    <td>Sample Call Rate</td><td>≥ 0.97</td>
                    <td><span class="badge badge-pass">Pass/Fail</span></td>
                    <td>Low call rate indicates poor DNA quality or failed hybridization (Anderson 2010; Marees 2018)</td>
                </tr>
                <tr>
                    <td>LRR SD</td><td>≤ 0.35</td>
                    <td><span class="badge badge-pass">Pass/Fail</span></td>
                    <td>High LRR SD indicates noisy intensity data unsuitable for CNV/GWAS (Turner 2011)</td>
                </tr>
                <tr>
                    <td>BAF SD (het sites)</td><td>≤ 0.15</td>
                    <td><span class="badge badge-warn">Flag</span></td>
                    <td>Elevated BAF SD may indicate contamination or noisy heterozygous intensity clusters (Turner 2011; Marees 2018)</td>
                </tr>
                <tr>
                    <td>Heterozygosity Rate</td><td>Within ± 3 SD</td>
                    <td><span class="badge badge-warn">Flag</span></td>
                    <td>Outliers may indicate contamination (high) or inbreeding (low); interpret against cohort ancestry structure (Anderson 2010; Marees 2018)</td>
                </tr>
                <tr>
                    <td>Variant Call Rate</td><td>≥ 0.98</td>
                    <td><span class="badge badge-info">Variant QC</span></td>
                    <td>Poorly genotyped variants produce spurious associations (Anderson 2010; Turner 2011)</td>
                </tr>
                <tr>
                    <td>HWE p-value</td><td>≥ 1 × 10⁻⁶</td>
                    <td><span class="badge badge-info">Variant QC</span></td>
                    <td>Extreme HWE deviation indicates genotyping errors (Anderson 2010; Marees 2018)</td>
                </tr>
                <tr>
                    <td>MAF</td><td>≥ 0.01</td>
                    <td><span class="badge badge-info">Variant QC</span></td>
                    <td>Rare variants have low power and higher error sensitivity in GWAS; PCA step uses MAF ≥ 0.05 for stable loadings (Marees 2018)</td>
                </tr>
                <tr>
                    <td>Inbreeding F</td><td>|F| ≤ 0.05</td>
                    <td><span class="badge badge-warn">Flag</span></td>
                    <td>F-statistic outliers may indicate sample quality or relatedness issues (Turner 2011)</td>
                </tr>
            </table>
            <p style="font-size:0.78rem;color:#64748b;margin-top:0.65rem">
                Key references: Anderson et al. 2010 (doi:10.1038/nprot.2010.116),
                Marees et al. 2018 (doi:10.1002/mpr.1608),
                Turner et al. 2011 (doi:10.1002/0471142905.hg0119s68).
            </p>
        </div>
    </section>

    <section id="qc-dashboard">
        <h2>📈 Interactive Sample QC Dashboard</h2>
        <div class="card">
            <div id="plot-scatter" class="plot-box" style="min-height:400px"></div>
        </div>
        <div class="plot-grid">
            <div class="card"><div id="plot-cr" class="plot-box"></div></div>
            <div class="card"><div id="plot-lrr" class="plot-box"></div></div>
        </div>
        <div class="plot-grid" id="baf-container" style="display:none">
            <div class="card"><div id="plot-baf" class="plot-box"></div></div>
            <div class="card" id="het-container" style="display:none"><div id="plot-het" class="plot-box"></div></div>
        </div>
        <noscript>
            {_fig_block('qc_dashboard', 'Sample QC Dashboard (static)')}
        </noscript>
    </section>

    <section id="sample-table">
        <h2>📋 Per-Sample QC Details</h2>
        <div class="card">
            <div class="sample-controls">
                <input type="text" id="sample-search" class="search-input" placeholder="Search by sample ID…">
                <select id="sample-filter" class="plot-input" style="font-size:0.8rem">
                    <option value="all">Show all samples</option>
                    <option value="pass">Passing QC only</option>
                    <option value="fail">Failing QC only</option>
                    <option value="flag">Flagged (BAF/Het) only</option>
                </select>
                <span class="sample-count" id="shown-count">{stats.get('n_samples',0)} of {stats.get('n_samples',0)} samples</span>
                <button id="dl-pass" class="plot-btn" type="button" title="Download list of passing sample IDs">⬇ Pass list</button>
                <button id="dl-fail" class="plot-btn" type="button" title="Download list of failing sample IDs">⬇ Fail list</button>
                <button id="dl-flagged" class="plot-btn" type="button" title="Download list of flagged sample IDs (BAF or Het outlier)">⬇ Flagged list</button>
            </div>
            <div class="sample-table-wrap">
                <table class="sample-table">
                    <thead><tr>
                        <th data-col="0">Sample ID</th>
                        <th data-col="1">Call Rate</th>
                        <th data-col="2">LRR SD</th>
                        <th data-col="3">LRR Mean</th>
                        <th data-col="4">LRR Median</th>
                        <th data-col="5">BAF SD</th>
                        <th data-col="6">Het Rate</th>
                        <th data-col="7">Sex</th>
                        <th data-col="8">Status</th>
                        <th data-col="9">Flags</th>
                        <th data-col="10">Failure Reason</th>
                    </tr></thead>
                    <tbody id="sample-tbody"></tbody>
                </table>
            </div>
        </div>
    </section>

    <section id="comparison">
        <h2>🔄 Reclustering Comparison</h2>
        {_fig_block('comparison_samples', 'Sample-Level QC')}
        {_fig_block('comparison_variants', 'Variant-Level QC')}
        {_pre_block(comparison_text, 'Comparison Report')}
    </section>

    <section id="variant-qc">
        <h2>🧪 Variant QC</h2>
        <div class="card">
            <p style="font-size:0.82rem;color:#64748b;margin-bottom:0.75rem">
                Variant-level QC metrics are shown for the full cohort ("All") and for each
                ancestry group meeting the minimum sample threshold. Ancestry-stratified HWE
                testing avoids inflated deviation from the Wahlund effect in mixed-ancestry
                samples (Anderson et al. 2010; Marees et al. 2018).
            </p>
            <div class="vqc-tabs" id="vqc-tab-bar">
                <button class="vqc-tab active" data-tab="vqc-all" type="button">All</button>
                {''.join(f'<button class="vqc-tab" data-tab="vqc-{anc}" type="button">{anc}</button>' for anc in sorted(ancestry_vqc_texts.keys()))}
            </div>
        </div>
        <div id="vqc-all" class="vqc-panel" style="display:block">
            {_pre_block(variant_qc_text, 'Variant QC Summary (All Samples)')}
            <div class="plot-grid">
                <div class="card"><div id="plot-vqc-cr-all" class="plot-box"></div></div>
                <div class="card"><div id="plot-vqc-maf-all" class="plot-box"></div></div>
            </div>
        </div>
        {''.join(
            f"""<div id="vqc-{anc}" class="vqc-panel" style="display:none">
            {_pre_block(text, f'Variant QC Summary ({anc})')}
            <div class="plot-grid">
                <div class="card"><div id="plot-vqc-cr-{anc}" class="plot-box"></div></div>
                <div class="card"><div id="plot-vqc-maf-{anc}" class="plot-box"></div></div>
            </div>
        </div>"""
            for anc, text in sorted(ancestry_vqc_texts.items())
        )}
        {_pre_block(ancestry_strat_summary, 'Ancestry-Stratified QC Summary') if ancestry_strat_summary else ''}
    </section>

    <section id="sex-check">
        <h2>⚥ Sex Check</h2>
        <div class="card" style="margin-bottom:0.75rem;font-size:0.83rem;color:#374151">
            <strong>Interpretation guide:</strong> Females (XX) show higher median chrX LRR (near 0)
            and low chrY LRR (negative, reflecting no chrY). Males (XY) show lower median chrX LRR
            (slightly negative, hemizygous X) and higher chrY LRR (near 0 or slightly positive).
            Outliers may indicate sex discordance, sex chromosome aneuploidies (XXY: high X and Y;
            XO: very low X with absent Y signal), or sample swaps. Samples outside the expected
            sex-specific clusters warrant manual review.
        </div>
        <div class="plot-controls">
            <div class="btn-group">
                <button id="sex-scatter-btn" class="plot-btn active" type="button">Scatter</button>
                <button id="sex-density-btn" class="plot-btn" type="button">Density</button>
            </div>
            <div class="control-sep"></div>
            <label for="sex-density-bins">Bins: <span id="sex-density-bins-value">35</span></label>
            <input id="sex-density-bins" class="plot-slider" type="range" min="10" max="80" value="35" step="1">
        </div>
        <div class="card"><div id="plot-sex-check" class="plot-box"></div></div>
        <noscript>
        {_fig_block('sex_check', 'Median chrX vs chrY LRR by Predicted Sex')}
        </noscript>
    </section>

    <section id="pca">
        <h2>🌍 Ancestry PCA</h2>
        <div class="plot-controls">
            <div class="btn-group">
                <button id="pca-scatter-btn" class="plot-btn active" type="button">Scatter</button>
                <button id="pca-density-btn" class="plot-btn" type="button">Density</button>
            </div>
            <label for="pca-density-bins">Bins: <span id="pca-density-bins-value">35</span></label>
            <input id="pca-density-bins" class="plot-slider" type="range" min="10" max="80" value="35" step="1">
            <div class="control-sep"></div>
            <button id="pca-view-toggle" class="plot-btn" type="button">Show 3D view</button>
            <div class="control-sep"></div>
            <label for="pca-color-mode">Color by</label>
            <select id="pca-color-mode" class="plot-input">
                <option value="sex">Sex</option>
                <option value="cluster">Cluster</option>
            </select>
            <div class="control-sep"></div>
            <label for="pca-cluster-k">k</label>
            <input id="pca-cluster-k" class="plot-input" type="number" min="2" max="12" value="3">
            <label for="pca-cluster-pcs">PCs</label>
            <input id="pca-cluster-pcs" class="plot-input" type="number" min="2" max="20" value="4">
            <button id="pca-cluster-run" class="plot-btn" type="button">Run k-means</button>
            <button id="pca-cluster-export" class="plot-btn" type="button" disabled>Export clusters</button>
            <div class="control-sep"></div>
            <button id="peddy-overlay-toggle" class="plot-btn" type="button">Color by peddy ancestry</button>
            <button id="peddy-pcs-toggle" class="plot-btn" type="button">Show peddy PCs</button>
            <span id="pca-cluster-status" class="plot-status">No clusters computed.</span>
        </div>
        <div class="plot-grid" id="pca-2d-grid">
            <div class="card"><div id="plot-pca12" class="plot-box"></div></div>
            <div class="card" id="pca34-container" style="display:none"><div id="plot-pca34" class="plot-box"></div></div>
        </div>
        <div class="card" id="pca-3d-container" style="display:none"><div id="plot-pca3d" class="plot-box"></div></div>
        <div class="plot-controls-note">Toggle Scatter/Density views, switch 2D/3D, and run customizable k-means clustering for ancestry grouping. &ldquo;Color by peddy ancestry&rdquo; recolors the pipeline&rsquo;s own samples by their peddy-predicted ancestry (no coordinate mixing). &ldquo;Show peddy PCs&rdquo; reveals a separate section with peddy&rsquo;s own PCA plots (peddy coordinate space, not comparable to pipeline PCs).</div>
        <noscript>
            {_fig_block('pca', 'Principal Components (colored by predicted sex)')}
        </noscript>
        {_pre_block(pca_summary, 'PCA QC Summary')}
    </section>

    <section id="peddy">
        <h2>🔬 Peddy QC (Pedigree/Sex/Ancestry)</h2>
        <p style="font-size:0.85rem;color:#64748b;margin-bottom:0.75rem">
            Peddy validates sex, ancestry, and relatedness from VCF genotypes.
            Ancestry predictions are colored by predicted population (EUR, AFR, EAS, AMR, SAS).
            <strong>Note:</strong> Peddy PCs are projected onto the 1000 Genomes reference panel
            and are in a separate coordinate space from the pipeline&rsquo;s own PCA &mdash;
            the two are not directly comparable.
        </p>
        <div class="plot-grid">
            <div class="card"><div id="plot-peddy-pca12" class="plot-box"></div></div>
            <div class="card"><div id="plot-peddy-pca34" class="plot-box"></div></div>
        </div>
        <div class="plot-controls-note">
            Peddy ancestry PCA projections colored by predicted population (peddy&rsquo;s own coordinate space).
            Use &ldquo;Color by peddy ancestry&rdquo; in the Ancestry PCA section above to recolor
            the pipeline&rsquo;s own samples by their peddy ancestry predictions without mixing coordinate spaces.
        </div>
    </section>

    <section id="realignment">
        <h2>🧭 Manifest Realignment</h2>
        {_pre_block(realign_text, 'Realignment Summary')}
    </section>

    <section id="diagnostics">
        <h2>🔍 QC Diagnostics</h2>
        {_pre_block(diag_text, 'Diagnostic Report')}
    </section>

    <section id="citations">
        <h2>📚 GWAS Methods and Best-Practice Citations</h2>
        <div class="card">
            <p style="font-size:0.82rem;color:#64748b;margin-bottom:0.75rem">
                Compiled citation summary for methods, tools, and algorithms used in this
                pipeline. The same content is also exported as <code>citations_summary.tsv</code>.
            </p>
            <table class="threshold-table">
                <tr><th>Citation</th><th>DOI</th><th>Cited for</th><th>Pipeline application</th></tr>
                {citation_rows}
            </table>
        </div>
    </section>

    <section id="methods">
        <h2>📝 Methods Text (for publications)</h2>
        <div class="methods-text">
            <p>{methods_text.replace(chr(10)+chr(10), '</p><p>')}</p>
        </div>
    </section>

    <section id="tools">
        <h2>🛠️ Tool Versions</h2>
        <div class="card">
            {tools_html}
        </div>
    </section>

    </div><!-- /container -->

    <div class="footer">
        Report generated by the Illumina IDAT Processing Pipeline
    </div>
"""

    # Assemble complete document
    html = '<!DOCTYPE html>\n<html lang="en">\n<head>\n'
    html += '    <meta charset="UTF-8">\n'
    html += '    <meta name="viewport" content="width=device-width, initial-scale=1.0">\n'
    html += '    <title>Illumina IDAT Pipeline Report</title>\n'
    html += '    <script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>\n'
    html += f'    <style>{css}</style>\n'
    html += '</head>\n<body>\n'
    html += html_body
    html += f'\n<script type="application/json" id="qc-data">{sample_json}</script>\n'
    html += f'<script type="application/json" id="pca-data">{pca_json}</script>\n'
    html += f'<script type="application/json" id="sex-data">{sex_check_json}</script>\n'
    html += f'<script type="application/json" id="peddy-data">{peddy_json}</script>\n'
    html += f'<script type="application/json" id="collated-vqc-data">{collated_vqc_json}</script>\n'
    html += f'<script>{js}</script>\n'
    html += '</body>\n</html>'

    return html


def main():
    parser = argparse.ArgumentParser(
        description="Generate comprehensive HTML pipeline report"
    )
    parser.add_argument("--output-dir", required=True,
                        help="Pipeline output directory")
    parser.add_argument("--genome", default=None,
                        help="Genome build (CHM13, GRCh38, GRCh37). "
                             "Auto-detected from reference/ dir if not provided.")
    args = parser.parse_args()

    if not os.path.isdir(args.output_dir):
        print(f"Error: Output directory not found: {args.output_dir}",
              file=sys.stderr)
        sys.exit(1)

    report_path = generate_html_report(args.output_dir, genome=args.genome)

    if report_path:
        print(f"Pipeline report: {report_path}")
        print(f"Summary stats:   {os.path.join(args.output_dir, 'summary_statistics.tsv')}")
        print(f"Methods text:    {os.path.join(args.output_dir, 'methods_text.txt')}")
        print(f"Citations TSV:   {os.path.join(args.output_dir, 'citations_summary.tsv')}")
        print(f"Summary dir:     {os.path.join(args.output_dir, 'summary')}")
    else:
        print("Error: Could not generate report.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
