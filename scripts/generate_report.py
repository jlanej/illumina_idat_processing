#!/usr/bin/env python3
"""
generate_report.py

Generate a comprehensive HTML pipeline report with publication-quality
figures, summary statistics, and auto-generated methods text.

This is the final step of the pipeline, producing a single self-contained
HTML report with:
  - Pipeline configuration and runtime summary
  - Sample QC dashboard (call rate vs LRR SD, BAF SD, het rate distributions)
  - Manifest realignment summary
  - Reclustering comparison (if Stage 2 was run)
  - Sex check visualization
  - PCA scatter plots (PC1 vs PC2, PC3 vs PC4)
  - Variant QC summary
  - Auto-generated methods paragraph for publications
  - Summary statistics table in publication-ready format

Usage:
    python3 generate_report.py --output-dir /path/to/pipeline_output

All input files are auto-detected from the standard output directory structure.
"""

import argparse
import base64
import io
import os
import subprocess
import sys
import statistics
from datetime import datetime


def get_tool_versions():
    """Collect versions of key tools."""
    versions = {}
    for tool, cmd in [
        ('bcftools', ['bcftools', '--version']),
        ('plink2', ['plink2', '--version']),
        ('bwa', ['bwa']),  # bwa prints version to stderr
        ('flashpca', ['flashpca', '--version']),
    ]:
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
            output = (proc.stdout + proc.stderr).strip().split('\n')[0]
            versions[tool] = output
        except Exception:
            versions[tool] = 'not available'
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
    except ImportError:
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


def create_pca_plot(pca_file, qc_rows):
    """Create PCA scatter plots (PC1 vs PC2, PC3 vs PC4) colored by predicted sex."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    if not os.path.exists(pca_file):
        return None

    # Read PCA projections
    pca_data = {}
    with open(pca_file) as f:
        header = f.readline().strip().split('\t')
        if len(header) < 4:
            header = f.readline().strip().split()  # whitespace-separated
        for line in f:
            fields = line.strip().split('\t') if '\t' in line else line.strip().split()
            if len(fields) < 4:
                continue
            # IID is second column (FID, IID, PC1, PC2, ...)
            sample_id = fields[1] if len(fields) > 1 else fields[0]
            try:
                pcs = [float(x) for x in fields[2:]]
                pca_data[sample_id] = pcs
            except (ValueError, IndexError):
                pass

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
                 if safe_float(r.get('call_rate'), 0) >= 0.97 and
                 safe_float(r.get('lrr_sd'), 999) <= 0.35)
    stats['n_pass'] = n_pass
    stats['n_fail'] = n - n_pass

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

    bcftools_ver = tool_versions.get('bcftools', 'bcftools')
    plink2_ver = tool_versions.get('plink2', 'plink2')

    text = (
        f"Genotyping intensity data (IDAT files) were processed using the "
        f"Illumina IDAT Processing Pipeline. Raw IDAT files were converted to "
        f"GTC format using the GenCall algorithm (bcftools +idat2gtc, "
        f"{bcftools_ver.split()[0] if bcftools_ver else 'bcftools'} "
        f"{bcftools_ver.split()[1] if len(bcftools_ver.split()) > 1 else ''}), "
        f"then to VCF format with B Allele Frequency (BAF) and Log R Ratio (LRR) "
        f"intensities (bcftools +gtc2vcf). "
        f"Probe coordinates were validated by realigning CSV manifest flank "
        f"sequences against the {genome} reference genome using BWA-MEM. "
        f"\n\n"
        f"A two-stage genotyping approach was employed. Stage 1 genotype calls "
        f"were generated using manufacturer-provided EGT cluster definitions. "
        f"Per-sample QC metrics (call rate, LRR standard deviation, BAF standard "
        f"deviation at heterozygous sites, and heterozygosity rate) were computed "
        f"on autosomal variants. High-quality samples (call rate ≥ 0.97 and "
        f"LRR SD ≤ 0.35; n = {n_pass} of {n}) were used to recompute "
        f"study-specific genotype cluster definitions (EGT file) in Stage 2. "
        f"All samples were then re-genotyped using the study-specific clusters, "
        f"with additional BAF/LRR median adjustment (--adjust-clusters). "
        f"\n\n"
        f"Variant-level QC included per-variant missingness, Hardy-Weinberg "
        f"equilibrium tests (mid-p adjustment), allele frequency estimation, "
        f"per-sample inbreeding coefficients (F statistic), and transition/"
        f"transversion ratio computation ({plink2_ver.split()[0] if plink2_ver else 'plink2'}). "
        f"Ancestry principal components were computed using stringent variant QC "
        f"(missingness < 2%, HWE p ≥ 1e-6, MAF ≥ 5%), LD pruning "
        f"(window = 1000 kb, step = 100, r² < 0.1), and flashpca2. "
        f"After processing, the mean call rate was {cr:.4f} and the mean "
        f"LRR SD was {sd:.4f}."
    )
    return text


def generate_html_report(output_dir):
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

    # Diagnostic report
    diag_text = read_text(os.path.join(output_dir, 'qc_diagnostic_report.txt'))

    # PCA file
    pca_file = os.path.join(output_dir, 'ancestry_pca', 'pca_projections.tsv')
    pca_summary = read_text(os.path.join(output_dir, 'ancestry_pca', 'qc_summary.txt'))

    # Generate methods text
    genome = 'CHM13'  # Could parse from config, but CHM13 is default
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

    except ImportError:
        print("Warning: matplotlib not available. Report will have no figures.",
              file=sys.stderr)

    # ---- Build HTML ----
    html = _build_html(stats, stage1_stats,
                       figures, realign_text, comparison_text,
                       variant_qc_text, diag_text, methods_text,
                       pca_summary, tool_versions, stage_label)

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

    return report_path


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


def _build_html(stats, stage1_stats, figures, realign_text,
                comparison_text, variant_qc_text, diag_text,
                methods_text, pca_summary, tool_versions, stage_label):
    """Build the HTML report string."""

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
        return f'<h3>{title}</h3><pre class="report-pre">{escaped}</pre>'

    def _fig_block(key, title=''):
        if key not in figures:
            return ''
        return (f'<h3>{title}</h3>'
                f'<img src="data:image/png;base64,{figures[key]}" '
                f'style="max-width:100%; border:1px solid #ddd; border-radius:4px;">')

    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

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
        comparison_table = '<h3>Stage 1 → Stage 2 Comparison</h3><table class="stats-table">'
        comparison_table += '<tr><th>Metric</th><th>Stage 1 Mean</th><th>Stage 2 Mean</th><th>Change</th></tr>'
        for label, metric in [('Call Rate', 'call_rate'), ('LRR SD', 'lrr_sd')]:
            s1 = stage1_stats.get(f'{metric}_mean')
            s2 = stats.get(f'{metric}_mean')
            if s1 is not None and s2 is not None:
                diff = s2 - s1
                color = 'green' if (diff > 0 and metric == 'call_rate') or (diff < 0 and metric == 'lrr_sd') else 'red'
                comparison_table += (
                    f'<tr><td>{label}</td><td>{s1:.4f}</td><td>{s2:.4f}</td>'
                    f'<td style="color:{color};font-weight:bold">{diff:+.4f}</td></tr>'
                )
        comparison_table += '</table>'

    # Tool versions
    tools_html = '<ul>'
    for tool, ver in tool_versions.items():
        tools_html += f'<li><code>{tool}</code>: {ver}</li>'
    tools_html += '</ul>'

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Illumina IDAT Pipeline Report</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
               max-width: 1200px; margin: 0 auto; padding: 20px; color: #333;
               line-height: 1.6; background: #f8f9fa; }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #2c3e50; border-bottom: 1px solid #bdc3c7; padding-bottom: 5px;
              margin-top: 40px; }}
        h3 {{ color: #34495e; }}
        .summary-box {{ background: white; border: 1px solid #ddd; border-radius: 8px;
                        padding: 20px; margin: 15px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .stats-table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
        .stats-table th, .stats-table td {{ border: 1px solid #ddd; padding: 8px 12px;
                                            text-align: right; }}
        .stats-table th {{ background: #3498db; color: white; font-weight: 600; }}
        .stats-table td:first-child {{ text-align: left; font-weight: 500; }}
        .stats-table tr:nth-child(even) {{ background: #f2f2f2; }}
        .stats-table tr:hover {{ background: #e8f4f8; }}
        .report-pre {{ background: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 6px;
                       overflow-x: auto; font-size: 12px; line-height: 1.4; }}
        .methods-text {{ background: white; border-left: 4px solid #3498db; padding: 15px 20px;
                         margin: 15px 0; font-style: italic; border-radius: 0 8px 8px 0; }}
        .highlight {{ background: #fff3cd; padding: 2px 6px; border-radius: 3px; }}
        .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                        gap: 15px; margin: 15px 0; }}
        .metric-card {{ background: white; border: 1px solid #ddd; border-radius: 8px;
                        padding: 15px; text-align: center; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }}
        .metric-card .value {{ font-size: 24px; font-weight: 700; color: #2c3e50; }}
        .metric-card .label {{ font-size: 12px; color: #7f8c8d; margin-top: 5px; }}
        .toc {{ background: white; border: 1px solid #ddd; border-radius: 8px;
                padding: 15px 25px; margin: 20px 0; }}
        .toc a {{ text-decoration: none; color: #3498db; }}
        .toc a:hover {{ text-decoration: underline; }}
        .timestamp {{ color: #999; font-size: 12px; }}
        img {{ margin: 10px 0; }}
    </style>
</head>
<body>
    <h1>🧬 Illumina IDAT Processing Pipeline Report</h1>
    <p class="timestamp">Generated: {now}</p>

    <div class="toc">
        <strong>Contents:</strong>
        <a href="#overview">Overview</a> ·
        <a href="#qc-dashboard">QC Dashboard</a> ·
        <a href="#comparison">Stage Comparison</a> ·
        <a href="#variant-qc">Variant QC</a> ·
        <a href="#sex-check">Sex Check</a> ·
        <a href="#pca">Ancestry PCA</a> ·
        <a href="#realignment">Manifest Realignment</a> ·
        <a href="#diagnostics">Diagnostics</a> ·
        <a href="#methods">Methods</a> ·
        <a href="#tools">Tools</a>
    </div>

    <h2 id="overview">Overview</h2>
    <div class="metric-grid">
        <div class="metric-card">
            <div class="value">{stats.get('n_samples', 0)}</div>
            <div class="label">Total Samples</div>
        </div>
        <div class="metric-card">
            <div class="value">{stats.get('n_pass', 0)}</div>
            <div class="label">Passing QC</div>
        </div>
        <div class="metric-card">
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
    </div>

    <div class="summary-box">
        <h3>Summary Statistics ({stage_label})</h3>
        <table class="stats-table">
            <tr><th>Metric</th><th>Mean</th><th>Median</th><th>SD</th><th>Min</th><th>Max</th></tr>
            {stats_rows}
        </table>
    </div>

    {comparison_table}

    <h2 id="qc-dashboard">Sample QC Dashboard</h2>
    {_fig_block('qc_dashboard', '')}

    <h2 id="comparison">Reclustering Comparison</h2>
    {_fig_block('comparison_samples', 'Sample-Level QC')}
    {_fig_block('comparison_variants', 'Variant-Level QC')}
    {_pre_block(comparison_text, 'Comparison Report')}

    <h2 id="variant-qc">Variant QC</h2>
    {_pre_block(variant_qc_text, 'Variant QC Summary')}

    <h2 id="sex-check">Sex Check</h2>
    {_fig_block('sex_check', 'Median chrX vs chrY LRR by Predicted Sex')}

    <h2 id="pca">Ancestry PCA</h2>
    {_fig_block('pca', 'Principal Components (colored by predicted sex)')}
    {_pre_block(pca_summary, 'PCA QC Summary')}

    <h2 id="realignment">Manifest Realignment</h2>
    {_pre_block(realign_text, 'Realignment Summary')}

    <h2 id="diagnostics">QC Diagnostics</h2>
    {_pre_block(diag_text, 'Diagnostic Report')}

    <h2 id="methods">Methods Text (for publications)</h2>
    <div class="methods-text">
        <p>{methods_text.replace(chr(10)+chr(10), '</p><p>')}</p>
    </div>

    <h2 id="tools">Tool Versions</h2>
    {tools_html}

    <hr>
    <p class="timestamp">Report generated by the Illumina IDAT Processing Pipeline</p>
</body>
</html>"""

    return html


def main():
    parser = argparse.ArgumentParser(
        description="Generate comprehensive HTML pipeline report"
    )
    parser.add_argument("--output-dir", required=True,
                        help="Pipeline output directory")
    args = parser.parse_args()

    if not os.path.isdir(args.output_dir):
        print(f"Error: Output directory not found: {args.output_dir}",
              file=sys.stderr)
        sys.exit(1)

    report_path = generate_html_report(args.output_dir)

    if report_path:
        print(f"Pipeline report: {report_path}")
        print(f"Summary stats:   {os.path.join(args.output_dir, 'summary_statistics.tsv')}")
        print(f"Methods text:    {os.path.join(args.output_dir, 'methods_text.txt')}")
    else:
        print("Error: Could not generate report.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

