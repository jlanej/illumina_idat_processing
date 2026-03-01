#!/usr/bin/env python3
"""
plot_qc_comparison.py

Generate diagnostic plots and summary comparing QC metrics before and after
reclustering. Compares:
  - Per-sample call rate distributions
  - Per-sample LRR SD distributions
  - Variant-level QC (missingness, HWE, MAF) summaries

Usage:
    python3 plot_qc_comparison.py \
        --stage1-sample-qc stage1/qc/stage1_sample_qc.tsv \
        --stage2-sample-qc stage2/qc/stage2_sample_qc.tsv \
        --output-dir comparison_output

    # With variant QC:
    python3 plot_qc_comparison.py \
        --stage1-sample-qc stage1/qc/stage1_sample_qc.tsv \
        --stage2-sample-qc stage2/qc/stage2_sample_qc.tsv \
        --stage1-variant-qc-dir stage1/qc \
        --stage2-variant-qc-dir stage2/qc \
        --output-dir comparison_output
"""

import argparse
import os
import sys


def read_sample_qc(filepath):
    """Read sample QC TSV file. Returns dict of sample_id -> {call_rate, lrr_sd}."""
    samples = {}
    with open(filepath) as f:
        header = f.readline().strip().split('\t')
        cr_idx = header.index('call_rate') if 'call_rate' in header else 1
        sd_idx = header.index('lrr_sd') if 'lrr_sd' in header else 2
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            sample_id = fields[0]
            try:
                cr = float(fields[cr_idx]) if fields[cr_idx] != 'NA' else None
            except (ValueError, IndexError):
                cr = None
            try:
                sd = float(fields[sd_idx]) if fields[sd_idx] != 'NA' else None
            except (ValueError, IndexError):
                sd = None
            samples[sample_id] = {'call_rate': cr, 'lrr_sd': sd}
    return samples


def read_variant_missingness(filepath):
    """Read plink2 .vmiss file. Returns list of F_MISS values."""
    values = []
    if not os.path.exists(filepath):
        return values
    with open(filepath) as f:
        header = f.readline()  # skip header
        for line in f:
            fields = line.strip().split()
            if len(fields) >= 5:
                try:
                    values.append(float(fields[-1]))
                except ValueError:
                    pass
    return values


def read_hardy(filepath):
    """Read plink2 .hardy file. Returns list of p-values."""
    values = []
    if not os.path.exists(filepath):
        return values
    with open(filepath) as f:
        header = f.readline()  # skip header
        for line in f:
            fields = line.strip().split()
            if len(fields) >= 5:
                try:
                    values.append(float(fields[-1]))
                except ValueError:
                    pass
    return values


def read_afreq(filepath):
    """Read plink2 .afreq file. Returns list of MAF values."""
    values = []
    if not os.path.exists(filepath):
        return values
    with open(filepath) as f:
        header = f.readline()  # skip header
        for line in f:
            fields = line.strip().split()
            if len(fields) >= 5:
                try:
                    af = float(fields[4])
                    maf = min(af, 1.0 - af)
                    values.append(maf)
                except (ValueError, IndexError):
                    pass
    return values


def compute_stats(values):
    """Compute summary statistics for a list of values."""
    if not values:
        return None
    n = len(values)
    mean = sum(values) / n
    sorted_v = sorted(values)
    median = sorted_v[n // 2] if n % 2 == 1 else (sorted_v[n // 2 - 1] + sorted_v[n // 2]) / 2
    variance = sum((v - mean) ** 2 for v in values) / (n - 1) if n > 1 else 0
    return {
        'n': n,
        'mean': mean,
        'median': median,
        'min': min(values),
        'max': max(values),
        'sd': variance ** 0.5,
    }


def generate_comparison_report(s1_data, s2_data, s1_vqc_dir, s2_vqc_dir):
    """Generate a text comparison report."""
    lines = []
    lines.append("=" * 66)
    lines.append("  QC Comparison Report: Pre vs Post Reclustering")
    lines.append("=" * 66)
    lines.append("")

    # Find matched samples
    common = sorted(set(s1_data.keys()) & set(s2_data.keys()))
    lines.append(f"Matched samples: {len(common)}")
    lines.append(f"Stage 1 only:    {len(set(s1_data.keys()) - set(s2_data.keys()))}")
    lines.append(f"Stage 2 only:    {len(set(s2_data.keys()) - set(s1_data.keys()))}")
    lines.append("")

    # Call rate comparison
    s1_cr = [s1_data[s]['call_rate'] for s in common if s1_data[s]['call_rate'] is not None]
    s2_cr = [s2_data[s]['call_rate'] for s in common if s2_data[s]['call_rate'] is not None]
    matched_cr = [(s1_data[s]['call_rate'], s2_data[s]['call_rate']) for s in common
                  if s1_data[s]['call_rate'] is not None and s2_data[s]['call_rate'] is not None]

    lines.append("-" * 66)
    lines.append("  Sample Call Rate (autosomes)")
    lines.append("-" * 66)
    for label, vals in [("Stage 1", s1_cr), ("Stage 2", s2_cr)]:
        stats = compute_stats(vals)
        if stats:
            lines.append(f"  {label}:  mean={stats['mean']:.4f}  median={stats['median']:.4f}"
                         f"  min={stats['min']:.4f}  max={stats['max']:.4f}  n={stats['n']}")
    if matched_cr:
        diffs = [s2 - s1 for s1, s2 in matched_cr]
        diff_stats = compute_stats(diffs)
        improved = sum(1 for d in diffs if d > 0)
        lines.append(f"  Change:  mean={diff_stats['mean']:+.4f}  median={diff_stats['median']:+.4f}")
        lines.append(f"  Improved: {improved}/{len(diffs)} ({100 * improved / len(diffs):.1f}%) samples")
        for threshold in [0.95, 0.97, 0.98, 0.99]:
            n1 = sum(1 for cr in s1_cr if cr >= threshold)
            n2 = sum(1 for cr in s2_cr if cr >= threshold)
            lines.append(f"  Passing CR>={threshold}: {n1} -> {n2} ({n2 - n1:+d})")
    lines.append("")

    # LRR SD comparison
    s1_sd = [s1_data[s]['lrr_sd'] for s in common if s1_data[s]['lrr_sd'] is not None]
    s2_sd = [s2_data[s]['lrr_sd'] for s in common if s2_data[s]['lrr_sd'] is not None]
    matched_sd = [(s1_data[s]['lrr_sd'], s2_data[s]['lrr_sd']) for s in common
                  if s1_data[s]['lrr_sd'] is not None and s2_data[s]['lrr_sd'] is not None]

    lines.append("-" * 66)
    lines.append("  Sample LRR SD (autosomes)")
    lines.append("-" * 66)
    for label, vals in [("Stage 1", s1_sd), ("Stage 2", s2_sd)]:
        stats = compute_stats(vals)
        if stats:
            lines.append(f"  {label}:  mean={stats['mean']:.4f}  median={stats['median']:.4f}"
                         f"  min={stats['min']:.4f}  max={stats['max']:.4f}  n={stats['n']}")
    if matched_sd:
        diffs = [s2 - s1 for s1, s2 in matched_sd]
        diff_stats = compute_stats(diffs)
        improved = sum(1 for d in diffs if d < 0)
        lines.append(f"  Change:  mean={diff_stats['mean']:+.4f}  median={diff_stats['median']:+.4f}")
        lines.append(f"  Improved: {improved}/{len(diffs)} ({100 * improved / len(diffs):.1f}%) samples")
        for threshold in [0.20, 0.25, 0.30, 0.35]:
            n1 = sum(1 for sd in s1_sd if sd <= threshold)
            n2 = sum(1 for sd in s2_sd if sd <= threshold)
            lines.append(f"  Passing LRR_SD<={threshold}: {n1} -> {n2} ({n2 - n1:+d})")
    lines.append("")

    # Variant QC comparison
    for qc_name, reader, thresholds, direction in [
        ("Variant Missingness", read_variant_missingness, [0.01, 0.02, 0.05], "below"),
        ("HWE p-value", read_hardy, [1e-6, 1e-10, 1e-20], "above"),
        ("Minor Allele Frequency", read_afreq, [0.01, 0.05], "above"),
    ]:
        # Map to actual plink2 file extensions
        ext_map = {"Variant Missingness": "vmiss", "HWE p-value": "hardy",
                    "Minor Allele Frequency": "afreq"}
        ext = ext_map[qc_name]

        s1_vals = reader(os.path.join(s1_vqc_dir, f"variant_qc.{ext}")) if s1_vqc_dir else []
        s2_vals = reader(os.path.join(s2_vqc_dir, f"variant_qc.{ext}")) if s2_vqc_dir else []

        if s1_vals or s2_vals:
            lines.append("-" * 66)
            lines.append(f"  {qc_name}")
            lines.append("-" * 66)
            for label, vals in [("Stage 1", s1_vals), ("Stage 2", s2_vals)]:
                stats = compute_stats(vals)
                if stats:
                    lines.append(f"  {label}: n={stats['n']}  mean={stats['mean']:.6f}"
                                 f"  median={stats['median']:.6f}")
            for t in thresholds:
                n1 = sum(1 for v in s1_vals if (v <= t if direction == "below" else v >= t)) if s1_vals else 0
                n2 = sum(1 for v in s2_vals if (v <= t if direction == "below" else v >= t)) if s2_vals else 0
                t_str = f"{t}" if t >= 0.01 else f"{t:.0e}"
                op = "<=" if direction == "below" else ">="
                l1 = f"{n1}" if s1_vals else "N/A"
                l2 = f"{n2}" if s2_vals else "N/A"
                diff_str = f" ({n2 - n1:+d})" if s1_vals and s2_vals else ""
                lines.append(f"  Variants {op}{t_str}: {l1} -> {l2}{diff_str}")
            lines.append("")

    lines.append("=" * 66)
    return '\n'.join(lines)


def create_plots(s1_data, s2_data, output_dir):
    """Create diagnostic comparison plots."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("Warning: matplotlib not available. Skipping plots.", file=sys.stderr)
        return

    common = sorted(set(s1_data.keys()) & set(s2_data.keys()))

    # Gather matched data
    cr_s1, cr_s2 = [], []
    sd_s1, sd_s2 = [], []
    for s in common:
        if s1_data[s]['call_rate'] is not None and s2_data[s]['call_rate'] is not None:
            cr_s1.append(s1_data[s]['call_rate'])
            cr_s2.append(s2_data[s]['call_rate'])
        if s1_data[s]['lrr_sd'] is not None and s2_data[s]['lrr_sd'] is not None:
            sd_s1.append(s1_data[s]['lrr_sd'])
            sd_s2.append(s2_data[s]['lrr_sd'])

    if not cr_s1 and not sd_s1:
        print("Warning: No matched samples with valid metrics. Skipping plots.", file=sys.stderr)
        return

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('QC Comparison: Pre vs Post Reclustering', fontsize=14, fontweight='bold')

    # Row 1: Call Rate
    if cr_s1:
        # Histogram overlay
        ax = axes[0, 0]
        bins = 50
        ax.hist(cr_s1, bins=bins, alpha=0.6, label='Stage 1 (baseline)', color='steelblue')
        ax.hist(cr_s2, bins=bins, alpha=0.6, label='Stage 2 (reclustered)', color='darkorange')
        ax.set_xlabel('Call Rate')
        ax.set_ylabel('Count')
        ax.set_title('Call Rate Distribution')
        ax.legend(fontsize=8)
        ax.axvline(x=0.97, color='red', linestyle='--', alpha=0.5, label='Threshold')

        # Scatter: stage1 vs stage2
        ax = axes[0, 1]
        ax.scatter(cr_s1, cr_s2, alpha=0.3, s=5, color='steelblue')
        lims = [min(min(cr_s1), min(cr_s2)), max(max(cr_s1), max(cr_s2))]
        ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=1)
        ax.set_xlabel('Stage 1 Call Rate')
        ax.set_ylabel('Stage 2 Call Rate')
        ax.set_title('Call Rate: Stage 1 vs Stage 2')

        # Change distribution
        ax = axes[0, 2]
        cr_diffs = [s2 - s1 for s1, s2 in zip(cr_s1, cr_s2)]
        ax.hist(cr_diffs, bins=50, color='steelblue', alpha=0.7)
        ax.axvline(x=0, color='red', linestyle='--', alpha=0.5)
        mean_diff = sum(cr_diffs) / len(cr_diffs) if cr_diffs else 0
        ax.axvline(x=mean_diff, color='green', linestyle='-', alpha=0.7,
                   label=f'Mean: {mean_diff:+.4f}')
        ax.set_xlabel('Call Rate Change (Stage2 - Stage1)')
        ax.set_ylabel('Count')
        ax.set_title('Call Rate Change per Sample')
        ax.legend(fontsize=8)
    else:
        for j in range(3):
            axes[0, j].text(0.5, 0.5, 'No call rate data', ha='center', va='center',
                            transform=axes[0, j].transAxes)

    # Row 2: LRR SD
    if sd_s1:
        # Histogram overlay
        ax = axes[1, 0]
        ax.hist(sd_s1, bins=50, alpha=0.6, label='Stage 1 (baseline)', color='steelblue')
        ax.hist(sd_s2, bins=50, alpha=0.6, label='Stage 2 (reclustered)', color='darkorange')
        ax.set_xlabel('LRR SD')
        ax.set_ylabel('Count')
        ax.set_title('LRR SD Distribution')
        ax.legend(fontsize=8)
        ax.axvline(x=0.35, color='red', linestyle='--', alpha=0.5, label='Threshold')

        # Scatter: stage1 vs stage2
        ax = axes[1, 1]
        ax.scatter(sd_s1, sd_s2, alpha=0.3, s=5, color='steelblue')
        lims = [min(min(sd_s1), min(sd_s2)), max(max(sd_s1), max(sd_s2))]
        ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=1)
        ax.set_xlabel('Stage 1 LRR SD')
        ax.set_ylabel('Stage 2 LRR SD')
        ax.set_title('LRR SD: Stage 1 vs Stage 2')

        # Change distribution
        ax = axes[1, 2]
        sd_diffs = [s2 - s1 for s1, s2 in zip(sd_s1, sd_s2)]
        ax.hist(sd_diffs, bins=50, color='steelblue', alpha=0.7)
        ax.axvline(x=0, color='red', linestyle='--', alpha=0.5)
        mean_diff = sum(sd_diffs) / len(sd_diffs) if sd_diffs else 0
        ax.axvline(x=mean_diff, color='green', linestyle='-', alpha=0.7,
                   label=f'Mean: {mean_diff:+.4f}')
        ax.set_xlabel('LRR SD Change (Stage2 - Stage1)')
        ax.set_ylabel('Count')
        ax.set_title('LRR SD Change per Sample')
        ax.legend(fontsize=8)
    else:
        for j in range(3):
            axes[1, j].text(0.5, 0.5, 'No LRR SD data', ha='center', va='center',
                            transform=axes[1, j].transAxes)

    plt.tight_layout()
    plot_path = os.path.join(output_dir, 'qc_comparison_samples.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Sample QC comparison plot: {plot_path}")


def create_variant_plots(s1_vqc_dir, s2_vqc_dir, output_dir):
    """Create variant-level QC comparison plots."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return

    s1_miss = read_variant_missingness(os.path.join(s1_vqc_dir, "variant_qc.vmiss")) if s1_vqc_dir else []
    s2_miss = read_variant_missingness(os.path.join(s2_vqc_dir, "variant_qc.vmiss")) if s2_vqc_dir else []
    s1_hwe = read_hardy(os.path.join(s1_vqc_dir, "variant_qc.hardy")) if s1_vqc_dir else []
    s2_hwe = read_hardy(os.path.join(s2_vqc_dir, "variant_qc.hardy")) if s2_vqc_dir else []
    s1_maf = read_afreq(os.path.join(s1_vqc_dir, "variant_qc.afreq")) if s1_vqc_dir else []
    s2_maf = read_afreq(os.path.join(s2_vqc_dir, "variant_qc.afreq")) if s2_vqc_dir else []

    has_data = any([s1_miss, s2_miss, s1_hwe, s2_hwe, s1_maf, s2_maf])
    if not has_data:
        return

    import math

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle('Variant QC Comparison: Pre vs Post Reclustering', fontsize=14, fontweight='bold')

    # Variant missingness
    ax = axes[0]
    if s1_miss or s2_miss:
        if s1_miss:
            ax.hist(s1_miss, bins=50, alpha=0.6, label=f'Stage 1 (n={len(s1_miss)})',
                    color='steelblue')
        if s2_miss:
            ax.hist(s2_miss, bins=50, alpha=0.6, label=f'Stage 2 (n={len(s2_miss)})',
                    color='darkorange')
        ax.set_xlabel('Variant Missingness (F_MISS)')
        ax.set_ylabel('Count')
        ax.set_title('Variant Missingness')
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, 'No missingness data', ha='center', va='center',
                transform=ax.transAxes)

    # HWE p-values (log scale)
    ax = axes[1]
    if s1_hwe or s2_hwe:
        if s1_hwe:
            s1_log = [-math.log10(max(p, 1e-300)) for p in s1_hwe]
            ax.hist(s1_log, bins=50, alpha=0.6, label=f'Stage 1 (n={len(s1_hwe)})',
                    color='steelblue')
        if s2_hwe:
            s2_log = [-math.log10(max(p, 1e-300)) for p in s2_hwe]
            ax.hist(s2_log, bins=50, alpha=0.6, label=f'Stage 2 (n={len(s2_hwe)})',
                    color='darkorange')
        ax.set_xlabel('-log10(HWE p-value)')
        ax.set_ylabel('Count')
        ax.set_title('Hardy-Weinberg Equilibrium')
        ax.axvline(x=6, color='red', linestyle='--', alpha=0.5, label='p=1e-6')
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, 'No HWE data', ha='center', va='center',
                transform=ax.transAxes)

    # MAF distribution
    ax = axes[2]
    if s1_maf or s2_maf:
        if s1_maf:
            ax.hist(s1_maf, bins=50, alpha=0.6, label=f'Stage 1 (n={len(s1_maf)})',
                    color='steelblue')
        if s2_maf:
            ax.hist(s2_maf, bins=50, alpha=0.6, label=f'Stage 2 (n={len(s2_maf)})',
                    color='darkorange')
        ax.set_xlabel('Minor Allele Frequency')
        ax.set_ylabel('Count')
        ax.set_title('MAF Distribution')
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, 'No MAF data', ha='center', va='center',
                transform=ax.transAxes)

    plt.tight_layout()
    plot_path = os.path.join(output_dir, 'qc_comparison_variants.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Variant QC comparison plot: {plot_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate diagnostic plots comparing QC metrics pre vs post reclustering"
    )
    parser.add_argument("--stage1-sample-qc", required=True,
                        help="Stage 1 sample QC TSV file")
    parser.add_argument("--stage2-sample-qc", required=True,
                        help="Stage 2 sample QC TSV file")
    parser.add_argument("--stage1-variant-qc-dir", default=None,
                        help="Stage 1 variant QC directory (containing plink2 output)")
    parser.add_argument("--stage2-variant-qc-dir", default=None,
                        help="Stage 2 variant QC directory (containing plink2 output)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for plots and report")
    args = parser.parse_args()

    for f in [args.stage1_sample_qc, args.stage2_sample_qc]:
        if not os.path.exists(f):
            print(f"Error: File not found: {f}", file=sys.stderr)
            sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    print("Generating QC comparison report and plots...")

    # Read sample QC data
    s1_data = read_sample_qc(args.stage1_sample_qc)
    s2_data = read_sample_qc(args.stage2_sample_qc)

    # Generate comparison report
    report = generate_comparison_report(
        s1_data, s2_data,
        args.stage1_variant_qc_dir, args.stage2_variant_qc_dir
    )
    report_path = os.path.join(args.output_dir, 'qc_comparison_report.txt')
    with open(report_path, 'w') as f:
        f.write(report)
    print(report)
    print(f"\n  Report saved: {report_path}")

    # Generate plots
    create_plots(s1_data, s2_data, args.output_dir)

    # Generate variant QC plots if data available
    if args.stage1_variant_qc_dir or args.stage2_variant_qc_dir:
        create_variant_plots(
            args.stage1_variant_qc_dir, args.stage2_variant_qc_dir,
            args.output_dir
        )

    print("\n  QC comparison complete.")


if __name__ == "__main__":
    main()
