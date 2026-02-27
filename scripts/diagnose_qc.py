#!/usr/bin/env python3
"""
diagnose_qc.py

Diagnostic tool for troubleshooting Illumina IDAT processing QC issues.
Analyzes pipeline output to identify root causes of deflated call rates,
inflated LRR standard deviations, and other quality problems.

Can be run standalone against pipeline output or invoked by the pipeline
for automated post-run diagnostics.

Usage:
    python3 diagnose_qc.py --output-dir /path/to/pipeline_output
    python3 diagnose_qc.py --qc-file /path/to/sample_qc.tsv
    python3 diagnose_qc.py --output-dir /path/to/output --report /path/to/report.txt
"""

import argparse
import os
import sys
import statistics
from pathlib import Path


# Expected ranges for high-quality Illumina data
EXPECTED_CALL_RATE_MIN = 0.95
EXPECTED_CALL_RATE_TYPICAL = 0.98
EXPECTED_LRR_SD_MAX = 0.40
EXPECTED_LRR_SD_TYPICAL = 0.20


def parse_qc_file(filepath):
    """Parse a sample_qc.tsv file. Returns list of dicts."""
    samples = []
    with open(filepath) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            row = {}
            for i, col in enumerate(header):
                row[col] = fields[i] if i < len(fields) else 'NA'
            samples.append(row)
    return samples


def parse_norm_log(filepath):
    """Parse a bcftools norm warnings log file. Returns count of REF mismatches."""
    if not os.path.exists(filepath):
        return None
    count = 0
    with open(filepath) as f:
        for line in f:
            if 'REF_MISMATCH' in line:
                count += 1
    return count


def parse_realign_summary(filepath):
    """Parse realign_summary.txt for mapping stats."""
    if not os.path.exists(filepath):
        return None
    stats = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if 'Total flank' in line:
                parts = line.split(':')
                if len(parts) >= 2:
                    stats['total'] = int(parts[-1].strip())
            elif 'Mapped:' in line and 'Unmapped' not in line:
                parts = line.split(':')
                if len(parts) >= 2:
                    val = parts[-1].strip().split()[0]
                    stats['mapped'] = int(val)
            elif 'Unmapped:' in line:
                parts = line.split(':')
                if len(parts) >= 2:
                    val = parts[-1].strip().split()[0]
                    stats['unmapped'] = int(val)
            elif 'Unique' in line and 'MAPQ' in line:
                parts = line.split(':')
                if len(parts) >= 2:
                    val = parts[-1].strip().split()[0]
                    stats['unique'] = int(val)
            elif 'Multi-mapped' in line:
                parts = line.split(':')
                if len(parts) >= 2:
                    val = parts[-1].strip().split()[0]
                    stats['multi_mapped'] = int(val)
    return stats if stats else None


def diagnose(output_dir=None, qc_file=None, report_file=None):
    """Run diagnostic analysis on pipeline output."""
    findings = []
    warnings = []
    errors = []

    # ---------------------------------------------------------------
    # 1. Locate and parse QC files
    # ---------------------------------------------------------------
    qc_files = []
    if qc_file:
        qc_files.append(('provided', qc_file))
    elif output_dir:
        for stage in ['stage1', 'stage2']:
            for name in ['stage1_sample_qc.tsv', 'stage2_sample_qc.tsv']:
                path = os.path.join(output_dir, stage, 'qc', name)
                if os.path.exists(path):
                    qc_files.append((f'{stage}/{name}', path))
        # Also check pipeline_output subdirectory (process_1000g.sh layout)
        for stage in ['stage1', 'stage2']:
            for name in ['stage1_sample_qc.tsv', 'stage2_sample_qc.tsv']:
                path = os.path.join(output_dir, 'pipeline_output', stage, 'qc', name)
                if os.path.exists(path):
                    qc_files.append((f'pipeline_output/{stage}/{name}', path))

    if not qc_files:
        errors.append("No QC files found. Provide --qc-file or --output-dir with pipeline output.")
        return findings, warnings, errors

    for label, path in qc_files:
        findings.append(f"\n{'='*60}")
        findings.append(f"  QC Analysis: {label}")
        findings.append(f"{'='*60}")

        samples = parse_qc_file(path)
        if not samples:
            errors.append(f"No samples found in {path}")
            continue

        n_samples = len(samples)
        findings.append(f"  Samples: {n_samples}")

        # Call rate analysis
        call_rates = []
        for s in samples:
            try:
                cr = float(s.get('call_rate', 'NA'))
                call_rates.append(cr)
            except (ValueError, TypeError):
                pass

        if call_rates:
            cr_mean = statistics.mean(call_rates)
            cr_median = statistics.median(call_rates)
            cr_min = min(call_rates)
            cr_max = max(call_rates)
            cr_sd = statistics.stdev(call_rates) if len(call_rates) > 1 else 0

            findings.append(f"\n  Call Rate Statistics:")
            findings.append(f"    Mean:   {cr_mean:.4f}")
            findings.append(f"    Median: {cr_median:.4f}")
            findings.append(f"    Min:    {cr_min:.4f}")
            findings.append(f"    Max:    {cr_max:.4f}")
            findings.append(f"    StDev:  {cr_sd:.4f}")

            n_below_95 = sum(1 for cr in call_rates if cr < 0.95)
            n_below_97 = sum(1 for cr in call_rates if cr < 0.97)
            findings.append(f"    Samples below 0.95: {n_below_95} ({n_below_95*100/n_samples:.1f}%)")
            findings.append(f"    Samples below 0.97: {n_below_97} ({n_below_97*100/n_samples:.1f}%)")

            # Diagnose call rate issues
            if cr_mean < 0.80:
                errors.append(
                    f"CRITICAL: Mean call rate {cr_mean:.4f} is extremely low (expected >0.97).\n"
                    "  Likely causes:\n"
                    "    1. bcftools norm -c x was used instead of -c ws, setting genotypes\n"
                    "       to missing for REF-mismatched probes (common with cross-build data).\n"
                    "       FIX: Use 'bcftools norm -c ws' to swap REF/ALT instead of blanking.\n"
                    "    2. Manifest/reference build mismatch without proper realignment.\n"
                    "    3. Wrong manifest files for the array type.\n"
                    "    4. Corrupted IDAT files."
                )
            elif cr_mean < EXPECTED_CALL_RATE_MIN:
                warnings.append(
                    f"Mean call rate {cr_mean:.4f} is below expected ({EXPECTED_CALL_RATE_MIN}).\n"
                    "  Check bcftools norm mode (-c ws recommended) and manifest realignment."
                )

            # Uniformly low call rates suggest a systematic issue
            if cr_sd < 0.05 and cr_mean < 0.90:
                warnings.append(
                    f"Call rates are uniformly low (mean={cr_mean:.4f}, sd={cr_sd:.4f}).\n"
                    "  This pattern indicates a systematic issue (e.g., norm -c x, wrong manifest)\n"
                    "  rather than sample-specific quality problems."
                )

        # LRR SD analysis
        lrr_sds = []
        for s in samples:
            try:
                sd = float(s.get('lrr_sd', 'NA'))
                lrr_sds.append(sd)
            except (ValueError, TypeError):
                pass

        if lrr_sds:
            sd_mean = statistics.mean(lrr_sds)
            sd_median = statistics.median(lrr_sds)
            sd_min = min(lrr_sds)
            sd_max = max(lrr_sds)

            findings.append(f"\n  LRR SD Statistics:")
            findings.append(f"    Mean:   {sd_mean:.4f}")
            findings.append(f"    Median: {sd_median:.4f}")
            findings.append(f"    Min:    {sd_min:.4f}")
            findings.append(f"    Max:    {sd_max:.4f}")

            n_above_035 = sum(1 for sd in lrr_sds if sd > 0.35)
            findings.append(f"    Samples above 0.35: {n_above_035} ({n_above_035*100/n_samples:.1f}%)")

            if sd_mean > 1.0:
                errors.append(
                    f"CRITICAL: Mean LRR SD {sd_mean:.4f} is extremely high (expected <0.30).\n"
                    "  This is likely a downstream effect of the same issue causing low call rates.\n"
                    "  When many probes have missing genotypes (e.g., from norm -c x), the\n"
                    "  remaining LRR values are biased and have inflated variance.\n"
                    "  FIX: Resolve the call rate issue first (use norm -c ws)."
                )
            elif sd_mean > EXPECTED_LRR_SD_MAX:
                warnings.append(
                    f"Mean LRR SD {sd_mean:.4f} is above expected ({EXPECTED_LRR_SD_MAX})."
                )

    # ---------------------------------------------------------------
    # 2. Check norm warnings log
    # ---------------------------------------------------------------
    if output_dir:
        norm_log_paths = []
        for root, dirs, files in os.walk(output_dir):
            for f in files:
                if 'norm_warnings' in f and f.endswith('.log'):
                    norm_log_paths.append(os.path.join(root, f))

        for log_path in norm_log_paths:
            n_swaps = parse_norm_log(log_path)
            if n_swaps is not None:
                rel_path = os.path.relpath(log_path, output_dir)
                findings.append(f"\n  Norm warnings ({rel_path}): {n_swaps} REF swaps")
                if n_swaps > 10000:
                    findings.append(
                        f"    Note: {n_swaps} REF/ALT swaps is expected for cross-build data\n"
                        f"    (e.g., GRCh37 manifest on GRCh38 reference). These swaps are\n"
                        f"    handled correctly by -c ws mode — genotypes are preserved."
                    )

    # ---------------------------------------------------------------
    # 3. Check realignment summary
    # ---------------------------------------------------------------
    if output_dir:
        for search_dir in [output_dir, os.path.join(output_dir, 'pipeline_output')]:
            summary_path = os.path.join(search_dir, 'realigned_manifests', 'realign_summary.txt')
            if os.path.exists(summary_path):
                stats = parse_realign_summary(summary_path)
                if stats:
                    findings.append(f"\n  Manifest Realignment Summary:")
                    total = stats.get('total', 0)
                    mapped = stats.get('mapped', 0)
                    unmapped = stats.get('unmapped', 0)
                    unique = stats.get('unique', 0)
                    multi = stats.get('multi_mapped', 0)
                    findings.append(f"    Total probes:    {total}")
                    findings.append(f"    Mapped:          {mapped} ({mapped*100/total:.1f}%)" if total else "")
                    findings.append(f"    Unmapped:        {unmapped} ({unmapped*100/total:.1f}%)" if total else "")
                    findings.append(f"    Unique (MQ>=10): {unique}" if unique else "")
                    findings.append(f"    Multi-mapped:    {multi}" if multi else "")

                    if total > 0 and mapped * 100 / total < 90:
                        warnings.append(
                            f"Only {mapped*100/total:.1f}% of probes mapped during realignment.\n"
                            "  This may indicate a significant build mismatch or manifest issue."
                        )

    # ---------------------------------------------------------------
    # 4. Check for known issues
    # ---------------------------------------------------------------
    if output_dir:
        # Check if old -c x flag is present in any script or log
        for root, dirs, files in os.walk(output_dir):
            for f in files:
                if f.endswith('.log') or f.endswith('.sh'):
                    path = os.path.join(root, f)
                    try:
                        with open(path) as fh:
                            content = fh.read()
                            if 'norm' in content and '-c x' in content:
                                warnings.append(
                                    f"Found '-c x' in {os.path.relpath(path, output_dir)}.\n"
                                    "  This mode sets genotypes to missing on REF mismatch.\n"
                                    "  Use '-c ws' instead to swap REF/ALT and preserve genotypes."
                                )
                    except (IOError, UnicodeDecodeError):
                        pass

    return findings, warnings, errors


def format_report(findings, warnings, errors):
    """Format diagnostic results into a readable report."""
    lines = []
    lines.append("=" * 60)
    lines.append("  Illumina IDAT Pipeline Diagnostic Report")
    lines.append("=" * 60)

    for line in findings:
        if line:
            lines.append(line)

    if warnings:
        lines.append(f"\n{'='*60}")
        lines.append(f"  WARNINGS ({len(warnings)})")
        lines.append(f"{'='*60}")
        for i, w in enumerate(warnings, 1):
            lines.append(f"\n  [{i}] {w}")

    if errors:
        lines.append(f"\n{'='*60}")
        lines.append(f"  ERRORS ({len(errors)})")
        lines.append(f"{'='*60}")
        for i, e in enumerate(errors, 1):
            lines.append(f"\n  [{i}] {e}")

    if not warnings and not errors:
        lines.append(f"\n{'='*60}")
        lines.append("  No issues detected. QC metrics appear normal.")
        lines.append(f"{'='*60}")

    lines.append("")
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Diagnose QC issues in Illumina IDAT processing pipeline output"
    )
    parser.add_argument("--output-dir", help="Pipeline output directory")
    parser.add_argument("--qc-file", help="Path to a sample_qc.tsv file")
    parser.add_argument("--report", help="Write diagnostic report to file")
    args = parser.parse_args()

    if not args.output_dir and not args.qc_file:
        parser.error("Provide --output-dir or --qc-file")

    findings, warnings, errors = diagnose(
        output_dir=args.output_dir,
        qc_file=args.qc_file
    )

    report = format_report(findings, warnings, errors)
    print(report)

    if args.report:
        with open(args.report, 'w') as f:
            f.write(report)
        print(f"Report saved to: {args.report}")

    # Exit with non-zero status if errors found
    if errors:
        sys.exit(1)


if __name__ == "__main__":
    main()

