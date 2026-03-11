#!/usr/bin/env python3
"""
generate_mock_data.py

Generate deterministic mock pipeline output for integration testing and
GitHub Pages demo report.  All random values use a fixed seed so output
is reproducible across runs.

Usage:
    python3 tests/generate_mock_data.py --output-dir /tmp/mock_output
"""

import argparse
import math
import os
import random
import sys

# Fixed seed for reproducibility
SEED = 42


def _rng():
    """Return a seeded RNG."""
    return random.Random(SEED)


def generate_sample_qc(output_dir, stage, n_samples=30):
    """Write a stage sample_qc.tsv with deterministic mock data.

    The mock cohort is designed to exercise every QC feature:
      - Most samples pass all thresholds
      - 2-3 samples fail call rate (< 0.97)
      - 2-3 samples fail LRR SD (> 0.35)
      - 1-2 samples flagged for BAF SD (> 0.15)
      - 1-2 samples are het rate outliers (very high/low het)
      - Mix of M/F predicted sex
    """
    rng = _rng()

    qc_dir = os.path.join(output_dir, stage, 'qc')
    os.makedirs(qc_dir, exist_ok=True)

    rows = []
    for i in range(n_samples):
        sid = f"SAMPLE_{i + 1:03d}"

        # --- Call Rate ---
        if i == 3:
            cr = 0.9521  # fail: below 0.97
        elif i == 12:
            cr = 0.9388  # fail: below 0.97
        elif i == 22:
            cr = 0.9655  # fail: below 0.97
        else:
            cr = round(rng.uniform(0.9700, 0.9999), 4)

        # --- LRR SD ---
        if i == 5:
            lsd = 0.4123  # fail: above 0.35
        elif i == 14:
            lsd = 0.3891  # fail: above 0.35
        elif i == 24:
            lsd = 0.3601  # fail: above 0.35
        else:
            lsd = round(rng.uniform(0.08, 0.30), 4)

        # --- LRR Mean / Median ---
        lrr_mean = round(rng.uniform(-0.02, 0.02), 4)
        lrr_median = round(lrr_mean + rng.uniform(-0.005, 0.005), 4)

        # --- BAF SD (het sites) ---
        if i == 7:
            bsd = 0.1823  # flag: above 0.15
        elif i == 18:
            bsd = 0.1645  # flag: above 0.15
        else:
            bsd = round(rng.uniform(0.02, 0.09), 4)

        # --- Het Rate ---
        if i == 9:
            hr = 0.3412  # very high → outlier
        elif i == 20:
            hr = 0.1534  # very low → outlier
        else:
            hr = round(rng.uniform(0.22, 0.28), 4)

        # --- Gender ---
        gender = 'M' if i % 2 == 0 else 'F'

        rows.append({
            'sample_id': sid,
            'call_rate': f"{cr:.4f}",
            'lrr_sd': f"{lsd:.4f}",
            'lrr_mean': f"{lrr_mean:.4f}",
            'lrr_median': f"{lrr_median:.4f}",
            'baf_sd': f"{bsd:.4f}",
            'het_rate': f"{hr:.4f}",
            'computed_gender': gender,
        })

    filepath = os.path.join(qc_dir, f'{stage}_sample_qc.tsv')
    cols = ['sample_id', 'call_rate', 'lrr_sd', 'lrr_mean',
            'lrr_median', 'baf_sd', 'het_rate', 'computed_gender']
    with open(filepath, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(r[c] for c in cols) + '\n')

    return filepath


def generate_variant_qc_summary(output_dir, stage):
    """Write a mock variant QC summary."""
    vqc_dir = os.path.join(output_dir, stage, 'qc', 'variant_qc')
    os.makedirs(vqc_dir, exist_ok=True)

    text = """==================================================
  Variant QC Summary (autosomes)
==================================================

  Input variants:         650,427
  After call-rate filter: 642,185 (≥ 0.98)
  After HWE filter:       639,412 (p ≥ 1e-6)
  After MAF filter:       531,208 (≥ 0.01)

  Ti/Tv ratio:            2.07
  Mean inbreeding F:      0.0012
  Inbreeding F range:     -0.032 to 0.041

  Missingness histogram:
    0.00-0.01:  ████████████████████████  598,127
    0.01-0.02:  ████                       44,058
    0.02-0.05:  ██                          6,923
    0.05-0.10:  ▏                           1,062
    0.10-1.00:  ▏                             257

  MAF distribution:
    0.00-0.01:  ████████████████           108,977
    0.01-0.05:  ████████████████████       132,245
    0.05-0.10:  ██████████████             98,456
    0.10-0.25:  ████████████████████████  164,123
    0.25-0.50:  ████████████████████████  146,626
"""
    filepath = os.path.join(vqc_dir, 'variant_qc_summary.txt')
    with open(filepath, 'w') as f:
        f.write(text)
    return filepath


def generate_comparison_report(output_dir):
    """Write a mock QC comparison report for stage1→stage2."""
    comp_dir = os.path.join(output_dir, 'stage2', 'qc', 'comparison')
    os.makedirs(comp_dir, exist_ok=True)

    text = """==================================================
  Stage 1 → Stage 2 QC Comparison
==================================================

  Metric              Stage 1     Stage 2     Change
  ─────────────────── ─────────── ─────────── ──────────
  Mean Call Rate       0.9841      0.9893      +0.0052
  Mean LRR SD          0.2312      0.2045      -0.0267
  Mean BAF SD          0.0612      0.0534      -0.0078
  Mean Het Rate        0.2489      0.2478      -0.0011

  Variant Missingness  0.0134      0.0098      -0.0036
  HWE failures (1e-6)  1,245       823         -422
  Mean MAF             0.1723      0.1734      +0.0011
  Ti/Tv ratio          2.05        2.07        +0.02

  Reclustering improved call rate for 28/30 samples
  and reduced LRR SD for 26/30 samples.
"""
    filepath = os.path.join(comp_dir, 'qc_comparison_report.txt')
    with open(filepath, 'w') as f:
        f.write(text)
    return filepath


def generate_realign_summary(output_dir):
    """Write a mock manifest realignment summary."""
    realign_dir = os.path.join(output_dir, 'realigned_manifests')
    os.makedirs(realign_dir, exist_ok=True)

    text = """==================================================
  Manifest Realignment Summary
==================================================

  Manifest:         GSA-24v3-0_A1
  Reference:        chm13v2.0 (T2T-CHM13v2.0)
  Aligner:          bwa 0.7.18

  Total probes:     654,027
  Mapped:           650,427 (99.45%)
  Unmapped:         2,134 (0.33%)
  Multi-mapped:     1,466 (0.22%)

  Position changes: 3,812 (0.58%)
    Same chromosome, different position:  3,701
    Different chromosome:                    111

  Strand changes:   245
  New placements:   89 (previously unmapped)
"""
    filepath = os.path.join(realign_dir, 'realign_summary.txt')
    with open(filepath, 'w') as f:
        f.write(text)
    return filepath


def generate_diagnostic_report(output_dir):
    """Write a mock QC diagnostic report."""
    text = """==================================================
  QC Diagnostic Report
==================================================

  ✓ Call rate distribution appears normal
  ✓ LRR SD distribution is within expected range
  ⚠ 2 samples flagged for BAF SD > 0.15 (possible contamination)
  ⚠ 2 samples are heterozygosity rate outliers (± 3 SD)
  ✓ Gender predictions consistent (15 M, 15 F)
  ✓ Ti/Tv ratio 2.07 is within expected range (1.8-2.2)
  ✓ No build mismatch detected between manifest and reference

  Recommendations:
  - Review BAF SD-flagged samples (SAMPLE_008, SAMPLE_019) for
    possible contamination or mosaicism
  - Review het rate outliers (SAMPLE_010, SAMPLE_021) for
    population stratification or sample swaps
"""
    filepath = os.path.join(output_dir, 'qc_diagnostic_report.txt')
    with open(filepath, 'w') as f:
        f.write(text)
    return filepath


def generate_pca_data(output_dir, n_samples=30):
    """Write mock PCA projections and QC summary."""
    pca_dir = os.path.join(output_dir, 'ancestry_pca')
    os.makedirs(pca_dir, exist_ok=True)

    rng = _rng()

    # PCA projections (20 PCs)
    header = ['FID', 'IID'] + [f'PC{i}' for i in range(1, 21)]
    filepath = os.path.join(pca_dir, 'pca_projections.tsv')
    with open(filepath, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for i in range(n_samples):
            sid = f"SAMPLE_{i + 1:03d}"
            # Create 3 rough clusters
            cluster = i % 3
            base_pc1 = [-0.02, 0.01, 0.03][cluster]
            base_pc2 = [0.01, -0.02, 0.02][cluster]
            pcs = [round(base_pc1 + rng.gauss(0, 0.005), 6),
                   round(base_pc2 + rng.gauss(0, 0.005), 6)]
            for _ in range(18):
                pcs.append(round(rng.gauss(0, 0.003), 6))
            row = [sid, sid] + [str(p) for p in pcs]
            f.write('\t'.join(row) + '\n')

    # QC summary
    qc_text = """==================================================
  PCA QC Summary
==================================================

  Input samples:          30
  After call rate ≥ 0.98: 28
  After mind ≤ 0.02:      28

  Input variants:         650,427
  After geno ≤ 0.02:      642,185
  After HWE (1e-6):       639,412
  After MAF ≥ 0.05:       485,234
  After LD pruning:       98,412

  PCs computed: 20
  Eigenvalues: 12.34, 8.56, 5.23, 3.89, 2.67, ...
"""
    with open(os.path.join(pca_dir, 'qc_summary.txt'), 'w') as f:
        f.write(qc_text)

    return filepath


def generate_sex_check_data(output_dir, n_samples=30):
    """Write mock sex-check chrX/chrY median LRR table for interactive plotting."""
    sex_dir = os.path.join(output_dir, 'sex_check')
    os.makedirs(sex_dir, exist_ok=True)
    rng = _rng()

    filepath = os.path.join(sex_dir, 'sex_check_chrXY_lrr.tsv')
    with open(filepath, 'w') as f:
        f.write('sample_id\tchrx_lrr_median\tchry_lrr_median\tcomputed_gender\n')
        for i in range(n_samples):
            sid = f"SAMPLE_{i + 1:03d}"
            male = (i % 2 == 0)
            if male:
                chrx = -0.11 + rng.gauss(0, 0.012)
                chry = 0.08 + rng.gauss(0, 0.012)
                sex = 'M'
            else:
                chrx = 0.03 + rng.gauss(0, 0.012)
                chry = -0.07 + rng.gauss(0, 0.012)
                sex = 'F'
            f.write(f'{sid}\t{chrx:.4f}\t{chry:.4f}\t{sex}\n')

    return filepath


def generate_mock_notice(output_dir):
    """Write a marker text indicating this output is synthetic example data."""
    filepath = os.path.join(output_dir, 'mock_data_notice.txt')
    with open(filepath, 'w') as f:
        f.write("This report is an example generated from deterministic mock data for demonstration purposes only.")
    return filepath


def generate_peddy_data(output_dir, n_samples=30):
    """Write mock peddy output files (het_check.csv, sex_check.csv, ped_check.csv, peddy.ped).

    Simulates peddy's ancestry prediction, sex check, and relatedness analysis.
    """
    peddy_dir = os.path.join(output_dir, 'peddy')
    os.makedirs(peddy_dir, exist_ok=True)
    rng = _rng()

    ancestry_groups = ['EUR', 'AFR', 'EAS', 'AMR', 'SAS']

    # --- het_check.csv ---
    het_path = os.path.join(peddy_dir, 'peddy.het_check.csv')
    with open(het_path, 'w') as f:
        f.write('sample_id,sampled_sites,mean_depth,median_depth,depth_outlier,'
                'het_count,het_ratio,ratio_outlier,idr_baf,p10,p90,'
                'PC1,PC2,PC3,PC4,ancestry-prediction,ancestry-prob\n')
        for i in range(n_samples):
            sid = f"SAMPLE_{i + 1:03d}"
            cluster = i % 3
            ancestry = ancestry_groups[cluster % len(ancestry_groups)]
            prob = round(rng.uniform(0.85, 0.99), 4)
            sites = rng.randint(15000, 20000)
            depth = round(rng.uniform(25.0, 45.0), 1)
            het_count = rng.randint(3000, 5500)
            het_ratio = round(het_count / sites, 4)
            ratio_outlier = 'False'
            if i == 9:
                ratio_outlier = 'True'
            pc1 = round([-0.015, 0.008, 0.025][cluster] + rng.gauss(0, 0.004), 6)
            pc2 = round([0.008, -0.015, 0.018][cluster] + rng.gauss(0, 0.004), 6)
            pc3 = round(rng.gauss(0, 0.003), 6)
            pc4 = round(rng.gauss(0, 0.003), 6)
            idr = round(rng.uniform(0.01, 0.04), 4)
            p10 = round(rng.uniform(0.0, 0.05), 4)
            p90 = round(p10 + idr + rng.uniform(0.85, 0.95), 4)
            f.write(f'{sid},{sites},{depth},{depth},False,'
                    f'{het_count},{het_ratio},{ratio_outlier},{idr},{p10},{p90},'
                    f'{pc1},{pc2},{pc3},{pc4},{ancestry},{prob}\n')

    # --- sex_check.csv ---
    sex_path = os.path.join(peddy_dir, 'peddy.sex_check.csv')
    with open(sex_path, 'w') as f:
        f.write('sample_id,error,het_count,hom_alt_count,hom_ref_count,'
                'het_ratio,ped_sex,predicted_sex\n')
        for i in range(n_samples):
            sid = f"SAMPLE_{i + 1:03d}"
            male = (i % 2 == 0)
            predicted = 'male' if male else 'female'
            ped_sex = 'male' if male else 'female'
            error = 'False'
            if male:
                het_c = rng.randint(5, 30)
                hom_alt = rng.randint(200, 400)
                hom_ref = rng.randint(4000, 6000)
            else:
                het_c = rng.randint(1500, 3000)
                hom_alt = rng.randint(1000, 2000)
                hom_ref = rng.randint(2000, 4000)
            hr = round(het_c / (het_c + hom_alt + hom_ref), 4) if (het_c + hom_alt + hom_ref) > 0 else 0
            f.write(f'{sid},{error},{het_c},{hom_alt},{hom_ref},'
                    f'{hr},{ped_sex},{predicted}\n')

    # --- ped_check.csv ---
    ped_path = os.path.join(peddy_dir, 'peddy.ped_check.csv')
    with open(ped_path, 'w') as f:
        f.write('sample_a,sample_b,n,rel,pedigree_relatedness,rel_difference,'
                'ibs0,ibs2,shared_hets,hets_a,hets_b,pedigree_parents\n')
        # Write a few representative pairs
        for i in range(min(5, n_samples)):
            for j in range(i + 1, min(6, n_samples)):
                sid_a = f"SAMPLE_{i + 1:03d}"
                sid_b = f"SAMPLE_{j + 1:03d}"
                rel = round(rng.uniform(-0.05, 0.05), 4)
                ibs0 = rng.randint(1000, 3000)
                ibs2 = rng.randint(8000, 12000)
                shared = rng.randint(2000, 4000)
                ha = rng.randint(3000, 5000)
                hb = rng.randint(3000, 5000)
                f.write(f'{sid_a},{sid_b},{rng.randint(15000, 20000)},{rel},'
                        f'0.0,{rel},  {ibs0},{ibs2},{shared},{ha},{hb},False\n')

    # --- peddy.ped (enhanced pedigree) ---
    ped_path = os.path.join(peddy_dir, 'peddy.peddy.ped')
    with open(ped_path, 'w') as f:
        for i in range(n_samples):
            sid = f"SAMPLE_{i + 1:03d}"
            sex = '1' if i % 2 == 0 else '2'
            f.write(f'{sid}\t{sid}\t0\t0\t{sex}\t-9\n')

    # --- peddy_final.ped ---
    final_path = os.path.join(peddy_dir, 'peddy_final.ped')
    with open(final_path, 'w') as f:
        for i in range(n_samples):
            sid = f"SAMPLE_{i + 1:03d}"
            sex = '1' if i % 2 == 0 else '2'
            f.write(f'{sid}\t{sid}\t0\t0\t{sex}\t-9\n')

    # --- background_pca.json (ancestry PCA projection) ---
    bg_pca_path = os.path.join(peddy_dir, 'peddy.background_pca.json')
    import json
    bg_data = []
    for i in range(n_samples):
        sid = f"SAMPLE_{i + 1:03d}"
        cluster = i % 3
        bg_data.append({
            'sample_id': sid,
            'PC1': round([-0.015, 0.008, 0.025][cluster] + rng.gauss(0, 0.004), 6),
            'PC2': round([0.008, -0.015, 0.018][cluster] + rng.gauss(0, 0.004), 6),
            'PC3': round(rng.gauss(0, 0.003), 6),
            'PC4': round(rng.gauss(0, 0.003), 6),
            'ancestry': ancestry_groups[cluster % len(ancestry_groups)],
        })
    with open(bg_pca_path, 'w') as f:
        json.dump(bg_data, f)

    return het_path


def generate_ancestry_stratified_qc(output_dir, n_samples=30):
    """Write mock ancestry-stratified QC outputs.

    Creates per-ancestry variant QC summaries and a collated variant QC
    file for the ancestry groups that exceed the minimum sample threshold
    in the mock data.
    """
    rng = _rng()
    strat_dir = os.path.join(output_dir, 'ancestry_stratified_qc')
    os.makedirs(strat_dir, exist_ok=True)

    ancestry_groups = ['EUR', 'AFR', 'EAS', 'AMR', 'SAS']
    # In mock data, ancestry is assigned as i % 3 -> EUR, AFR, EAS cycle
    # so with 30 samples: EUR=10, AFR=10, EAS=10, AMR=0, SAS=0
    ancestry_counts = {}
    for i in range(n_samples):
        anc = ancestry_groups[i % 3]
        ancestry_counts[anc] = ancestry_counts.get(anc, 0) + 1

    # Write ancestry assignments
    assign_path = os.path.join(strat_dir, 'ancestry_assignments.tsv')
    with open(assign_path, 'w') as f:
        f.write('sample_id\tancestry\n')
        for i in range(n_samples):
            sid = f"SAMPLE_{i + 1:03d}"
            anc = ancestry_groups[i % 3]
            f.write(f'{sid}\t{anc}\n')

    # Generate per-ancestry variant QC for groups with >= 5 samples (mock threshold)
    mock_ancestries = [a for a, c in ancestry_counts.items() if c >= 5]

    for anc in mock_ancestries:
        anc_dir = os.path.join(strat_dir, anc)
        os.makedirs(os.path.join(anc_dir, 'variant_qc'), exist_ok=True)

        # Variant QC summary
        n_vars = rng.randint(600000, 650000)
        text = f"""==================================================
  Variant QC Summary ({anc}, autosomes)
==================================================

  Input variants:         {n_vars:,}
  After call-rate filter: {int(n_vars * 0.985):,} (≥ 0.98)
  After HWE filter:       {int(n_vars * 0.978):,} (p ≥ 1e-6)
  After MAF filter:       {int(n_vars * 0.82):,} (≥ 0.01)

  Ti/Tv ratio:            {round(rng.uniform(2.0, 2.15), 2)}
  Mean inbreeding F:      {round(rng.uniform(-0.005, 0.005), 4)}
"""
        with open(os.path.join(anc_dir, 'variant_qc', 'variant_qc_summary.txt'),
                  'w') as f:
            f.write(text)

        # Write mock .vmiss file
        vmiss_path = os.path.join(anc_dir, 'variant_qc', 'variant_qc.vmiss')
        with open(vmiss_path, 'w') as f:
            f.write('#CHROM\tID\tMISSING_CT\tOBS_CT\tF_MISS\n')
            n_count = ancestry_counts.get(anc, 10)
            for v in range(min(200, n_vars)):
                miss_ct = rng.randint(0, max(1, n_count // 5))
                f_miss = round(miss_ct / n_count, 6)
                f.write(f'chr1\trs{100000 + v}\t{miss_ct}\t{n_count}\t{f_miss}\n')

        # Write mock .hardy file
        hardy_path = os.path.join(anc_dir, 'variant_qc', 'variant_qc.hardy')
        with open(hardy_path, 'w') as f:
            f.write('#CHROM\tID\tA1\tA2\tCT\tP\n')
            for v in range(min(200, n_vars)):
                p = 10 ** rng.uniform(-12, 0)
                f.write(f'chr1\trs{100000 + v}\tA\tG\t{n_count}\t{p:.6e}\n')

        # Write mock .afreq file
        afreq_path = os.path.join(anc_dir, 'variant_qc', 'variant_qc.afreq')
        with open(afreq_path, 'w') as f:
            f.write('#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n')
            for v in range(min(200, n_vars)):
                af = round(rng.uniform(0.0, 0.5), 6)
                f.write(f'chr1\trs{100000 + v}\tA\tG\t{af}\t{n_count * 2}\n')

    # Generate collated variant QC
    collated_path = os.path.join(strat_dir, 'collated_variant_qc.tsv')
    with open(collated_path, 'w') as f:
        cols = ['variant_id', 'all_call_rate', 'all_hwe_p', 'all_maf']
        for anc in sorted(mock_ancestries):
            cols.extend([f'{anc}_call_rate', f'{anc}_hwe_p', f'{anc}_maf'])
        cols.extend([
            'all_ancestries_call_rate_pass',
            'all_ancestries_hwe_pass',
            'all_ancestries_maf_pass',
            'all_ancestries_qc_pass',
        ])
        f.write('\t'.join(cols) + '\n')
        for v in range(200):
            row = [f'rs{100000 + v}']
            # all
            cr = round(1 - rng.uniform(0, 0.05), 6)
            hwe = f'{10 ** rng.uniform(-8, 0):.6e}'
            maf = round(rng.uniform(0, 0.5), 6)
            row.extend([str(cr), hwe, str(maf)])
            cr_pass = True
            hwe_pass = True
            maf_pass = True
            for anc in sorted(mock_ancestries):
                cr = round(1 - rng.uniform(0, 0.06), 6)
                hwe = f'{10 ** rng.uniform(-10, 0):.6e}'
                maf = round(rng.uniform(0, 0.5), 6)
                row.extend([str(cr), hwe, str(maf)])
                cr_pass = cr_pass and cr >= 0.98
                hwe_pass = hwe_pass and float(hwe) >= 1e-6
                maf_pass = maf_pass and maf >= 0.01
            row.extend([
                '1' if cr_pass else '0',
                '1' if hwe_pass else '0',
                '1' if maf_pass else '0',
                '1' if (cr_pass and hwe_pass and maf_pass) else '0',
            ])
            f.write('\t'.join(row) + '\n')

    # Summary
    summary_path = os.path.join(strat_dir, 'ancestry_stratified_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("======================================================\n")
        f.write("  Ancestry-Stratified QC Summary\n")
        f.write("======================================================\n\n")
        f.write("Ancestry Groups Analyzed:\n")
        for anc in sorted(mock_ancestries):
            f.write(f"  {anc}: {ancestry_counts[anc]} samples\n")
        f.write("\n")
        f.write(f"Collated variant QC: {collated_path}\n")

    return summary_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate deterministic mock pipeline output for testing"
    )
    parser.add_argument("--output-dir", required=True,
                        help="Directory to write mock output")
    parser.add_argument("--n-samples", type=int, default=30,
                        help="Number of mock samples (default: 30)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    n = args.n_samples

    print(f"Generating mock pipeline output in: {args.output_dir}")
    print(f"  Samples: {n}")

    # Stage 1 QC (same structure, deterministic from fixed seed)
    s1 = generate_sample_qc(args.output_dir, 'stage1', n)
    print(f"  Stage 1 sample QC: {s1}")

    # Stage 2 QC (same samples, slightly better metrics)
    s2 = generate_sample_qc(args.output_dir, 'stage2', n)
    print(f"  Stage 2 sample QC: {s2}")

    # Variant QC
    vqc = generate_variant_qc_summary(args.output_dir, 'stage2')
    print(f"  Variant QC:        {vqc}")

    # Comparison report
    comp = generate_comparison_report(args.output_dir)
    print(f"  Comparison report: {comp}")

    # Realignment summary
    realign = generate_realign_summary(args.output_dir)
    print(f"  Realignment:       {realign}")

    # Diagnostic report
    diag = generate_diagnostic_report(args.output_dir)
    print(f"  Diagnostics:       {diag}")

    # PCA data
    pca = generate_pca_data(args.output_dir, n)
    print(f"  PCA projections:   {pca}")

    sex_check = generate_sex_check_data(args.output_dir, n)
    print(f"  Sex check table:   {sex_check}")

    peddy = generate_peddy_data(args.output_dir, n)
    print(f"  Peddy data:        {peddy}")

    note = generate_mock_notice(args.output_dir)
    print(f"  Mock notice:      {note}")

    ancestry_strat = generate_ancestry_stratified_qc(args.output_dir, n)
    print(f"  Ancestry strat:   {ancestry_strat}")

    print("\nDone. Run generate_report.py --output-dir to build the HTML report.")


if __name__ == "__main__":
    main()
