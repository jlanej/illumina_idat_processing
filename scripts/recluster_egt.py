#!/usr/bin/env python3
"""
recluster_egt.py

Recompute Illumina EGT cluster file from high-quality samples' GTC data.

Reads the original EGT as a template, extracts normalized intensity data
(theta, R) from GTC files of high-quality samples, recomputes cluster
statistics (mean/stdev of theta and R for each genotype class AA/AB/BB)
per probe, and writes a new EGT cluster file.

The new EGT can then be used with idat2gtc to re-call genotypes for all
samples, producing more accurate genotype calls tailored to the study
population rather than relying on Illumina's default training set.

Usage:
    python3 recluster_egt.py \\
        --egt original.egt \\
        --bpm manifest.bpm \\
        --gtc-dir stage1/gtc/ \\
        --hq-samples high_quality_samples.txt \\
        --output-egt reclustered.egt

Requires: numpy

References:
    - Illumina BeadArrayFiles: https://github.com/Illumina/BeadArrayFiles
    - EGT format: ClusterFile.py in BeadArrayFiles
    - GTC format: GenotypeCalls.py in BeadArrayFiles
"""

import argparse
import os
import struct
import sys
from glob import glob

import numpy as np


# =============================================================================
# Binary I/O helpers (matching Illumina BeadArrayFiles conventions)
# =============================================================================

def read_byte(f):
    return struct.unpack("<B", f.read(1))[0]

def read_int(f):
    return struct.unpack("<i", f.read(4))[0]

def read_float(f):
    return np.frombuffer(f.read(4), dtype=np.float32)[0]

def read_string(f):
    total_length = 0
    partial_length = read_byte(f)
    num_bytes = 0
    while partial_length & 0x80 > 0:
        total_length += (partial_length & 0x7F) << (7 * num_bytes)
        partial_length = struct.unpack("<B", f.read(1))[0]
        num_bytes += 1
    total_length += partial_length << (7 * num_bytes)
    result = f.read(total_length)
    return result.decode("utf-8")

def write_byte(f, val):
    f.write(struct.pack("<B", val))

def write_int(f, val):
    f.write(struct.pack("<i", val))

def write_float(f, val):
    f.write(struct.pack("<f", float(val)))

def write_string(f, s):
    """Write a .NET-style length-prefixed string."""
    encoded = s.encode("utf-8")
    length = len(encoded)
    # Write length as variable-length integer
    while length >= 0x80:
        write_byte(f, (length & 0x7F) | 0x80)
        length >>= 7
    write_byte(f, length)
    f.write(encoded)


# =============================================================================
# EGT (Cluster File) Reader / Writer
# =============================================================================

class ClusterStats:
    """Statistics for a single genotype cluster (AA, AB, or BB)."""
    def __init__(self, theta_mean=0.0, theta_dev=0.0, r_mean=0.0, r_dev=0.0, N=0):
        self.theta_mean = theta_mean
        self.theta_dev = theta_dev
        self.r_mean = r_mean
        self.r_dev = r_dev
        self.N = N


class ClusterRecord:
    """Cluster data for a single locus/probe."""
    def __init__(self):
        self.aa = ClusterStats()
        self.ab = ClusterStats()
        self.bb = ClusterStats()
        self.intensity_threshold = 0.0
        self.cluster_separation = 0.0
        self.total_score = 0.0
        self.original_score = 0.0
        self.edited = False
        self.address = 0
        self.locus_name = ""
        self.genotype_string = ""
        # Extra floats in the data block
        self.extra_floats = []


class EGTFile:
    """Read and write Illumina EGT cluster files."""

    def __init__(self):
        self.gencall_version = ""
        self.cluster_version = ""
        self.call_version = ""
        self.normalization_version = ""
        self.date_created = ""
        self.manifest_name = ""
        self.data_block_version = 0
        self.opa = ""
        self.records = []  # list of ClusterRecord

    @staticmethod
    def read(filepath):
        egt = EGTFile()
        with open(filepath, "rb") as f:
            version = read_int(f)
            if version != 3:
                raise ValueError(f"EGT version {version} not supported (expected 3)")

            egt.gencall_version = read_string(f)
            egt.cluster_version = read_string(f)
            egt.call_version = read_string(f)
            egt.normalization_version = read_string(f)
            egt.date_created = read_string(f)

            is_wgt = read_byte(f)
            if is_wgt != 1:
                raise ValueError("Only WGT cluster files supported")

            egt.manifest_name = read_string(f)
            egt.data_block_version = read_int(f)
            egt.opa = read_string(f)
            num_records = read_int(f)

            # Read cluster records
            records = []
            for _ in range(num_records):
                rec = ClusterRecord()
                # N counts
                aa_n, ab_n, bb_n = read_int(f), read_int(f), read_int(f)
                # R deviations
                aa_r_dev, ab_r_dev, bb_r_dev = read_float(f), read_float(f), read_float(f)
                # R means
                aa_r_mean, ab_r_mean, bb_r_mean = read_float(f), read_float(f), read_float(f)
                # Theta deviations
                aa_t_dev, ab_t_dev, bb_t_dev = read_float(f), read_float(f), read_float(f)
                # Theta means
                aa_t_mean, ab_t_mean, bb_t_mean = read_float(f), read_float(f), read_float(f)

                rec.aa = ClusterStats(aa_t_mean, aa_t_dev, aa_r_mean, aa_r_dev, aa_n)
                rec.ab = ClusterStats(ab_t_mean, ab_t_dev, ab_r_mean, ab_r_dev, ab_n)
                rec.bb = ClusterStats(bb_t_mean, bb_t_dev, bb_r_mean, bb_r_dev, bb_n)

                # Read extra fields based on data block version
                extra = []
                if egt.data_block_version == 9:
                    rec.intensity_threshold = read_float(f)
                    for _ in range(14):
                        extra.append(read_float(f))
                elif egt.data_block_version >= 6:
                    rec.intensity_threshold = 0.0
                    for _ in range(15):
                        extra.append(read_float(f))
                elif egt.data_block_version < 5:
                    rec.intensity_threshold = 0.0
                    extra.append(read_float(f))
                    extra.append(read_float(f))
                else:
                    raise ValueError(f"Unsupported data block version {egt.data_block_version}")

                rec.extra_floats = extra
                records.append(rec)

            # Read cluster scores
            for rec in records:
                rec.cluster_separation = read_float(f)
                rec.total_score = read_float(f)
                rec.original_score = read_float(f)
                rec.edited = read_byte(f) != 0

            # Read genotype strings
            for rec in records:
                rec.genotype_string = read_string(f)

            # Read locus names
            for rec in records:
                rec.locus_name = read_string(f)

            # Read addresses
            for rec in records:
                rec.address = read_int(f)

            # Read cluster counts (should match N values)
            for rec in records:
                counts = [read_int(f), read_int(f), read_int(f)]
                assert rec.aa.N == counts[0], f"AA count mismatch for {rec.locus_name}"
                assert rec.ab.N == counts[1], f"AB count mismatch for {rec.locus_name}"
                assert rec.bb.N == counts[2], f"BB count mismatch for {rec.locus_name}"

            egt.records = records
        return egt

    def write(self, filepath):
        with open(filepath, "wb") as f:
            write_int(f, 3)  # version
            write_string(f, self.gencall_version)
            write_string(f, self.cluster_version)
            write_string(f, self.call_version)
            write_string(f, self.normalization_version)
            write_string(f, self.date_created)
            write_byte(f, 1)  # is_wgt
            write_string(f, self.manifest_name)
            write_int(f, self.data_block_version)
            write_string(f, self.opa)
            write_int(f, len(self.records))

            # Write cluster records
            for rec in self.records:
                # N counts
                write_int(f, rec.aa.N)
                write_int(f, rec.ab.N)
                write_int(f, rec.bb.N)
                # R deviations
                write_float(f, rec.aa.r_dev)
                write_float(f, rec.ab.r_dev)
                write_float(f, rec.bb.r_dev)
                # R means
                write_float(f, rec.aa.r_mean)
                write_float(f, rec.ab.r_mean)
                write_float(f, rec.bb.r_mean)
                # Theta deviations
                write_float(f, rec.aa.theta_dev)
                write_float(f, rec.ab.theta_dev)
                write_float(f, rec.bb.theta_dev)
                # Theta means
                write_float(f, rec.aa.theta_mean)
                write_float(f, rec.ab.theta_mean)
                write_float(f, rec.bb.theta_mean)

                # Extra fields
                if self.data_block_version == 9:
                    write_float(f, rec.intensity_threshold)
                    for val in rec.extra_floats[:14]:
                        write_float(f, val)
                    # Pad if needed
                    for _ in range(14 - len(rec.extra_floats)):
                        write_float(f, 0.0)
                elif self.data_block_version >= 6:
                    for val in rec.extra_floats[:15]:
                        write_float(f, val)
                    for _ in range(15 - len(rec.extra_floats)):
                        write_float(f, 0.0)
                elif self.data_block_version < 5:
                    for val in rec.extra_floats[:2]:
                        write_float(f, val)
                    for _ in range(2 - len(rec.extra_floats)):
                        write_float(f, 0.0)

            # Write cluster scores
            for rec in self.records:
                write_float(f, rec.cluster_separation)
                write_float(f, rec.total_score)
                write_float(f, rec.original_score)
                write_byte(f, 1 if rec.edited else 0)

            # Write genotype strings
            for rec in self.records:
                write_string(f, rec.genotype_string)

            # Write locus names
            for rec in self.records:
                write_string(f, rec.locus_name)

            # Write addresses
            for rec in self.records:
                write_int(f, rec.address)

            # Write cluster counts
            for rec in self.records:
                write_int(f, rec.aa.N)
                write_int(f, rec.ab.N)
                write_int(f, rec.bb.N)


# =============================================================================
# GTC Reader (minimal, for extracting intensities and genotypes)
# =============================================================================

class GTCFile:
    """Minimal GTC reader for extracting genotypes and normalized intensities."""

    # Table of contents entry IDs
    _ID_NUM_SNPS = 1
    _ID_SAMPLE_NAME = 10
    _ID_NORMALIZATION_TRANSFORMS = 400
    _ID_RAW_X = 1000
    _ID_RAW_Y = 1001
    _ID_GENOTYPES = 1002

    def __init__(self, filepath):
        self.filepath = filepath
        self.toc = {}
        self.version = 0
        with open(filepath, "rb") as f:
            identifier = f.read(3).decode("utf-8")
            if identifier != "gtc":
                raise ValueError(f"Not a GTC file: {filepath}")
            self.version = read_byte(f)
            num_toc = read_int(f)
            for _ in range(num_toc):
                entry_id, offset = struct.unpack("<hI", f.read(6))
                self.toc[entry_id] = offset

    @property
    def num_snps(self):
        with open(self.filepath, "rb") as f:
            f.seek(self.toc[self._ID_NUM_SNPS])
            return read_int(f)

    @property
    def sample_name(self):
        with open(self.filepath, "rb") as f:
            f.seek(self.toc[self._ID_SAMPLE_NAME])
            return read_string(f)

    def get_genotypes(self):
        """Return array of genotype codes (0=NC, 1=AA, 2=AB, 3=BB)."""
        with open(self.filepath, "rb") as f:
            f.seek(self.toc[self._ID_GENOTYPES])
            n = read_int(f)
            return np.frombuffer(f.read(n), dtype=np.uint8)

    def get_raw_x(self):
        """Return array of raw X intensities."""
        with open(self.filepath, "rb") as f:
            f.seek(self.toc[self._ID_RAW_X])
            n = read_int(f)
            return np.frombuffer(f.read(n * 2), dtype=np.uint16)

    def get_raw_y(self):
        """Return array of raw Y intensities."""
        with open(self.filepath, "rb") as f:
            f.seek(self.toc[self._ID_RAW_Y])
            n = read_int(f)
            return np.frombuffer(f.read(n * 2), dtype=np.uint16)

    def get_normalization_transforms(self):
        """Return list of normalization transforms."""
        with open(self.filepath, "rb") as f:
            f.seek(self.toc[self._ID_NORMALIZATION_TRANSFORMS])
            n = read_int(f)
            transforms = []
            for _ in range(n):
                # Each transform: version(int), offset_x(f), offset_y(f),
                # scale_x(f), scale_y(f), shear(f), theta(f)
                _ver = read_int(f)
                offset_x = read_float(f)
                offset_y = read_float(f)
                scale_x = read_float(f)
                scale_y = read_float(f)
                shear = read_float(f)
                theta = read_float(f)
                # Reserved bytes (depends on version)
                transforms.append({
                    "offset_x": offset_x, "offset_y": offset_y,
                    "scale_x": scale_x, "scale_y": scale_y,
                    "shear": shear, "theta": theta,
                })
            return transforms


# =============================================================================
# BPM Reader (minimal, for normalization lookup table)
# =============================================================================

def _read_bpm_string(f):
    """Read a .NET-style length-prefixed string from BPM, tolerating binary data.

    Some BPM fields (e.g. SourceSeq, TopGenomicSeq) may contain raw bytes
    that aren't valid UTF-8. Fall back to latin-1 for those.
    """
    total_length = 0
    partial_length = read_byte(f)
    num_bytes = 0
    while partial_length & 0x80 > 0:
        total_length += (partial_length & 0x7F) << (7 * num_bytes)
        partial_length = struct.unpack("<B", f.read(1))[0]
        num_bytes += 1
    total_length += partial_length << (7 * num_bytes)
    raw = f.read(total_length)
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


class BPMFile:
    """Minimal BPM reader for normalization ID lookups.

    Supports BPM versions 1–4+ with per-locus entry versions 4–8.
    Based on the Illumina BeadArrayFiles reference implementation
    (BeadPoolManifest.py / LocusEntry.py).
    """

    def __init__(self, filepath):
        self.filepath = filepath
        self.names = []
        self.norm_ids = []
        self.num_loci = 0
        self._read()

    def _read(self):
        with open(self.filepath, "rb") as f:
            magic = f.read(3).decode("utf-8")
            if magic != "BPM":
                raise ValueError(f"Not a BPM file: {self.filepath}")
            _version = read_byte(f)
            _version_int = read_int(f)
            manifest_name = _read_bpm_string(f)
            _controls = _read_bpm_string(f)
            self.num_loci = read_int(f)

            print(f"  BPM version: {_version}.{_version_int}, manifest: {manifest_name}",
                  file=sys.stderr)

            # Read locus entries — field layout depends on per-entry version.
            # Reference: Illumina BeadArrayFiles LocusEntry.py
            #
            # Per-entry version field layout:
            #   All versions (≥4):
            #     version(int), ilmn_id(str), name(str),
            #     3 strings (snp, ilmn_strand, customer_strand),
            #     address_a(int), allele_a_probe_seq(str),
            #     address_b(int), allele_b_probe_seq(str),
            #     genome_build(str), chr(str), map_info(str),
            #     ploidy(str), species(str), source(str),
            #     source_version(str), source_strand(str),
            #     source_seq(str), top_genomic_seq(str),
            #     bead_set_id(int)
            #   Version ≥ 7: + exp_clusters(byte), intensity_only(byte), assay_type(int)
            #   Version = 6: + exp_clusters(byte), intensity_only(byte), assay_type(byte)
            #   Version < 6 (4,5): + assay_type(byte)
            #   All versions: + norm_id(byte) [inline per-locus]
            self.norm_ids = []
            self.names = []
            for idx in range(self.num_loci):
                entry_version = read_int(f)

                _ilmn_id = _read_bpm_string(f)
                name = _read_bpm_string(f)
                self.names.append(name)

                # 3 strings: SNP, ILMN strand, Customer strand
                for __ in range(3):
                    _read_bpm_string(f)

                read_int(f)              # address_a
                _read_bpm_string(f)      # allele_a_probe_seq
                read_int(f)              # address_b (0 for Infinium II)
                _read_bpm_string(f)      # allele_b_probe_seq

                _read_bpm_string(f)      # genome_build
                _read_bpm_string(f)      # chr
                _read_bpm_string(f)      # map_info (MapInfo / position)
                _read_bpm_string(f)      # ploidy
                _read_bpm_string(f)      # species
                _read_bpm_string(f)      # source
                _read_bpm_string(f)      # source_version
                _read_bpm_string(f)      # source_strand
                _read_bpm_string(f)      # source_seq  (may be binary/long)
                _read_bpm_string(f)      # top_genomic_seq (may be binary/long)
                read_int(f)              # bead_set_id

                if entry_version >= 7:
                    read_byte(f)         # exp_clusters
                    read_byte(f)         # intensity_only
                    read_int(f)          # assay_type (int in v7+)
                elif entry_version == 6:
                    read_byte(f)         # exp_clusters
                    read_byte(f)         # intensity_only
                    read_byte(f)         # assay_type (byte in v6)
                else:
                    # entry_version 4, 5
                    read_byte(f)         # assay_type

                # Normalization ID is stored inline per locus entry
                norm_id = read_byte(f)
                self.norm_ids.append(norm_id)

                if (idx + 1) % 500000 == 0 or idx == 0:
                    if idx == 0:
                        print(f"    Locus entry version: {entry_version}",
                              file=sys.stderr)
                    print(f"    Read {idx + 1}/{self.num_loci} locus entries...",
                          file=sys.stderr)

            print(f"  Read all {self.num_loci} locus entries", file=sys.stderr)
            print(f"  Read {len(self.norm_ids)} normalization IDs (inline)",
                  file=sys.stderr)


def normalize_intensities(raw_x, raw_y, genotypes, norm_ids, transforms):
    """
    Apply Illumina normalization transforms to raw intensities.

    Returns normalized X, Y arrays and derived theta, R arrays.
    """
    n = len(raw_x)
    norm_x = np.zeros(n, dtype=np.float64)
    norm_y = np.zeros(n, dtype=np.float64)

    for i in range(n):
        nid = norm_ids[i] if i < len(norm_ids) else 0
        if nid >= len(transforms):
            nid = 0
        t = transforms[nid]

        # Apply affine normalization (Illumina's approach)
        x = float(raw_x[i]) - t["offset_x"]
        y = float(raw_y[i]) - t["offset_y"]

        # Apply shear and scale
        temp_x = x + t["shear"] * y
        temp_y = y

        norm_x[i] = temp_x * t["scale_x"]
        norm_y[i] = temp_y * t["scale_y"]

        # Clamp to non-negative
        if norm_x[i] < 0:
            norm_x[i] = 0.0
        if norm_y[i] < 0:
            norm_y[i] = 0.0

    # Compute theta and R (Illumina polar coordinates)
    # theta = (2/pi) * arctan2(Y, X), so theta in [0, 1]
    # R = X + Y (total intensity)
    with np.errstate(divide="ignore", invalid="ignore"):
        theta = (2.0 / np.pi) * np.arctan2(norm_y, norm_x)
        r = norm_x + norm_y

    # Handle NaN/inf
    theta = np.nan_to_num(theta, nan=0.0)
    r = np.nan_to_num(r, nan=0.0)

    return norm_x, norm_y, theta, r


def recluster(egt, gtc_files, norm_ids, min_samples_per_cluster=5):
    """
    Recompute cluster statistics from high-quality GTC samples.

    For each probe, aggregates theta and R values across HQ samples
    grouped by genotype, then computes new cluster means and stdevs.
    Falls back to original clusters if insufficient samples per genotype.

    Args:
        egt: Original EGTFile (used as template and fallback)
        gtc_files: List of GTCFile objects for high-quality samples
        norm_ids: List of normalization IDs per probe (from BPM)
        min_samples_per_cluster: Minimum samples needed to recompute a cluster

    Returns:
        EGTFile: New EGT with updated cluster statistics
    """
    num_probes = len(egt.records)
    num_samples = len(gtc_files)

    print(f"  Reclustering {num_probes} probes using {num_samples} samples...", file=sys.stderr)

    # Collect per-probe data across all HQ samples
    # Arrays: [num_probes, num_samples]
    all_theta = np.zeros((num_probes, num_samples), dtype=np.float64)
    all_r = np.zeros((num_probes, num_samples), dtype=np.float64)
    all_geno = np.zeros((num_probes, num_samples), dtype=np.uint8)

    for sidx, gtc in enumerate(gtc_files):
        if (sidx + 1) % 50 == 0 or sidx == 0:
            print(f"    Reading sample {sidx + 1}/{num_samples}...", file=sys.stderr)

        genotypes = gtc.get_genotypes()
        raw_x = gtc.get_raw_x()
        raw_y = gtc.get_raw_y()
        transforms = gtc.get_normalization_transforms()

        _, _, theta, r = normalize_intensities(raw_x, raw_y, genotypes, norm_ids, transforms)

        n = min(num_probes, len(theta))
        all_theta[:n, sidx] = theta[:n]
        all_r[:n, sidx] = r[:n]
        all_geno[:n, sidx] = genotypes[:n]

    # Create new EGT from template
    new_egt = EGTFile()
    new_egt.gencall_version = egt.gencall_version
    new_egt.cluster_version = "reclustered"
    new_egt.call_version = egt.call_version
    new_egt.normalization_version = egt.normalization_version
    new_egt.date_created = egt.date_created
    new_egt.manifest_name = egt.manifest_name
    new_egt.data_block_version = egt.data_block_version
    new_egt.opa = egt.opa

    updated_count = 0
    fallback_count = 0

    for pidx in range(num_probes):
        old_rec = egt.records[pidx]
        new_rec = ClusterRecord()
        new_rec.locus_name = old_rec.locus_name
        new_rec.address = old_rec.address
        new_rec.genotype_string = old_rec.genotype_string
        new_rec.intensity_threshold = old_rec.intensity_threshold
        new_rec.extra_floats = old_rec.extra_floats[:]

        # Preserve original scores
        new_rec.original_score = old_rec.original_score
        new_rec.edited = True  # Mark as reclustered

        genos = all_geno[pidx, :]
        thetas = all_theta[pidx, :]
        rs = all_r[pidx, :]

        probe_updated = False

        # Compute new stats for each genotype class
        for geno_code, old_stats, label in [
            (1, old_rec.aa, "AA"),
            (2, old_rec.ab, "AB"),
            (3, old_rec.bb, "BB"),
        ]:
            mask = genos == geno_code
            n = int(np.sum(mask))

            if n >= min_samples_per_cluster:
                t_vals = thetas[mask]
                r_vals = rs[mask]
                new_stats = ClusterStats(
                    theta_mean=float(np.mean(t_vals)),
                    theta_dev=float(np.std(t_vals, ddof=0)),
                    r_mean=float(np.mean(r_vals)),
                    r_dev=float(np.std(r_vals, ddof=0)),
                    N=n,
                )
                probe_updated = True
            else:
                # Fall back to original cluster stats
                new_stats = ClusterStats(
                    theta_mean=old_stats.theta_mean,
                    theta_dev=old_stats.theta_dev,
                    r_mean=old_stats.r_mean,
                    r_dev=old_stats.r_dev,
                    N=old_stats.N,
                )

            if label == "AA":
                new_rec.aa = new_stats
            elif label == "AB":
                new_rec.ab = new_stats
            else:
                new_rec.bb = new_stats

        # Recompute cluster separation score
        if probe_updated:
            new_rec.cluster_separation = compute_cluster_separation(
                new_rec.aa, new_rec.ab, new_rec.bb
            )
            new_rec.total_score = old_rec.total_score  # Keep GenTrain score
            updated_count += 1
        else:
            new_rec.cluster_separation = old_rec.cluster_separation
            new_rec.total_score = old_rec.total_score
            fallback_count += 1

        new_egt.records.append(new_rec)

    print(f"  Probes updated: {updated_count}", file=sys.stderr)
    print(f"  Probes using original clusters (insufficient data): {fallback_count}", file=sys.stderr)

    return new_egt


def compute_cluster_separation(aa, ab, bb):
    """
    Compute cluster separation score based on distance between clusters
    in theta space, normalized by their spread.
    """
    separations = []
    for c1, c2 in [(aa, ab), (ab, bb), (aa, bb)]:
        if c1.N > 0 and c2.N > 0:
            # Distance between means / sum of deviations
            theta_dist = abs(c1.theta_mean - c2.theta_mean)
            theta_spread = c1.theta_dev + c2.theta_dev
            if theta_spread > 0:
                separations.append(theta_dist / theta_spread)
    if separations:
        return float(np.min(separations))
    return 0.0


def main():
    parser = argparse.ArgumentParser(
        description="Recompute EGT cluster file from high-quality GTC samples"
    )
    parser.add_argument("--egt", required=True,
                        help="Original EGT cluster file")
    parser.add_argument("--bpm", required=True,
                        help="BPM manifest file (for normalization lookups)")
    parser.add_argument("--gtc-dir", required=True,
                        help="Directory containing GTC files from Stage 1")
    parser.add_argument("--hq-samples", required=True,
                        help="File listing high-quality sample IDs (one per line)")
    parser.add_argument("--output-egt", required=True,
                        help="Output path for reclustered EGT file")
    parser.add_argument("--min-cluster-samples", type=int, default=5,
                        help="Minimum samples per genotype cluster to recompute (default: 5)")
    args = parser.parse_args()

    # Read original EGT
    print(f"Reading original EGT: {args.egt}", file=sys.stderr)
    egt = EGTFile.read(args.egt)
    print(f"  {len(egt.records)} probes, manifest: {egt.manifest_name}", file=sys.stderr)

    # Read BPM for normalization IDs
    print(f"Reading BPM: {args.bpm}", file=sys.stderr)
    bpm = BPMFile(args.bpm)
    print(f"  {bpm.num_loci} loci", file=sys.stderr)

    # Read HQ sample list
    with open(args.hq_samples) as f:
        hq_sample_ids = {line.strip() for line in f if line.strip()}
    print(f"High-quality samples: {len(hq_sample_ids)}", file=sys.stderr)

    # Find matching GTC files
    all_gtcs = glob(os.path.join(args.gtc_dir, "*.gtc"))
    gtc_files = []
    for gtc_path in sorted(all_gtcs):
        # Extract sample identifier from GTC filename
        basename = os.path.basename(gtc_path).replace(".gtc", "")
        # Try matching by basename or by reading sample name from GTC
        if basename in hq_sample_ids:
            gtc_files.append(GTCFile(gtc_path))
        else:
            # Try reading sample name from the GTC file
            try:
                gtc = GTCFile(gtc_path)
                if gtc.sample_name in hq_sample_ids:
                    gtc_files.append(gtc)
            except Exception:
                pass

    # If no matches found by ID, try matching all GTC files (user may have
    # provided a list based on barcode-position identifiers)
    if not gtc_files:
        print("  No exact sample ID matches. Trying barcode matching...", file=sys.stderr)
        for gtc_path in sorted(all_gtcs):
            basename = os.path.basename(gtc_path).replace(".gtc", "")
            # Check if any HQ sample ID is a substring of the GTC filename
            for sid in hq_sample_ids:
                if sid in basename or basename in sid:
                    gtc_files.append(GTCFile(gtc_path))
                    break

    print(f"  Matched {len(gtc_files)} GTC files for reclustering", file=sys.stderr)

    if len(gtc_files) < 10:
        print(f"WARNING: Only {len(gtc_files)} GTC files matched. "
              "Reclustering may not be reliable.", file=sys.stderr)

    if len(gtc_files) == 0:
        print("ERROR: No GTC files matched high-quality sample list. "
              "Cannot recluster.", file=sys.stderr)
        sys.exit(1)

    # Perform reclustering
    print("Reclustering...", file=sys.stderr)
    new_egt = recluster(egt, gtc_files, bpm.norm_ids, args.min_cluster_samples)

    # Write new EGT
    print(f"Writing reclustered EGT: {args.output_egt}", file=sys.stderr)
    new_egt.write(args.output_egt)

    # Verify by re-reading
    verify = EGTFile.read(args.output_egt)
    assert len(verify.records) == len(new_egt.records), "Verification failed: record count mismatch"
    print(f"Verification passed: {len(verify.records)} probes", file=sys.stderr)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
