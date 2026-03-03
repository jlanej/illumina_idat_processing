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
import time
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
            if egt.data_block_version not in (5, 7, 8, 9):
                raise ValueError(
                    f"Data block version {egt.data_block_version} in cluster "
                    f"file not supported (expected 5, 7, 8, or 9)"
                )
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
                # Matches clusterrecord_read() from gtc2vcf/idat2gtc.c:
                #   version >= 7: intensity_threshold + 14 extra floats
                #   version < 5:  no extra fields
                #   version 5:    no extra fields
                #   version 6:    unsupported in reference
                extra = []
                if egt.data_block_version >= 7:
                    rec.intensity_threshold = read_float(f)
                    for _ in range(14):
                        extra.append(read_float(f))
                else:
                    rec.intensity_threshold = 0.0

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

                # Extra fields (matching corrected reader)
                if self.data_block_version >= 7:
                    write_float(f, rec.intensity_threshold)
                    for val in rec.extra_floats[:14]:
                        write_float(f, val)
                    # Pad if needed
                    for _ in range(14 - len(rec.extra_floats)):
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
        """Return list of normalization transforms.

        Each XForm in the GTC file is 52 bytes:
          version(int32) + offset_x(f) + offset_y(f) + scale_x(f) +
          scale_y(f) + shear(f) + theta(f) + cvx(f) + cvy(f) +
          nn12(f) + rr12(f) + taa(f) + tbb(f)
        We must read all 52 bytes per transform to stay aligned.
        """
        with open(self.filepath, "rb") as f:
            f.seek(self.toc[self._ID_NORMALIZATION_TRANSFORMS])
            n = read_int(f)
            transforms = []
            for _ in range(n):
                _ver = read_int(f)
                offset_x = read_float(f)
                offset_y = read_float(f)
                scale_x = read_float(f)
                scale_y = read_float(f)
                shear = read_float(f)
                theta = read_float(f)
                # Read remaining 6 floats (cvx, cvy, nn12, rr12, taa, tbb)
                # to stay aligned for the next transform
                f.read(6 * 4)  # skip 24 bytes
                transforms.append({
                    "offset_x": offset_x, "offset_y": offset_y,
                    "scale_x": scale_x, "scale_y": scale_y,
                    "shear": shear, "theta": theta,
                })
            return transforms

    def read_all(self):
        """Read genotypes, raw intensities, and transforms in a single file open.

        Returns (genotypes, raw_x, raw_y, transforms) — same as calling
        get_genotypes(), get_raw_x(), get_raw_y(), and
        get_normalization_transforms() individually, but avoids opening
        and closing the file four times per sample.
        """
        with open(self.filepath, "rb") as f:
            # Genotypes
            f.seek(self.toc[self._ID_GENOTYPES])
            n = read_int(f)
            genotypes = np.frombuffer(f.read(n), dtype=np.uint8)

            # Raw X
            f.seek(self.toc[self._ID_RAW_X])
            n = read_int(f)
            raw_x = np.frombuffer(f.read(n * 2), dtype=np.uint16)

            # Raw Y
            f.seek(self.toc[self._ID_RAW_Y])
            n = read_int(f)
            raw_y = np.frombuffer(f.read(n * 2), dtype=np.uint16)

            # Normalization transforms (52 bytes each)
            f.seek(self.toc[self._ID_NORMALIZATION_TRANSFORMS])
            n = read_int(f)
            transforms = []
            for _ in range(n):
                _ver = read_int(f)
                offset_x = read_float(f)
                offset_y = read_float(f)
                scale_x = read_float(f)
                scale_y = read_float(f)
                shear = read_float(f)
                theta = read_float(f)
                # Skip remaining 6 floats (cvx, cvy, nn12, rr12, taa, tbb)
                f.read(6 * 4)
                transforms.append({
                    "offset_x": offset_x, "offset_y": offset_y,
                    "scale_x": scale_x, "scale_y": scale_y,
                    "shear": shear, "theta": theta,
                })

        return genotypes, raw_x, raw_y, transforms


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
    """BPM reader for normalization ID lookups and assay type extraction.

    Reads the Illumina BPM (Bead Pool Manifest) binary format to extract
    locus names, normalization IDs, and assay types needed for intensity
    normalization during reclustering.

    Critical implementation detail (from gtc2vcf reference):
      The raw norm_ids from the BPM are modified by assay_type before use:
        norm_ids[idx] += 100 * assay_type
      This separates Infinium I (A/T), Infinium I (G/C), and Infinium II
      probes into distinct normalization bins. A norm_lookups table then
      maps these modified IDs to compact sequential transform indices
      (0, 1, 2, ...) which index into the normalization transforms stored
      in GTC files.

    BPM binary format (verified against HumanOmni2.5-4v1-Multi_B.bpm):
      Header:
        - magic: 3 bytes ("BPM")
        - version: 1 byte
        - version_flag: 4 bytes (int32)
        - manifest_name: .NET length-prefixed string
        - controls: .NET length-prefixed string (control probe config)
        - num_loci: 4 bytes (int32)
      Locus index array:
        - indices: num_loci × 4 bytes (int32, 1-based sequential)
      Locus name strings:
        - names: num_loci × .NET length-prefixed strings
      Normalization IDs:
        - norm_ids: num_loci bytes (one byte per locus)
      Locus entries:
        - per-locus records with assay_type, addresses, etc.
    """

    def __init__(self, filepath):
        self.filepath = filepath
        self.names = []
        self.norm_ids = []        # raw norm_ids (before assay_type modification)
        self.resolved_norm_ids = []  # compact transform indices (after norm_lookups)
        self.norm_lookups = {}    # modified norm_id -> compact transform index
        self.num_loci = 0
        self._read()

    @staticmethod
    def _build_norm_lookups(modified_norm_ids):
        """Build norm_lookups: map modified norm_ids to compact sequential indices.

        Mirrors bpm_norm_lookups() from gtc2vcf/idat2gtc.c:
          - Collect unique modified norm_ids
          - Sort them
          - Assign sequential indices 0, 1, 2, ...
        """
        unique_ids = sorted(set(modified_norm_ids))
        return {nid: idx for idx, nid in enumerate(unique_ids)}

    def _read(self):
        with open(self.filepath, "rb") as f:
            # --- Header ---
            magic = f.read(3).decode("utf-8")
            if magic != "BPM":
                raise ValueError(f"Not a BPM file: {self.filepath}")
            _version = read_byte(f)
            _version_int = read_int(f)
            if _version_int & 0x1000:
                _version_int ^= 0x1000
            manifest_name = _read_bpm_string(f)
            if _version_int > 1:
                _controls = _read_bpm_string(f)
            self.num_loci = read_int(f)

            print(f"  BPM version: {_version}.{_version_int}, manifest: {manifest_name}",
                  file=sys.stderr)
            print(f"  Num loci: {self.num_loci}", file=sys.stderr)

            if self.num_loci <= 0 or self.num_loci > 10_000_000:
                raise ValueError(
                    f"Invalid num_loci={self.num_loci} — BPM header may be "
                    f"corrupt or in an unrecognized format"
                )

            # --- Locus index array ---
            # Array of num_loci int32 values (1-based sequential indices).
            index_data = struct.unpack(f"<{self.num_loci}i",
                                       f.read(self.num_loci * 4))
            print(f"  Read {self.num_loci} index entries", file=sys.stderr)

            # --- Locus names ---
            # num_loci .NET length-prefixed strings (probe/locus names).
            self.names = []
            try:
                for idx in range(self.num_loci):
                    name = _read_bpm_string(f)
                    self.names.append(name)

                    if (idx + 1) % 500000 == 0 or idx == 0:
                        print(f"    Read {idx + 1}/{self.num_loci} locus names...",
                              file=sys.stderr)

                print(f"  Read all {self.num_loci} locus names", file=sys.stderr)
                if self.names:
                    print(f"    First: {self.names[0]}, Last: {self.names[-1]}",
                          file=sys.stderr)
            except Exception as e:
                print(f"  WARNING: Failed reading locus names at entry "
                      f"{len(self.names)}: {e}", file=sys.stderr)
                # Fill remaining names with placeholders
                while len(self.names) < self.num_loci:
                    self.names.append(f"locus_{len(self.names)}")

            # --- Normalization IDs ---
            # Stored as a contiguous array of num_loci bytes immediately
            # after the locus names block.
            raw_norm_ids = [0] * self.num_loci
            try:
                norm_bytes = f.read(self.num_loci)
                if len(norm_bytes) == self.num_loci:
                    raw_norm_ids = list(norm_bytes)
                    max_nid = max(raw_norm_ids) if raw_norm_ids else 0
                    print(f"  Read {len(raw_norm_ids)} normalization IDs "
                          f"(range: 0-{max_nid})", file=sys.stderr)
                else:
                    print(f"  WARNING: Expected {self.num_loci} norm ID bytes, "
                          f"got {len(norm_bytes)}. Using defaults.",
                          file=sys.stderr)
            except Exception as e:
                print(f"  WARNING: Could not read normalization IDs: {e}. "
                      f"Using defaults.", file=sys.stderr)

            # --- Locus entries ---
            # Parse locus entries to extract assay_type and index per probe.
            # This follows the locusentry_read() function from idat2gtc.c.
            # assay_type: 0 = Infinium II, 1 = Infinium I (A/T), 2 = Infinium I (G/C)
            modified_norm_ids = list(raw_norm_ids)  # will be modified in-place
            assay_type_counts = {0: 0, 1: 0, 2: 0}
            locus_entries_parsed = 0

            try:
                for _entry_idx in range(self.num_loci):
                    # Read locus entry version
                    entry_version = read_int(f)
                    if entry_version < 4 or entry_version == 5 or entry_version > 8:
                        print(f"  WARNING: Unsupported locus entry version "
                              f"{entry_version} at entry {_entry_idx}. "
                              f"Stopping locus entry parsing.", file=sys.stderr)
                        break

                    # ilmn_id
                    _read_bpm_string(f)
                    # name
                    _read_bpm_string(f)
                    # 3 skipped strings
                    _read_bpm_string(f)
                    _read_bpm_string(f)
                    _read_bpm_string(f)
                    # index (1-based)
                    locus_index = read_int(f)
                    idx = locus_index - 1
                    # 1 skipped string
                    _read_bpm_string(f)
                    # ilmn_strand
                    _read_bpm_string(f)
                    # snp
                    _read_bpm_string(f)
                    # chrom
                    _read_bpm_string(f)
                    # ploidy
                    _read_bpm_string(f)
                    # species
                    _read_bpm_string(f)
                    # map_info
                    _read_bpm_string(f)
                    # top_genomic_seq (version 4 only, but present in all versions we support)
                    _read_bpm_string(f)
                    # customer_strand
                    _read_bpm_string(f)
                    # address_a (int32)
                    _address_a = read_int(f)
                    # address_b (int32)
                    _address_b = read_int(f)
                    # allele_a_probe_seq
                    _read_bpm_string(f)
                    # allele_b_probe_seq
                    _read_bpm_string(f)
                    # genome_build
                    _read_bpm_string(f)
                    # source
                    _read_bpm_string(f)
                    # source_version
                    _read_bpm_string(f)
                    # source_strand
                    _read_bpm_string(f)
                    # source_seq
                    _read_bpm_string(f)

                    assay_type = 0
                    if entry_version >= 6:
                        # 1 skipped byte
                        f.read(1)
                        # exp_clusters (1 byte)
                        f.read(1)
                        # intensity_only (1 byte)
                        f.read(1)
                        # assay_type (1 byte)
                        assay_type = read_byte(f)

                        if assay_type > 2:
                            print(f"  WARNING: Invalid assay_type {assay_type} "
                                  f"at entry {_entry_idx}. Defaulting to 0.",
                                  file=sys.stderr)
                            assay_type = 0
                    else:
                        # For older BPM versions, infer assay_type from address_b
                        # address_b == 0 means Infinium II (assay_type=0)
                        # address_b != 0 means Infinium I (assay_type inferred)
                        if _address_b != 0:
                            assay_type = 1  # default to Infinium I (A/T)

                    if entry_version >= 7:
                        # frac_a, frac_c, frac_g, frac_t (4 floats = 16 bytes)
                        f.read(16)
                    if entry_version == 8:
                        # ref_strand string
                        _read_bpm_string(f)

                    # Modify norm_id by assay_type (matching idat2gtc.c)
                    if 0 <= idx < self.num_loci:
                        modified_norm_ids[idx] += 100 * assay_type
                        assay_type_counts[assay_type] = \
                            assay_type_counts.get(assay_type, 0) + 1

                    locus_entries_parsed += 1
                    if (locus_entries_parsed) % 500000 == 0:
                        print(f"    Parsed {locus_entries_parsed}/{self.num_loci} "
                              f"locus entries...", file=sys.stderr)

                print(f"  Parsed {locus_entries_parsed} locus entries", file=sys.stderr)
                print(f"  Assay types: Infinium II={assay_type_counts.get(0, 0)}, "
                      f"Infinium I (A/T)={assay_type_counts.get(1, 0)}, "
                      f"Infinium I (G/C)={assay_type_counts.get(2, 0)}",
                      file=sys.stderr)

            except Exception as e:
                print(f"  WARNING: Failed parsing locus entries at entry "
                      f"{locus_entries_parsed}: {e}", file=sys.stderr)
                print(f"  Will use raw norm_ids without assay_type correction. "
                      f"This may produce incorrect normalization.", file=sys.stderr)
                if locus_entries_parsed == 0:
                    modified_norm_ids = list(raw_norm_ids)

            # Build norm_lookups mapping and resolved_norm_ids
            self.norm_ids = raw_norm_ids
            self.norm_lookups = self._build_norm_lookups(modified_norm_ids)
            self.resolved_norm_ids = [
                self.norm_lookups.get(nid, 0) for nid in modified_norm_ids
            ]

            n_unique = len(self.norm_lookups)
            print(f"  Normalization bins: {n_unique} unique transform indices "
                  f"(from {len(set(raw_norm_ids))} raw + assay_type expansion)",
                  file=sys.stderr)


def normalize_intensities(raw_x, raw_y, genotypes, norm_ids, transforms):
    """
    Apply Illumina normalization transforms to raw intensities.

    Matches the reference implementation raw_x_y2norm_x_y() and
    norm_x_y2ilmn_theta_r() from gtc2vcf/idat2gtc.c:

      1. Subtract offset_x, offset_y
      2. Rotate by theta (using cos/sin)
      3. Apply shear correction
      4. Divide by scale_x, scale_y
      5. Clamp to non-negative

    The norm_ids parameter should be the resolved/compact transform indices
    (from BPMFile.resolved_norm_ids), NOT the raw normalization IDs.

    Returns normalized X, Y arrays and derived theta, R arrays.

    Vectorized: uses fancy indexing for bulk operations.
    """
    n = len(raw_x)

    # Build per-probe transform index array
    nids = np.array(norm_ids[:n], dtype=np.int32) if len(norm_ids) >= n else \
        np.zeros(n, dtype=np.int32)
    num_transforms = len(transforms)
    nids = np.clip(nids, 0, num_transforms - 1)

    # Pre-extract transform parameters into arrays for vectorized lookup
    t_offset_x = np.array([t["offset_x"] for t in transforms], dtype=np.float64)
    t_offset_y = np.array([t["offset_y"] for t in transforms], dtype=np.float64)
    t_shear = np.array([t["shear"] for t in transforms], dtype=np.float64)
    t_scale_x = np.array([t["scale_x"] for t in transforms], dtype=np.float64)
    t_scale_y = np.array([t["scale_y"] for t in transforms], dtype=np.float64)
    t_theta = np.array([t["theta"] for t in transforms], dtype=np.float64)

    # Precompute cos/sin of rotation angle per transform
    t_cos = np.cos(t_theta)
    t_sin = np.sin(t_theta)

    # Step 1: Subtract offsets
    temp_x = raw_x.astype(np.float64) - t_offset_x[nids]
    temp_y = raw_y.astype(np.float64) - t_offset_y[nids]

    # Step 2: Rotation
    cos_vals = t_cos[nids]
    sin_vals = t_sin[nids]
    rot_x = cos_vals * temp_x + sin_vals * temp_y
    rot_y = -sin_vals * temp_x + cos_vals * temp_y

    # Step 3: Shear correction
    sheared_x = rot_x - t_shear[nids] * rot_y

    # Step 4: Scale (divide, not multiply)
    # Protect against division by zero
    sx = t_scale_x[nids]
    sy = t_scale_y[nids]
    sx = np.where(sx == 0, 1.0, sx)
    sy = np.where(sy == 0, 1.0, sy)
    norm_x = sheared_x / sx
    norm_y = rot_y / sy

    # Step 5: Clamp to non-negative
    np.maximum(norm_x, 0.0, out=norm_x)
    np.maximum(norm_y, 0.0, out=norm_y)

    # Compute theta and R (Illumina polar coordinates)
    # Matches norm_x_y2ilmn_theta_r() from idat2gtc.c:
    #   theta = (2/pi) * arctan2(Y, X)
    #   R = |X| + |Y|  (but if both < 0, R = FLT_MIN * FLT_EPSILON)
    # After clamping, X and Y are >= 0, so R = X + Y
    # Special case: if both are 0, theta and R are NaN
    both_zero = (norm_x == 0.0) & (norm_y == 0.0)

    with np.errstate(divide="ignore", invalid="ignore"):
        theta = (2.0 / np.pi) * np.arctan2(norm_y, norm_x)
        r = np.abs(norm_x) + np.abs(norm_y)

    # Set both-zero probes to NaN (matching reference)
    theta[both_zero] = np.nan
    r[both_zero] = np.nan

    # For reclustering, replace NaN with 0 so they don't contribute
    theta = np.nan_to_num(theta, nan=0.0)
    r = np.nan_to_num(r, nan=0.0)

    return norm_x, norm_y, theta, r


def _read_single_gtc(args):
    """Worker function for parallel GTC reading.

    Reads a single GTC file and returns normalized theta, R, and genotypes.
    Used by multiprocessing.Pool.map() to parallelize I/O-bound GTC reads.
    """
    gtc_path, norm_ids, num_probes, sidx = args
    try:
        gtc = GTCFile(gtc_path)
        genotypes, raw_x, raw_y, transforms = gtc.read_all()
        _, _, theta, r = normalize_intensities(raw_x, raw_y, genotypes, norm_ids, transforms)
        n = min(num_probes, len(theta))
        return sidx, theta[:n].astype(np.float32), r[:n].astype(np.float32), genotypes[:n], None
    except Exception as e:
        return sidx, None, None, None, str(e)


def recluster(egt, gtc_files, norm_ids, min_samples_per_cluster=5, n_workers=0):
    """
    Recompute cluster statistics from high-quality GTC samples.

    For each probe, aggregates theta and R values across HQ samples
    grouped by genotype, then computes new cluster means and stdevs.
    Falls back to original clusters if insufficient samples per genotype.

    Args:
        egt: Original EGTFile (used as template and fallback)
        gtc_files: List of GTCFile objects for high-quality samples
        norm_ids: List of resolved/compact normalization transform indices
                  per probe (from BPMFile.resolved_norm_ids). These map
                  directly into the normalization transforms array stored
                  in each GTC file.
        min_samples_per_cluster: Minimum samples needed to recompute a cluster

    Returns:
        EGTFile: New EGT with updated cluster statistics
    """
    num_probes = len(egt.records)
    num_samples = len(gtc_files)

    print(f"  Reclustering {num_probes} probes using {num_samples} samples...", file=sys.stderr)
    print(f"  Allocating intensity arrays ({num_probes} x {num_samples})...", file=sys.stderr)

    # Collect per-probe data across all HQ samples
    # Arrays: [num_probes, num_samples]
    # float32 is sufficient (Illumina intensities are single-precision origin)
    # and halves memory vs float64 — critical for thousands of samples.
    all_theta = np.zeros((num_probes, num_samples), dtype=np.float32)
    all_r = np.zeros((num_probes, num_samples), dtype=np.float32)
    all_geno = np.zeros((num_probes, num_samples), dtype=np.uint8)

    print(f"  Reading intensities from {num_samples} GTC files...", file=sys.stderr)
    read_start = time.time()

    # Determine number of workers for parallel GTC reading
    import multiprocessing
    if n_workers <= 0:
        n_workers = min(multiprocessing.cpu_count(), num_samples, 24)
    n_workers = max(1, n_workers)

    if n_workers > 1 and num_samples >= 4:
        print(f"  Using {n_workers} parallel workers for GTC reading", file=sys.stderr)

        # Prepare arguments for parallel workers
        worker_args = [
            (gtc.filepath, norm_ids, num_probes, sidx)
            for sidx, gtc in enumerate(gtc_files)
        ]

        completed = 0
        with multiprocessing.Pool(n_workers) as pool:
            # Use imap_unordered for best throughput (results arrive as completed)
            for sidx, theta, r, geno, error in pool.imap_unordered(
                _read_single_gtc, worker_args, chunksize=max(1, num_samples // (n_workers * 4))
            ):
                if error:
                    print(f"    WARNING: Failed to read GTC {sidx}: {error}", file=sys.stderr)
                    continue
                n = min(num_probes, len(theta))
                all_theta[:n, sidx] = theta
                all_r[:n, sidx] = r
                all_geno[:n, sidx] = geno
                completed += 1
                if completed % 50 == 0 or completed == num_samples:
                    elapsed = time.time() - read_start
                    avg = elapsed / completed
                    remaining = avg * (num_samples - completed) / 60.0
                    print(f"    {completed}/{num_samples} samples read "
                          f"({elapsed:.1f}s elapsed, ~{remaining:.1f}min remaining)",
                          file=sys.stderr)
    else:
        # Serial fallback for small sample counts or single worker
        for sidx, gtc in enumerate(gtc_files):
            sample_start = time.time()

            genotypes, raw_x, raw_y, transforms = gtc.read_all()

            _, _, theta, r = normalize_intensities(raw_x, raw_y, genotypes, norm_ids, transforms)

            n = min(num_probes, len(theta))
            all_theta[:n, sidx] = theta[:n]
            all_r[:n, sidx] = r[:n]
            all_geno[:n, sidx] = genotypes[:n]

            sample_elapsed = time.time() - sample_start
            total_elapsed = time.time() - read_start
            if (sidx + 1) % 10 == 0 or sidx == 0 or sidx == num_samples - 1:
                avg_per_sample = total_elapsed / (sidx + 1)
                remaining = avg_per_sample * (num_samples - sidx - 1)
                remaining_min = remaining / 60.0
                print(f"    Sample {sidx + 1}/{num_samples} "
                      f"({sample_elapsed:.1f}s this sample, "
                      f"{total_elapsed:.1f}s elapsed, "
                      f"~{remaining_min:.1f}min remaining) "
                      f"[{os.path.basename(gtc.filepath)}]", file=sys.stderr)

    total_read = time.time() - read_start
    print(f"  Finished reading all {num_samples} samples in {total_read:.1f}s "
          f"({total_read / num_samples:.1f}s/sample avg)", file=sys.stderr)

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

    print(f"  Computing new cluster statistics for {num_probes} probes...", file=sys.stderr)
    cluster_start = time.time()

    # --- Vectorized cluster statistics computation ---
    # Precompute per-probe, per-genotype counts, means, and stdevs in bulk
    # instead of looping probe-by-probe with NumPy calls on tiny arrays.

    # Per-genotype masks: (num_probes, num_samples) bool arrays
    geno_counts = {}  # geno_code -> (num_probes,) int array
    geno_theta_mean = {}
    geno_theta_std = {}
    geno_r_mean = {}
    geno_r_std = {}

    for geno_code in (1, 2, 3):
        mask = all_geno == geno_code  # (num_probes, num_samples)
        counts = mask.sum(axis=1)     # (num_probes,)
        geno_counts[geno_code] = counts

        # Compute means and stds using masked operations
        # Use NaN-based approach: set non-matching to NaN, then nanmean/nanstd
        theta_masked = np.where(mask, all_theta, np.nan)
        r_masked = np.where(mask, all_r, np.nan)

        with np.errstate(all="ignore"):
            geno_theta_mean[geno_code] = np.nanmean(theta_masked, axis=1)
            geno_theta_std[geno_code] = np.nanstd(theta_masked, axis=1, ddof=0)
            geno_r_mean[geno_code] = np.nanmean(r_masked, axis=1)
            geno_r_std[geno_code] = np.nanstd(r_masked, axis=1, ddof=0)

        # Replace NaN results (from all-NaN slices) with 0
        geno_theta_mean[geno_code] = np.nan_to_num(geno_theta_mean[geno_code], nan=0.0)
        geno_theta_std[geno_code] = np.nan_to_num(geno_theta_std[geno_code], nan=0.0)
        geno_r_mean[geno_code] = np.nan_to_num(geno_r_mean[geno_code], nan=0.0)
        geno_r_std[geno_code] = np.nan_to_num(geno_r_std[geno_code], nan=0.0)

    vectorized_elapsed = time.time() - cluster_start
    print(f"    Bulk statistics computed in {vectorized_elapsed:.1f}s", file=sys.stderr)

    # Determine which probes have at least one genotype class with enough samples
    probe_has_update = np.zeros(num_probes, dtype=bool)
    for geno_code in (1, 2, 3):
        probe_has_update |= geno_counts[geno_code] >= min_samples_per_cluster

    # --- Build cluster records (still a loop, but now only simple assignments) ---
    record_start = time.time()

    for pidx in range(num_probes):
        if (pidx + 1) % 500000 == 0 or pidx == num_probes - 1:
            rec_elapsed = time.time() - record_start
            if pidx > 0:
                avg_per_probe = rec_elapsed / (pidx + 1)
                remaining = avg_per_probe * (num_probes - pidx - 1)
                print(f"    Building record {pidx + 1}/{num_probes} "
                      f"({rec_elapsed:.1f}s elapsed, "
                      f"~{remaining:.1f}s remaining)",
                      file=sys.stderr)

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

        probe_updated = probe_has_update[pidx]

        # Assign cluster stats for each genotype class
        for geno_code, old_stats, label in [
            (1, old_rec.aa, "aa"),
            (2, old_rec.ab, "ab"),
            (3, old_rec.bb, "bb"),
        ]:
            n = int(geno_counts[geno_code][pidx])
            if n >= min_samples_per_cluster:
                new_stats = ClusterStats(
                    theta_mean=float(geno_theta_mean[geno_code][pidx]),
                    theta_dev=float(geno_theta_std[geno_code][pidx]),
                    r_mean=float(geno_r_mean[geno_code][pidx]),
                    r_dev=float(geno_r_std[geno_code][pidx]),
                    N=n,
                )
            else:
                # Fall back to original cluster stats
                new_stats = ClusterStats(
                    theta_mean=old_stats.theta_mean,
                    theta_dev=old_stats.theta_dev,
                    r_mean=old_stats.r_mean,
                    r_dev=old_stats.r_dev,
                    N=old_stats.N,
                )

            setattr(new_rec, label, new_stats)

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
    cluster_total = time.time() - cluster_start
    print(f"  Cluster computation completed in {cluster_total:.1f}s", file=sys.stderr)

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
    parser.add_argument("--workers", type=int, default=0,
                        help="Number of parallel workers for GTC reading (default: auto)")
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
    new_egt = recluster(egt, gtc_files, bpm.resolved_norm_ids,
                        args.min_cluster_samples, n_workers=args.workers)

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
