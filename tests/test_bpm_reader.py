#!/usr/bin/env python3
"""
Test BPM reader against the actual HumanOmni2.5-4v1-Multi_B.bpm file.

Validates the BPM binary format:
  Header: magic(3) + version(byte) + version_flag(int) + manifest_name(str)
          + controls(str) + num_loci(int)
  Index array: num_loci × int32 (1-based sequential)
  Name strings: num_loci × .NET-prefixed string
  Normalization IDs: num_loci × byte

Usage:
    python3 test_bpm_reader.py /path/to/HumanOmni2.5-4v1-Multi_B.bpm
"""
import os
import struct
import sys


def read_byte(f):
    return struct.unpack("<B", f.read(1))[0]

def read_int(f):
    return struct.unpack("<i", f.read(4))[0]

def read_string(f):
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


def test_bpm(filepath):
    results = []
    passed = 0
    failed = 0

    def check(label, condition, detail=""):
        nonlocal passed, failed
        if condition:
            results.append(f"  PASS: {label}")
            passed += 1
        else:
            results.append(f"  FAIL: {label} — {detail}")
            failed += 1

    results.append(f"Testing BPM: {filepath}")
    results.append(f"File size: {os.path.getsize(filepath)} bytes")
    results.append("")

    with open(filepath, "rb") as f:
        # --- Header ---
        magic = f.read(3).decode("utf-8")
        check("Magic is BPM", magic == "BPM", f"got {magic!r}")

        ver = read_byte(f)
        ver_int = read_int(f)
        check("Version byte is 1", ver == 1, f"got {ver}")
        check("Version flag is 4", ver_int == 4, f"got {ver_int}")

        manifest_name = read_string(f)
        check("Manifest name is correct",
              manifest_name == "HumanOmni2.5-4v1-Multi_B.bpm",
              f"got {manifest_name!r}")

        controls = read_string(f)
        check("Controls string is non-empty", len(controls) > 0, f"len={len(controls)}")

        num_loci = read_int(f)
        check("Num loci is 2450000", num_loci == 2450000, f"got {num_loci}")

        # --- Index array ---
        results.append("")
        results.append("  Checking index array...")
        first_idx = read_int(f)
        second_idx = read_int(f)
        check("First index is 1", first_idx == 1, f"got {first_idx}")
        check("Second index is 2", second_idx == 2, f"got {second_idx}")

        # Seek back and skip the full index array
        f.seek(-8, 1)
        f.seek(num_loci * 4, 1)

        # --- Name strings ---
        results.append("")
        results.append("  Reading locus names...")
        names = []
        for i in range(num_loci):
            names.append(read_string(f))

        check("Read correct number of names", len(names) == num_loci,
              f"got {len(names)}")
        check("First name is non-empty", len(names[0]) > 0 and names[0].isprintable(),
              f"got {names[0]!r}")
        check("Last name is non-empty", len(names[-1]) > 0 and names[-1].isprintable(),
              f"got {names[-1]!r}")
        results.append(f"    First: {names[0]}, Last: {names[-1]}")

        # --- Normalization IDs ---
        results.append("")
        results.append("  Reading normalization IDs...")
        norm_bytes = f.read(num_loci)
        check("Norm ID block is correct size", len(norm_bytes) == num_loci,
              f"got {len(norm_bytes)} bytes")

        if len(norm_bytes) == num_loci:
            max_nid = max(norm_bytes)
            min_nid = min(norm_bytes)
            check("Norm IDs are reasonable (max <= 100)", max_nid <= 100,
                  f"max={max_nid}")
            results.append(f"    Norm ID range: {min_nid}-{max_nid}")

    results.append("")
    results.append(f"Results: {passed} passed, {failed} failed")
    return "\n".join(results), failed == 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 test_bpm_reader.py <file.bpm>")
        sys.exit(1)

    output, success = test_bpm(sys.argv[1])
    print(output)
    sys.exit(0 if success else 1)
