#!/usr/bin/env python3
"""
Diagnostic tool to dump raw BPM header bytes and test different parsing strategies.
Run this against a BPM file to determine the correct header format.

Usage:
    python3 bpm_header_probe.py /path/to/file.bpm
"""
import struct
import sys


def read_byte(f):
    data = f.read(1)
    if len(data) < 1:
        raise EOFError("Unexpected EOF")
    return struct.unpack("<B", data)[0]

def read_int(f):
    data = f.read(4)
    if len(data) < 4:
        raise EOFError("Unexpected EOF")
    return struct.unpack("<i", data)[0]

def read_string(f):
    """Read .NET variable-length prefixed string."""
    total_length = 0
    partial_length = read_byte(f)
    num_bytes = 0
    while partial_length & 0x80 > 0:
        total_length += (partial_length & 0x7F) << (7 * num_bytes)
        partial_length = read_byte(f)
        num_bytes += 1
    total_length += partial_length << (7 * num_bytes)
    raw = f.read(total_length)
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")

def peek_bytes(f, n=32):
    """Read n bytes from current position and seek back."""
    pos = f.tell()
    data = f.read(n)
    f.seek(pos)
    return data

def show_bytes(data, label=""):
    """Display bytes as hex and ASCII."""
    hex_str = " ".join(f"{b:02x}" for b in data)
    ascii_str = "".join(chr(b) if 32 <= b < 127 else "." for b in data)
    print(f"  {label:30s} hex: {hex_str}")
    print(f"  {'':30s} asc: {ascii_str}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 bpm_header_probe.py <file.bpm>")
        sys.exit(1)

    filepath = sys.argv[1]

    with open(filepath, "rb") as f:
        # Show first 200 bytes raw
        print(f"=== BPM Header Probe: {filepath} ===")
        file_size = f.seek(0, 2)
        print(f"File size: {file_size} bytes")
        f.seek(0)

        raw_header = f.read(200)
        print(f"\nFirst 200 bytes:")
        for i in range(0, min(200, len(raw_header)), 16):
            chunk = raw_header[i:i+16]
            hex_str = " ".join(f"{b:02x}" for b in chunk)
            ascii_str = "".join(chr(b) if 32 <= b < 127 else "." for b in chunk)
            print(f"  {i:04d}: {hex_str:48s}  {ascii_str}")

        # Parse header
        print(f"\n=== Header: magic(3) + version(byte) + version_flag(int) ===")
        f.seek(0)
        magic = f.read(3).decode("utf-8")
        print(f"  magic: {magic}")
        ver_byte = read_byte(f)
        print(f"  version (byte): {ver_byte}")
        ver_int = read_int(f)
        print(f"  version_flag (int): {ver_int}")

        name = read_string(f)
        print(f"  manifest_name: '{name}' (len={len(name)})")

        controls = read_string(f)
        print(f"  controls: len={len(controls)}")

        num_loci = read_int(f)
        print(f"  num_loci: {num_loci}")

        pos_after_header = f.tell()
        print(f"  Position after header: {pos_after_header}")

        # Verified BPM format:
        #   index_array: num_loci × int32 (1-based sequential)
        #   name_strings: num_loci × .NET-prefixed string
        #   norm_ids: num_loci × byte
        print(f"\n=== Locus data (index array + names + norm IDs) ===")

        # Skip index array
        index_block_size = num_loci * 4
        # Verify indices are sequential
        first_idx = read_int(f)
        second_idx = read_int(f)
        f.seek(pos_after_header)
        print(f"  First two indices: {first_idx}, {second_idx}")
        if first_idx == 1 and second_idx == 2:
            print(f"  Index array looks sequential (1, 2, ...) — OK")
        else:
            print(f"  WARNING: Index array may not be sequential")

        f.seek(index_block_size, 1)
        print(f"  Skipped {index_block_size} bytes of index array")
        print(f"  Position after indices: {f.tell()}")

        # Read first few names
        show_bytes(peek_bytes(f, 40), "next 40 bytes (names start):")
        names = []
        for i in range(min(5, num_loci)):
            n = read_string(f)
            names.append(n)
        print(f"  First 5 names: {names}")

        # Skip remaining names
        print(f"  Skipping remaining {num_loci - 5} names...")
        for i in range(5, num_loci):
            read_string(f)

        pos_after_names = f.tell()
        print(f"  Position after all names: {pos_after_names}")

        # Read norm IDs
        norm_bytes = f.read(num_loci)
        print(f"  Read {len(norm_bytes)} norm ID bytes")
        if len(norm_bytes) == num_loci:
            print(f"  Norm ID range: {min(norm_bytes)}-{max(norm_bytes)}")

        pos_after_norms = f.tell()
        remaining = file_size - pos_after_norms
        print(f"  Position after norm IDs: {pos_after_norms}")
        print(f"  Remaining bytes: {remaining}")
        print(f"\n=== PARSE COMPLETE ===")


if __name__ == "__main__":
    main()

