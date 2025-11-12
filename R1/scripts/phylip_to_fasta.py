#!/usr/bin/env python3
"""
Convert a PHYLIP format alignment (strict or relaxed) to FASTA format.

Usage:
    python phylip_to_fasta.py input.phy output.fasta

Notes:
- Automatically detects strict vs relaxed PHYLIP format.
- Handles multi-line sequence blocks.
- Works with both interleaved and sequential PHYLIP formats.
"""

import sys
import re


def read_phylip(filename):
    """Reads a PHYLIP file and returns a dict of {name: sequence}."""
    with open(filename) as f:
        lines = [line.strip() for line in f if line.strip()]

    # First line gives counts
    header = lines[0]
    parts = header.split()
    if len(parts) < 2:
        raise ValueError("Invalid PHYLIP header line. Should contain number of sequences and sites.")
    nseq, nsites = map(int, parts[:2])

    lines = lines[1:]

    # Detect interleaved vs sequential
    # If number of sequences <= number of lines before next name repeats, it’s sequential
    seq_dict = {}
    seq_order = []
    relaxed_format = False

    # Guess relaxed format (if name seems longer than 10 chars)
    if len(lines) > 0 and len(lines[0].split()[0]) > 10:
        relaxed_format = True

    # Try to parse as sequential first
    name_pattern = re.compile(r"^(\S+)\s+([A-Za-z\-?]+)$")

    for line in lines:
        match = name_pattern.match(line)
        if match:
            name, seq = match.groups()
            seq_dict[name] = seq
            seq_order.append(name)
        else:
            # continuation of previous sequences (interleaved mode)
            break

    # If not all sequences found initially, must be interleaved
    if len(seq_dict) < nseq:
        # Interleaved mode
        seq_dict = {name: "" for name in seq_order}
        block_lines = lines

        while block_lines:
            for i, name in enumerate(seq_order):
                if not block_lines:
                    break
                line = block_lines.pop(0)
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) == 1:
                    seq = parts[0]
                else:
                    seq = parts[1]
                seq_dict[name] += seq

    # Validate
    for name, seq in seq_dict.items():
        if len(seq) != nsites:
            print(f"Warning: sequence '{name}' length {len(seq)} differs from expected {nsites}")

    return seq_dict


def write_fasta(seq_dict, outfile):
    """Writes sequences in FASTA format."""
    with open(outfile, "w") as out:
        for name, seq in seq_dict.items():
            out.write(f">{name}\n")
            for chunk in [seq[i:i+60] for i in range(0, len(seq), 60)]:
                out.write(chunk + "\n")


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]

    seqs = read_phylip(infile)
    write_fasta(seqs, outfile)

    print(f"Converted {len(seqs)} sequences to FASTA format -> {outfile}")


if __name__ == "__main__":
    main()
