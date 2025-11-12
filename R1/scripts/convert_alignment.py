#!/usr/bin/env python3
"""
Convert between FASTA and PHYLIP alignment formats automatically.

Usage:
    python convert_alignment.py input_file output_file [--relaxed]

Description:
    Automatically detects the input format:
      - FASTA (.fasta, .fa, .aln, .afa) → PHYLIP (.phy)
      - PHYLIP (.phy, .phylip) → FASTA (.fasta)

Options:
    --relaxed   Use relaxed PHYLIP format (long names, for HyPhy)
"""

import sys
import os
import re

# ------------------------
# FASTA → PHYLIP conversion
# ------------------------

def read_fasta(filename):
    """Read a FASTA file and return {name: sequence}."""
    sequences = {}
    with open(filename) as f:
        name, seq = None, []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.replace(" ", ""))
        if name:
            sequences[name] = "".join(seq)
    return sequences


def write_phylip(sequences, outfile, relaxed=False):
    """Write sequences to PHYLIP format (strict or relaxed)."""
    names = list(sequences.keys())
    seqs = list(sequences.values())
    nseq, length = len(seqs), len(seqs[0])

    # Check equal length
    for name, seq in sequences.items():
        if len(seq) != length:
            raise ValueError(f"Error: sequence '{name}' length {len(seq)} != expected {length}")

    with open(outfile, "w") as out:
        out.write(f"{nseq} {length}\n")
        for name, seq in sequences.items():
            if relaxed:
                out.write(f"{name.ljust(15)} {seq}\n")
            else:
                out.write(f"{name[:10].ljust(10)} {seq}\n")


# ------------------------
# PHYLIP → FASTA conversion
# ------------------------

def read_phylip(filename):
    """Read PHYLIP (strict or relaxed, sequential or interleaved)."""
    with open(filename) as f:
        lines = [line.rstrip() for line in f if line.strip()]

    header = lines[0]
    parts = header.split()
    if len(parts) < 2:
        raise ValueError("Invalid PHYLIP header line")
    nseq, nsites = map(int, parts[:2])
    lines = lines[1:]

    seq_dict, seq_order = {}, []
    name_pattern = re.compile(r"^(\S+)\s+([A-Za-z\-?]+)$")

    # Try sequential mode
    for line in lines:
        match = name_pattern.match(line)
        if match:
            name, seq = match.groups()
            seq_dict[name] = seq
            seq_order.append(name)
        else:
            break

    # If not all sequences found, try interleaved mode
    if len(seq_dict) < nseq:
        seq_dict = {name: "" for name in seq_order}
        block_lines = lines.copy()
        while block_lines:
            for name in seq_order:
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

    for name, seq in seq_dict.items():
        if len(seq) != nsites:
            print(f"Warning: sequence '{name}' length {len(seq)} != expected {nsites}")

    return seq_dict


def write_fasta(seq_dict, outfile):
    """Write sequences in FASTA format."""
    with open(outfile, "w") as out:
        for name, seq in seq_dict.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")


# ------------------------
# Auto-detection and conversion
# ------------------------

def detect_format(infile):
    """Guess file format based on extension or content."""
    ext = os.path.splitext(infile)[1].lower()
    if ext in [".fasta", ".fa", ".aln", ".afa"]:
        return "fasta"
    elif ext in [".phy", ".phylip"]:
        return "phylip"

    # Fallback: inspect file content
    with open(infile) as f:
        first = f.readline()
        if first.startswith(">"):
            return "fasta"
        elif re.match(r"^\s*\d+\s+\d+", first):
            return "phylip"
    raise ValueError("Could not detect input format (expected FASTA or PHYLIP).")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]
    relaxed = "--relaxed" in sys.argv

    fmt = detect_format(infile)

    if fmt == "fasta":
        seqs = read_fasta(infile)
        write_phylip(seqs, outfile, relaxed=relaxed)
        print(f"Converted FASTA → PHYLIP ({'relaxed' if relaxed else 'strict'}) → {outfile}")
    elif fmt == "phylip":
        seqs = read_phylip(infile)
        write_fasta(seqs, outfile)
        print(f"Converted PHYLIP → FASTA → {outfile}")
    else:
        raise ValueError("Unrecognized format.")


if __name__ == "__main__":
    main()
