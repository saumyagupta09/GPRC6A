#!/usr/bin/env python3
"""
Convert a multiple sequence alignment in FASTA format to PHYLIP format
compatible with HyPhy and PAML.

Usage:
    python fasta_to_phylip.py input.fasta output.phy [--relaxed]

Options:
    --relaxed   Allow longer sequence names (for HyPhy).
"""

import sys
import textwrap

def read_fasta(filename):
    """Reads a FASTA file and returns an ordered dict of {name: sequence}."""
    sequences = {}
    with open(filename) as f:
        name = None
        seqs = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(seqs)
                name = line[1:].split()[0]  # take first word after '>'
                seqs = []
            else:
                seqs.append(line.replace(" ", ""))
        if name:
            sequences[name] = "".join(seqs)
    return sequences


def write_phylip(sequences, outfile, relaxed=False):
    """Writes sequences to a PHYLIP file."""
    names = list(sequences.keys())
    seqs = list(sequences.values())

    nseq = len(seqs)
    length = len(seqs[0])

    # Validate alignment
    for name, seq in sequences.items():
        if len(seq) != length:
            raise ValueError(f"Sequence '{name}' has length {len(seq)}, expected {length}")

    with open(outfile, "w") as out:
        out.write(f"{nseq} {length}\n")

        for name, seq in sequences.items():
            if relaxed:
                out.write(f"{name.ljust(15)} {seq}\n")
            else:
                # Truncate to 10 characters for strict PHYLIP (PAML)
                out.write(f"{name[:10].ljust(10)} {seq}\n")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]
    relaxed = "--relaxed" in sys.argv

    seqs = read_fasta(infile)
    write_phylip(seqs, outfile, relaxed=relaxed)
    print(f"Converted {len(seqs)} sequences to PHYLIP format -> {outfile}")
    if relaxed:
        print("Note: Using relaxed name format (HyPhy-compatible).")
    else:
        print("Note: Using strict 10-char name format (PAML-compatible).")


if __name__ == "__main__":
    main()
