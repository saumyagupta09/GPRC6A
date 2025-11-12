#!/usr/bin/env python3
import argparse
import random
from pathlib import Path

# Accept IUPAC DNA + gap, both cases
IUPAC_DNA_GAP = set("ACGTRYKMSWBDHVN?-acgtrykmswbdhvn")

def read_sequential_phy(filepath):
    """
    Reads a PHYLIP-like sequential file where the first line is 'ntaxa nsites',
    then for each taxon: one line with the name, followed by one or more lines
    of sequence until the next taxon name line.
    Returns: taxa (list of names in order), seqs (dict name->sequence string), wrap_len (int)
    """
    lines = Path(filepath).read_text().splitlines()
    if not lines:
        raise ValueError("Empty file.")

    # Parse header
    parts = lines[0].split()
    if len(parts) < 2:
        raise ValueError("Header line must have at least two integers: 'ntaxa nsites'.")
    try:
        ntaxa_hdr, nsites_hdr = int(parts[0]), int(parts[1])
    except ValueError:
        raise ValueError("Header line does not contain two integers.")

    taxa = []
    seqs = {}
    cur = None
    wrap_len = None

    # Detect the typical line wrap length from the first sequence line we encounter
    def is_seq_line(s: str) -> bool:
        s = s.strip()
        return len(s) > 0 and set(s) <= IUPAC_DNA_GAP

    for idx, raw in enumerate(lines[1:], start=2):
        s = raw.strip()
        if not s:
            continue
        if is_seq_line(s):
            if cur is None:
                raise ValueError(f"Sequence data encountered before a taxon name (line {idx}).")
            seqs[cur] = seqs.get(cur, "") + s
            if wrap_len is None:
                wrap_len = len(s)
        else:
            cur = s
            taxa.append(cur)
            if cur in seqs:
                raise ValueError(f"Duplicate taxon name encountered: {cur}")
            seqs[cur] = ""

    if len(taxa) != ntaxa_hdr:
        # Not fatal, but warn the user
        print(f"Warning: header ntaxa={ntaxa_hdr}, but parsed {len(taxa)} taxa.")

    # Validate lengths
    lens = {len(seqs[t]) for t in taxa}
    if len(lens) != 1:
        raise ValueError("Not all sequences have the same length.")
    nsites = lens.pop()
    if nsites_hdr != nsites:
        print(f"Warning: header nsites={nsites_hdr}, but parsed {nsites} sites.")

    # Fallback wrap length if none detected
    if wrap_len is None or wrap_len <= 0:
        wrap_len = 60

    return taxa, seqs, wrap_len

def write_sequential_phy(taxa, seqs, wrap_len, out_path):
    """
    Writes the alignment in the same sequential style:
    - header line: 'ntaxa nsites'
    - one line per taxon name
    - sequence wrapped at wrap_len
    """
    ntaxa = len(taxa)
    nsites = len(seqs[taxa[0]]) if taxa else 0
    with open(out_path, "w") as fh:
        fh.write(f"{ntaxa:4d}{nsites:8d}\n")
        for name in taxa:
            fh.write(f"{name}\n")
            s = seqs[name]
            for i in range(0, len(s), wrap_len):
                fh.write(s[i:i+wrap_len] + "\n")

def sample_codons_without_replacement(full_len_nt, x_codons, rng):
    """
    Returns sorted nucleotide indices for the chosen codons (keeps original order).
    full_len_nt must be divisible by 3.
    """
    if full_len_nt % 3 != 0:
        raise ValueError("Alignment length is not divisible by 3; cannot sample whole codons.")
    n_codons = full_len_nt // 3
    if x_codons > n_codons:
        raise ValueError(f"Requested {x_codons} codons but only {n_codons} codons available.")
    codon_indices = rng.sample(range(n_codons), x_codons)  # without replacement
    codon_indices.sort()  # preserve original left-to-right order in the replicate
    # Expand to nucleotide column indices (0-based positions)
    nt_cols = []
    for c in codon_indices:
        start = 3 * c
        nt_cols.extend([start, start + 1, start + 2])
    return nt_cols

def build_replicate_alignment(taxa, seqs, nt_cols):
    """
    Extracts the specified nucleotide columns for each taxon to create a new alignment.
    """
    out = {}
    for name in taxa:
        s = seqs[name]
        # Efficient join by list comprehension
        out[name] = "".join(s[i] for i in nt_cols)
    return out

def main():
    ap = argparse.ArgumentParser(
        description="Create bootstrap-like replicates by sampling x codons WITHOUT replacement "
                    "from a PHYLIP-like sequential alignment (name line + wrapped sequence lines)."
    )
    ap.add_argument("input", help="Input alignment file (like your .phy).")
    ap.add_argument("--n", type=int, required=True, help="Number of replicate files to create.")
    ap.add_argument("--x", type=int, required=True, help="Number of codons per replicate (without replacement).")
    ap.add_argument("--out_prefix", required=True, help="Output prefix for replicate files.")
    ap.add_argument("--seed", type=int, default=None, help="Random seed (optional, for reproducibility).")
    args = ap.parse_args()

    taxa, seqs, wrap_len = read_sequential_phy(args.input)
    aln_len = len(seqs[taxa[0]])
    rng = random.Random(args.seed)

    # Zero-padding width for replicate numbering
    pad = len(str(max(1, args.n)))

    for i in range(1, args.n + 1):
        nt_cols = sample_codons_without_replacement(aln_len, args.x, rng)
        rep = build_replicate_alignment(taxa, seqs, nt_cols)
        out_path = f"{args.out_prefix}{i:0{pad}d}.phy"
        write_sequential_phy(taxa, rep, wrap_len, out_path)

if __name__ == "__main__":
    main()
