#!/usr/bin/env python3

import re
import argparse
from collections import namedtuple, defaultdict

SeqRec = namedtuple("SeqRec", ["header", "seq", "species", "acc"])

ORG_RE = re.compile(r"\[organism=([^\]]+)\]", re.IGNORECASE)

def parse_fasta(path):
    """Yield (header, sequence) tuples; header includes '>'."""
    header = None
    chunks = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield (header.rstrip(), "".join(chunks).replace(" ", "").replace("\r", "").replace("\t", ""))
                header = line.strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield (header.rstrip(), "".join(chunks).replace(" ", "").replace("\r", "").replace("\t", ""))

def first_token(header):
    """Return the first whitespace-delimited token of the header (without '>')."""
    h = header[1:].strip() if header.startswith(">") else header.strip()
    return h.split()[0] if h else "seq"

def get_species(header):
    """Extract species from [organism=...] tag; fallback to 'UNKNOWN' if missing."""
    m = ORG_RE.search(header)
    return m.group(1).strip() if m else "UNKNOWN"

def seq_length_for_compare(seq):
    """Length used for comparison (ignore gaps and terminal stop)."""
    return len(seq.replace("-", "").replace("*", ""))

def unknowns_count(seq):
    """Count of unknown/ambiguous residues."""
    return sum(ch in {"X","B","J","Z","*","?","U"} for ch in seq.upper())

def gaps_count(seq):
    return seq.count("-")

def choose_better(a: SeqRec, b: SeqRec):
    """
    Return the better record according to:
      longer length -> fewer unknowns -> fewer gaps -> smaller accession
    """
    la = seq_length_for_compare(a.seq)
    lb = seq_length_for_compare(b.seq)
    if la != lb:
        return a if la > lb else b
    ua = unknowns_count(a.seq)
    ub = unknowns_count(b.seq)
    if ua != ub:
        return a if ua < ub else b
    ga = gaps_count(a.seq)
    gb = gaps_count(b.seq)
    if ga != gb:
        return a if ga < gb else b
    return a if a.acc < b.acc else b

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--infile", required=True, help="Input multi-FASTA")
    ap.add_argument("-o", "--outfile", required=True, help="Output FASTA with one sequence per species (longest)")
    ap.add_argument("--minlen", type=int, default=1, help="Discard sequences shorter than this length (default: 1)")
    ap.add_argument("--report", default=None, help="Optional TSV report: species,kept_acc,kept_len,discarded_count")
    args = ap.parse_args()

    best_by_species = {}
    discarded_counts = defaultdict(int)

    for hdr, seq in parse_fasta(args.infile):
        sp = get_species(hdr)
        acc = first_token(hdr)
        if seq_length_for_compare(seq) < args.minlen:
            discarded_counts[sp] += 1
            continue
        rec = SeqRec(header=hdr, seq=seq, species=sp, acc=acc)
        if sp not in best_by_species:
            best_by_species[sp] = rec
        else:
            best_by_species[sp] = choose_better(best_by_species[sp], rec)
            # Count only the one that loses when species already had a candidate
            # (purely informational; not used for selection)
            if best_by_species[sp] is rec:
                discarded_counts[sp] += 1
            else:
                discarded_counts[sp] += 1

    # Write output FASTA
    with open(args.outfile, "w") as out:
        for sp, rec in sorted(best_by_species.items(), key=lambda x: x[0].lower()):
            out.write(f"{rec.header}\n")
            s = rec.seq
            for i in range(0, len(s), 60):
                out.write(s[i:i+60] + "\n")

    # Optional report
    if args.report:
        with open(args.report, "w") as rep:
            rep.write("species\tkept_acc\tkept_len\tdiscarded_count\n")
            for sp, rec in sorted(best_by_species.items(), key=lambda x: x[0].lower()):
                rep.write(f"{sp}\t{rec.acc}\t{seq_length_for_compare(rec.seq)}\t{discarded_counts.get(sp,0)}\n")

    print(f"Kept {len(best_by_species)} species ? {args.outfile}")
    if args.report:
        print(f"Wrote report ? {args.report}")

if __name__ == "__main__":
    main()
