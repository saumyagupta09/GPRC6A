#!/bin/bash

# Usage: ./split_alignment.sh input.fasta 10000
# input.fasta = your large alignment
# 10000 = number of columns per chunk

INPUT=$1
CHUNK_SIZE=$2
PREFIX="chunk_"

# Check if input exists
if [ ! -f "$INPUT" ]; then
    echo "Input file not found!"
    exit 1
fi

# Python script to split alignment
python3 - <<EOF
from Bio import SeqIO

input_fasta = "$INPUT"
chunk_size = int("$CHUNK_SIZE")
prefix = "$PREFIX"

seqs = list(SeqIO.parse(input_fasta, "fasta"))
num_chunks = (len(seqs[0].seq) + chunk_size - 1) // chunk_size

for i in range(num_chunks):
    start = i * chunk_size
    end = min((i + 1) * chunk_size, len(seqs[0].seq))
    records = []
    for s in seqs:
        chunk_seq = s.seq[start:end]
        records.append(s[:0])
        records[-1].seq = chunk_seq
    out_file = f"{prefix}{i+1}.fa"
    SeqIO.write(records, out_file, "fasta")
    print(f"Wrote {out_file} ({start}-{end} columns)")
EOF

