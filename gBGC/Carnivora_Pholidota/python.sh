#!/bin/bash

# Usage: ./split_alignment.sh input.fasta 10000
INPUT=$1
CHUNK_SIZE=$2
PREFIX="chunk_"

python3 - <<EOF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

input_fasta = "$INPUT"
chunk_size = int("$CHUNK_SIZE")
prefix = "$PREFIX"

seqs = list(SeqIO.parse(input_fasta, "fasta"))
seq_length = len(seqs[0].seq)
num_chunks = (seq_length + chunk_size - 1) // chunk_size

for i in range(num_chunks):
    start = i * chunk_size
    end = min((i + 1) * chunk_size, seq_length)
    records = []
    for s in seqs:
        chunk_seq = s.seq[start:end]
        records.append(SeqRecord(Seq(str(chunk_seq)), id=s.id, description=""))
    out_file = f"{prefix}{i+1}.fa"
    SeqIO.write(records, out_file, "fasta")
    print(f"Wrote {out_file} ({start}-{end} columns)")
EOF
