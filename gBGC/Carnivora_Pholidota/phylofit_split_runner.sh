#!/bin/bash

# Usage:
# ./phylofit_split_runner.sh alignment.fasta tree.nwk 500
# alignment.fasta = your large alignment
# tree.nwk = your phylogenetic tree
# 500 = number of columns per chunk

INPUT_ALIGNMENT=$1
TREE_FILE=$2
CHUNK_SIZE=$3
PREFIX="chunk_"

if [ -z "$INPUT_ALIGNMENT" ] || [ -z "$TREE_FILE" ] || [ -z "$CHUNK_SIZE" ]; then
    echo "Usage: $0 alignment.fasta tree.nwk chunk_size"
    exit 1
fi

echo "Splitting alignment into chunks of $CHUNK_SIZE columns..."

# Step 1: Split alignment using Python
python3 - <<EOF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

input_fasta = "$INPUT_ALIGNMENT"
chunk_size = int("$CHUNK_SIZE")
prefix = "$PREFIX"

seqs = list(SeqIO.parse(input_fasta, "fasta"))
seq_length = len(seqs[0].seq)
num_chunks = (seq_length + chunk_size - 1) // chunk_size
print(f"Total length: {seq_length}, number of chunks: {num_chunks}")

for i in range(num_chunks):
    start = i * chunk_size
    end = min((i + 1) * chunk_size, seq_length)
    chunk_records = [SeqRecord(Seq(str(s.seq[start:end])), id=s.id, description="") for s in seqs]
    out_file = f"{prefix}{i+1}.fa"
    SeqIO.write(chunk_records, out_file, "fasta")
    print(f"Wrote {out_file}: columns {start}-{end}")
EOF

# Step 2: Run phyloFit on each chunk
echo "Running phyloFit on each chunk..."
for chunk_file in ${PREFIX}*.fa; do
    base=$(basename $chunk_file .fa)
    out_root="neutral_${base}"
    echo "Running phyloFit on $chunk_file -> $out_root"
    phyloFit \
        --msa $chunk_file \
        --tree $TREE_FILE \
        --subst-mod JC69 \
        --out-root $out_root \
        --msa-format FASTA
done

echo "Done. Use 'neutral_chunk_1.mod' as the neutral model for phastBias."

