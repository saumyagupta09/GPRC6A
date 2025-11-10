from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

input_fasta = "Carnivora_Pholidota.aln"
chunk_size = 10000
prefix = "chunk_"

seqs = list(SeqIO.parse(input_fasta, "fasta"))
seq_length = len(seqs[0].seq)
num_chunks = (seq_length + chunk_size - 1) // chunk_size
print(f"Total length: {seq_length}, splitting into {num_chunks} chunks")

for i in range(num_chunks):
    start = i * chunk_size
    end = min((i + 1) * chunk_size, seq_length)
    chunk_records = [SeqRecord(Seq(str(s.seq[start:end])), id=s.id, description="") for s in seqs]
    out_file = f"{prefix}{i+1}.fa"
    SeqIO.write(chunk_records, out_file, "fasta")
    print(f"Wrote {out_file}: columns {start}-{end}")

