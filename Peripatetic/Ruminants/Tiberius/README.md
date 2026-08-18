ACC="NC_007307.6"

curl -s \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&strand=1&seq_start=34413890&seq_stop=34414375&rettype=fasta&retmode=text" \
> exon1.fa

curl -s \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&strand=1&seq_start=34414381&seq_stop=34414582&rettype=fasta&retmode=text" \
> exon2.fa

curl -s \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&strand=1&seq_start=34436396&seq_stop=34436496&rettype=fasta&retmode=text" \
> exon3.fa

{
    echo ">g7493.t1 CDS NC_007307.6"
    grep -v '^>' exon1.fa | tr -d '\n'
    grep -v '^>' exon2.fa | tr -d '\n'
    grep -v '^>' exon3.fa | tr -d '\n'
    echo
} > g7493.t1.cds.fa

grep -v '^>' g7493.t1.cds.fa | tr -d '\n' | wc -c

python - <<'PY'
from Bio import SeqIO

r = SeqIO.read("g7493.t1.cds.fa", "fasta")
cds = r.seq.upper()

print("CDS length:", len(cds))
print("Divisible by 3:", len(cds) % 3 == 0)
print("Start codon:", cds[:3])
print("Terminal codon:", cds[-3:])

protein = cds.translate()

print("Translated codons:", len(protein))
print("Internal stops:", [
    i + 1 for i, aa in enumerate(protein[:-1]) if aa == "*"
])

with open("g7493.t1.protein.fa", "w") as f:
    f.write(">g7493.t1 Tiberius predicted protein\n")
    f.write(str(protein).rstrip("*") + "\n")

print("\nCDS:")
print(cds)
print("\nProtein:")
print(str(protein))
PY


CDS length: 789
Divisible by 3: True
Start codon: ATG
Terminal codon: TGA
Translated codons: 263
Internal stops: []

CDS:
ATGTGCTTTGAAAAGGAAATGGAGTACCTTGACTCCTTGGCTATCTTGCTCCTGGCCCTCTCTCTGCTGGGAATTCTGTTTGTTCTGGCCATTGGCATAATATTTACAAGAAACCTGAACACACCCGTTGTGAAATCATCTGGGGAATTGATGGTCCGCTACGTAATCCTCTTCTGTCATTTCCTCAACTTCGCTGGCACAGGCTTTTTCATTAGAGAACCACAAAGCTTCACGTGTAAAACCCGGCAGACACTAATCTTGCATGAGCTTTACTCTCTGCATCTCCTACATTCTGATGAAGTCCCTGAAAATTCTGCTGGCCTTCAGCTGTCCAAGCTGCAGAACTTCCTGAAGTGCTTCTATAAACCCATCCCCATCATTTTCACTTGCACGGGCATCTGGTTGTCGTTTGCACCCTCTAGCTCATCTTTGCAGCCCCTGCTGTGGGACAGAATGTCTCCTTGCCCAGAGTCATTATCTTCGAATGAGGGGTCCATTCTTGCATTTGGCTCCATGCTGGGCTATGCTGCCATCCTGGCCTTCATGTGCTTCATTTGTGCCTTCAAAGGCAGGAAATTCCCTGAGAATTACAATGAGGCCAAATTCATAACATTTGGCATGCTCATTTACTTCATAGCTTGGATCACCTTCATCCCCATCTACCACGTTTGGCAAATATATGCTGGTTGGTGCGCACGTCAGACTCCAGCTGCGCATCCCGGGCTCCCCTGGGAGCGGCGCGGGGCGGGGCGGCAGGGGTCCAGGGAGCTCCGGCCCTGGCTGGGCTGA

Protein:
MCFEKEMEYLDSLAILLLALSLLGILFVLAIGIIFTRNLNTPVVKSSGELMVRYVILFCHFLNFAGTGFFIREPQSFTCKTRQTLILHELYSLHLLHSDEVPENSAGLQLSKLQNFLKCFYKPIPIIFTCTGIWLSFAPSSSSLQPLLWDRMSPCPESLSSNEGSILAFGSMLGYAAILAFMCFICAFKGRKFPENYNEAKFITFGMLIYFIAWITFIPIYHVWQIYAGWCARQTPAAHPGLPWERRGAGRQGSRELRPWLG*
