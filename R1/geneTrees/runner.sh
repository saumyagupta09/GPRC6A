for clade in Afrotheria Artiodactyla Carnivora_Pholidota Chiroptera Eulipotyphla Perissodactyla Primates_Scandantia_Dermoptera Rodentia_Lagomorpha
do
cd geneTrees/"$clade"/

awk 'NR==1{next} 
     /^[ACGTUNRYMKSWHBVDacgtunrymkswhbvd\.\-\?]+$/{ printf "%s", $0; next } 
     { if (NR>2) print ""; printf(">%s\n",$0) } 
     END{ print "" }' codon_alignment.phy > codon_alignment.fasta


# 1. Infer the gene tree with codon-aware model
iqtree2 -s codon_alignment.fasta -st CODON -m MFP -bb 1000 -alrt 1000 -nt AUTO

# 2. Compute site concordance factors
iqtree2 -s codon_alignment.fasta -st CODON -m GY -t codon_alignment.fasta.treefile --scf 1000 -nt 10

cd /home/morpheus/GPRC6A_hyphy
done


python3 compare_trees.py --gene_tree Afrotheria/codon_alignment.fasta.treefile --species_tree Afrotheria/Afrotheria_Xenarthra.nwk --output_prefix Afrotheria/GPRC6A_Afrotheria

python3 compare_trees.py --gene_tree Artiodactyla/codon_alignment.fasta.treefile --species_tree Artiodactyla/Artiodactyla.nwk --output_prefix Artiodactyla/GPRC6A_Artiodactyla
python3 compare_trees.py --gene_tree Chiroptera/codon_alignment.fasta.treefile --species_tree Chiroptera/Chiroptera.nwk --output_prefix Chiroptera/GPRC6A_Chiroptera
python3 compare_trees.py --gene_tree Eulipotyphla/codon_alignment.fasta.treefile --species_tree Eulipotyphla/Eulipotyphla.nwk --output_prefix Eulipotyphla/GPRC6A_Eulipotyphla
python3 compare_trees.py --gene_tree Perissodactyla/codon_alignment.fasta.treefile --species_tree Perissodactyla/Perissodactyla.nwk --output_prefix Perissodactyla/GPRC6A_Perissodactyla
python3 compare_trees.py --gene_tree Primates_Scandantia_Dermoptera/codon_alignment.fasta.treefile --species_tree Primates_Scandantia_Dermoptera/Primates_Scandantia_Dermoptera.nwk --output_prefix Primates_Scandantia_Dermoptera/GPRC6A_Primates_Scandantia_Dermoptera
python3 compare_trees.py --gene_tree Rodentia_Lagomorpha/codon_alignment.fasta.treefile --species_tree Rodentia_Lagomorpha/Rodentia_Lagomorpha.nwk --output_prefix Rodentia_Lagomorpha/GPRC6A_Rodentia_Lagomorpha
python3 compare_trees.py --gene_tree Carnivora_Pholidota/codon_alignment.fasta.treefile --species_tree Carnivora_Pholidota/Carnivora_Pholidota.nwk --output_prefix Carnivora_Pholidota/GPRC6A_Carnivora_Pholidota



python ../compare_trees2.py codon_alignment.fasta.treefile Afrotheria_Xenarthra.nwk GPRC6A_Afrotheria --outgroup Elephas_maximus
python ../compare_trees2.py codon_alignment.fasta.treefile Artiodactyla.nwk GPRC6A_Artiodactyla --outgroup Sus_scrofa
python ../compare_trees2.py codon_alignment.fasta.treefile Carnivora_Pholidota.nwk GPRC6A_Carnivora_Pholidota --outgroup Manis_javanica
python ../compare_trees2.py codon_alignment.fasta.treefile Eulipotyphla.nwk GPRC6A_Eulipotyphla --outgroup Condylura_cristata
python ../compare_trees2.py codon_alignment.fasta.treefile Perissodactyla.nwk GPRC6A_Perissodactyla --outgroup Diceros_bicornis
python ../compare_trees2.py codon_alignment.fasta.treefile Chiroptera.nwk GPRC6A_Chiroptera --outgroup Pteropus_alecto
python ../compare_trees2.py codon_alignment.fasta.treefile Rodentia_Lagomorpha.nwk GPRC6A_Rodentia_Lagomorpha --outgroup Oryctolagus_cuniculus
python ../compare_trees2.py codon_alignment.fasta.treefile Primates_Scandentia_Dermoptera.nwk GPRC6A_Primates_Scandentia_Dermoptera --outgroup Galeopterus_variegatus
