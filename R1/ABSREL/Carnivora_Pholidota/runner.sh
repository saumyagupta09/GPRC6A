hyphy absrel --tree Carnivora_Pholidota.nwk --alignment codon_alignment.phy --output Carnivora.out > Carnivora.out.txt

hyphy absrel --srv Yes --tree Carnivora_Pholidota.nwk --alignment codon_alignment.phy --output Carnivora_srv.out > Carnivora_srv.out.txt

hyphy absrel --multiple-hits Double --srv Yes --tree Carnivora_Pholidota.nwk --alignment codon_alignment.phy --output Carnivora_srv_multi_Double.out > Carnivora_srv_multi_Double.out.txt

hyphy absrel --multiple-hits Double+Triple --srv Yes --tree Carnivora_Pholidota.nwk --alignment codon_alignment.phy --output Carnivora_srv_multi_Double_Trip.out > Carnivora_srv_multi_Double_Trip.out.txt
