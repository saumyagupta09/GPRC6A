hyphy absrel --tree Chiroptera.nwk --alignment codon_alignment.phy --output Chiroptera.out > Chiroptera.out.txt

hyphy absrel --srv Yes --tree Chiroptera.nwk --alignment codon_alignment.phy --output Chiroptera_srv.out > Chiroptera_srv.out.txt

hyphy absrel --multiple-hits Double --srv Yes --tree Chiroptera.nwk --alignment codon_alignment.phy --output Chiroptera_srv_multi_Double.out > Chiroptera_srv_multi_Double.out.txt

hyphy absrel --multiple-hits Double+Triple --srv Yes --tree Chiroptera.nwk --alignment codon_alignment.phy --output Chiroptera_srv_multi_Double_Trip.out > Chiroptera_srv_multi_Double_Trip.out.txt
