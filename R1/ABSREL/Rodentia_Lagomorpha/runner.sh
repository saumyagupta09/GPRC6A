hyphy absrel --tree Rodentia_Lagomorpha.nwk --alignment codon_alignment.phy --output Rodentia.out > Rodentia.out.txt

hyphy absrel --srv Yes --tree Rodentia_Lagomorpha.nwk --alignment codon_alignment.phy --output Rodentia_srv.out > Rodentia_srv.out.txt

hyphy absrel --multiple-hits Double --srv Yes --tree Rodentia_Lagomorpha.nwk --alignment codon_alignment.phy --output Rodentia_srv_multi_Double.out > Rodentia_srv_multi_Double.out.txt

hyphy absrel --multiple-hits Double+Triple --srv Yes --tree Rodentia_Lagomorpha.nwk --alignment codon_alignment.phy --output Rodentia_srv_multi_Double_Trip.out > Rodentia_srv_multi_Double_Trip.out.txt
