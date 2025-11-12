hyphy absrel --tree Primates_Scandentia_Dermoptera.nwk --alignment codon_alignment.phy --output Primates.out > Primates.out.txt

hyphy absrel --srv Yes --tree Primates_Scandentia_Dermoptera.nwk --alignment codon_alignment.phy --output Primates_srv.out > Primates_srv.out.txt

hyphy absrel --multiple-hits Double --srv Yes --tree Primates_Scandentia_Dermoptera.nwk --alignment codon_alignment.phy --output Primates_srv_multi_Double.out > Primates_srv_multi_Double.out.txt

hyphy absrel --multiple-hits Double+Triple --srv Yes --tree Primates_Scandentia_Dermoptera.nwk --alignment codon_alignment.phy --output Primates_srv_multi_Double_Trip.out > Primates_srv_multi_Double_Trip.out.txt
