hyphy absrel --tree Eulipotyphla.nwk --alignment codon_alignment.phy --output Eulipotyphla.out > Eulipotyphla.out.txt

hyphy absrel --srv Yes --tree Eulipotyphla.nwk --alignment codon_alignment.phy --output Eulipotyphla_srv.out > Eulipotyphla_srv.out.txt

hyphy absrel --multiple-hits Double --srv Yes --tree Eulipotyphla.nwk --alignment codon_alignment.phy --output Eulipotyphla_srv_multi_Double.out > Eulipotyphla_srv_multi_Double.out.txt

hyphy absrel --multiple-hits Double+Triple --srv Yes --tree Eulipotyphla.nwk --alignment codon_alignment.phy --output Eulipotyphla_srv_multi_Double_Trip.out > Eulipotyphla_srv_multi_Double_Trip.out.txt
