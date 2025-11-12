hyphy absrel --tree Artiodactyla.nwk --alignment codon_alignment.phy --output Artiodactyla.out > Artiodactyla.out.txt

hyphy absrel --srv Yes --tree Artiodactyla.nwk --alignment codon_alignment.phy --output Artiodactyla_srv.out > Artiodactyla_srv.out.txt

hyphy absrel --multiple-hits Double --srv Yes --tree Artiodactyla.nwk --alignment codon_alignment.phy --output Artiodactyla_srv_multi_Double.out > Artiodactyla_srv_multi_Double.out.txt

hyphy absrel --multiple-hits Double+Triple --srv Yes --tree Artiodactyla.nwk --alignment codon_alignment.phy --output Artiodactyla_srv_multi_Double_Trip.out > Artiodactyla_srv_multi_Double_Trip.out.txt
