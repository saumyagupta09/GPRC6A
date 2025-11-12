hyphy absrel --tree Afrotheria_Xenarthra.nwk --alignment codon_alignment.phy --output Afrotheria.out > Afrotheria.out.txt

hyphy absrel --srv Yes --tree Afrotheria_Xenarthra.nwk --alignment codon_alignment.phy --output Afrotheria_srv.out > Afrotheria_srv.out.txt

hyphy absrel --multiple-hits Double --srv Yes --tree Afrotheria_Xenarthra.nwk --alignment codon_alignment.phy --output Afrotheria_srv_multi_Double.out > Afrotheria_srv_multi_Double.out.txt

hyphy absrel --multiple-hits Double+Triple --srv Yes --tree Afrotheria_Xenarthra.nwk --alignment codon_alignment.phy --output Afrotheria_srv_multi_Double_Trip.out > Afrotheria_srv_multi_Double_Trip.out.txt
