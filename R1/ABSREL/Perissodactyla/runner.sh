hyphy absrel --tree Perissodactyla.nwk --alignment codon_alignment.phy --output Perissodactyla.out > Perissodactyla.out.txt

hyphy absrel --srv Yes --tree Perissodactyla.nwk --alignment codon_alignment.phy --output Perissodactyla_srv.out > Perissodactyla_srv.out.txt

hyphy absrel --multiple-hits Double --srv Yes --tree Perissodactyla.nwk --alignment codon_alignment.phy --output Perissodactyla_srv_multi_Double.out > Perissodactyla_srv_multi_Double.out.txt

hyphy absrel --multiple-hits Double+Triple --srv Yes --tree Perissodactyla.nwk --alignment codon_alignment.phy --output Perissodactyla_srv_multi_Double_Trip.out > Perissodactyla_srv_multi_Double_Trip.out.txt
