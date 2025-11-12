for species in `cat Eulipotyphla.txt`
do
for replicate in 1 2 3 4 5 6 7 8 9 10
do
echo $species $replicate
hyphy relax --alignment codon_alignment.phy --tree "$species"_labelled.nwk --test fg --output "$species"_"$replicate"_codon_srv_out --srv Yes > "$species"_"$replicate"_codon_srv.out.txt
hyphy relax --alignment codon_alignment.phy --tree "$species"_labelled.nwk --test fg --output "$species"_"$replicate"_codon_srv_multDoub_out --srv Yes --multiple-hits Double > "$species"_"$replicate"_codon_srv_multDoub.out.txt
hyphy relax --alignment codon_alignment.phy --tree "$species"_labelled.nwk --test fg --output "$species"_"$replicate"_codon_srv_multDoubTrip_out --srv Yes --multiple-hits Double+Triple > "$species"_"$replicate"_codon_srv_multDoubTrip.out.txt
done
done
