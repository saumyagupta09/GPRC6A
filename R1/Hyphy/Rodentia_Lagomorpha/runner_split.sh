
#for spe in `cat Rodentia_Lagomorpha.txt`
#do
#echo $spe
#python ../split_phylogeny_into_groups_fg.py --tree "$spe"_labelled.nwk --aln codon_alignment.phy --n_groups 5 --outdir splits_"$spe" --id_match_mode relaxed
#done

for species in `cat Rodentia_Lagomorpha.txt|head -n 39|tail -n 4`
do
for replicate in 1 2 3 4 5 6 7 8 9 10
do
echo $species $replicate
hyphy relax --alignment splits_"$species"/group_*.phy --tree splits_"$species"/group_*.tree.nwk --test fg --output splits/"$species"_"$replicate"_codon_srv_out --srv Yes > splits/"$species"_"$replicate"_codon_srv.out.txt
hyphy relax --alignment splits_"$species"/group_*.phy --tree splits_"$species"/group_*.tree.nwk --test fg --output splits/"$species"_"$replicate"_codon_srv_multDoub_out --srv Yes --multiple-hits Double > splits/"$species"_"$replicate"_codon_srv_multDoub.out.txt
hyphy relax --alignment splits_"$species"/group_*.phy --tree splits_"$species"/group_*.tree.nwk --test fg --output splits/"$species"_"$replicate"_codon_srv_multDoubTrip_out --srv Yes --multiple-hits Double+Triple > splits/"$species"_"$replicate"_codon_srv_multDoubTrip.out.txt
done
done

