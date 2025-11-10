#Marsupials selection test

#Using RELAX from Hyphy

#Rename fasta files
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {print "";}' all_Marsupialia.fa > all_Marsupialia_renamed.fa
while read i
do
i1=$(echo $i | awk '{print$1}')
i2=$(echo $i | awk '{print$2}')
sed "/$i2/ s/.*/>$i1/g" all_Marsupialia_renamed.fa -i
done < name_id.txt

#Filter sequence
awk -v RS=">" '!/\<Phascolarctos_cinereus\>/ {print">"$0}' all_Marsupialia_renamed.fa | grep -v "^>$" | awk NF > tmp
mv tmp all_Marsupialia_renamed.fa

#Gene tree
time mafft --auto --maxiterate 1000 all_Marsupialia_renamed.fa > all_Marsupialia_renamed.aln
time ~/aswin/programs/iqtree-2.4.0-Linux-intel/bin/iqtree2 --ufboot 1000 --mem 80G -T AUTO -ntmax 32 -m MFP -s all_Marsupialia_renamed.aln


#Codon alignment

#Run these sequences through pre-msa.bf in order to correct frame-shift mutations and translate the resulting sequences to proteins.
/home/ceglab27/aswin/programs/hyphy-2.5.73/hyphy ~/aswin/programs/hyphy-analyses/codon-msa/pre-msa.bf --input all_Marsupialia_renamed.fa
#Take the output of step 2 and run in through the general MSA program to generate a protein MSA
muscle -in all_Marsupialia_renamed.fa_protein.fas -out all_Marsupialia_renamed.fa_protein.msa
#Run the protein MSA and the frameshift corrected nucleotide sequences from step 2 through post-msa.bf to obtain a nucleotide msa.
/home/ceglab27/aswin/programs/hyphy-2.5.73/hyphy ~/aswin/programs/hyphy-analyses/codon-msa/post-msa.bf --protein-msa all_Marsupialia_renamed.fa_protein.msa \
  --nucleotide-sequences all_Marsupialia_renamed.fa_nuc.fas --output all_Marsupialia_renamed.fa_nuc.msa
sed '/^>/ s/_1//g' all_Marsupialia_renamed.fa_nuc.msa -i

#Label tree
total_start_time=$SECONDS
for i in $(grep ">" all_Marsupialia_renamed.fa | tr -d ">")
do
echo ">" $i
grep -v $i <(grep ">" all_Marsupialia_renamed.fa | tr -d ">") > "$i"_bg
grep $i <(grep ">" all_Marsupialia_renamed.fa | tr -d ">") > "$i"_fg
/home/ceglab27/aswin/programs/hyphy-2.5.73/hyphy ~/aswin/programs/hyphy-analyses/LabelTrees/label-tree.bf --tree all_Marsupialia_renamed.aln.treefile --label bg --list "$i"_bg --output "$i"_all_Marsupialia_renamed.aln.treefile_tmp &> /dev/null
/home/ceglab27/aswin/programs/hyphy-2.5.73/hyphy ~/aswin/programs/hyphy-analyses/LabelTrees/label-tree.bf --tree "$i"_all_Marsupialia_renamed.aln.treefile_tmp --label fg --list "$i"_fg --output "$i"_as_fg_all_Marsupialia_renamed.aln.treefile &> /dev/null
rm "$i"_fg "$i"_bg "$i"_all_Marsupialia_renamed.aln.treefile_tmp
#Run relax
/home/ceglab27/aswin/programs/hyphy-2.5.73/hyphy relax --alignment all_Marsupialia_renamed.fa_nuc.msa --tree "$i"_as_fg_all_Marsupialia_renamed.aln.treefile --test fg --reference bg --starting-points 5 --output "$i"_as_fg_all_Marsupialia_renamed.aln.relax.json > "$i"_as_fg_all_Marsupialia_renamed.aln.relax.out &
done
wait
total_elapsed_time=$(( SECONDS - total_start_time ))
echo -e "\n Total time taken:" && echo $total_elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e


#Summary of relax
for r in $(find . -maxdepth 1 -name "*relax.out" -type f)
do
ka=`egrep "Relaxation/intensification parameter" $r | awk '{print$NF}' | paste -s -d ","`
k1=$(egrep "Relaxation/intensification parameter" $r | awk '{print$NF}' | tail -2 | head -1)
pa=$(egrep "Likelihood ratio test" $r | awk '{print$NF}' | sed 's/\.$//g' | tr -d "*" | paste -s -d ",")
p1=$(egrep "Likelihood ratio test" $r | awk '{print$NF}' | sed 's/\.$//g' | tr -d "*" | tail -1)
di=$(awk -v a="$k1" -v b="$p1" 'BEGIN{if(b>0.05) print"NS"; else if(a==1) print"N"; else if(a>1) print "I"; else print "R" }')
i=$(cat $r | grep -i "evidence for" | sed 's/ /_/g' | sed 's/_among_.*//g' | awk '{if($0~"No_significant") print"NS"; else if($0~"Evidence_for") print$0}' | sed 's/>Evidence_for_//g' | sed 's/_of_selection//g' | tr -d "*" | paste -s -d ",")
rs=$(echo $r | awk -F "/" '{print$NF}')
echo $rs $ka $k1 $pa $p1 $di $i
unset ka k1 pa p1 di i rs
done | sed '1i Test_file all_k-values final_k-value all_p-values final_p-value Interpretation_based_on_final_k Interpretation_from_test' | column -t > relax_summary

#Since one species had extremely high k-value probably due to 2 rates rerun that species
/home/ceglab27/aswin/programs/hyphy-2.5.73/hyphy relax --alignment all_Marsupialia_renamed.fa_nuc.msa --tree Dromiciops_gliroides_as_fg_all_Marsupialia_renamed.aln.treefile --test fg --reference bg --starting-points 5 --rates 2 --output Dromiciops_gliroides_as_fg_all_Marsupialia_renamed.aln.relax.json > Dromiciops_gliroides_as_fg_all_Marsupialia_renamed.aln.relax.out

#make summary again

cd /media/ceglab27/Saumya/GPRC6A/HYPHY/gene_tree
for r in $(find . -name "*.out" -type f)
do
ka=`egrep "Relaxation/intensification parameter" $r | awk '{print$NF}' | paste -s -d ","`
k1=$(egrep "Relaxation/intensification parameter" $r | awk '{print$NF}' | tail -2 | head -1)
pa=$(egrep "Likelihood ratio test" $r | awk '{print$NF}' | sed 's/\.$//g' | tr -d "*" | paste -s -d ",")
p1=$(egrep "Likelihood ratio test" $r | awk '{print$NF}' | sed 's/\.$//g' | tr -d "*" | tail -1)
di=$(awk -v a="$k1" -v b="$p1" 'BEGIN{if(b>0.05) print"NS"; else if(a==1) print"N"; else if(a>1) print "I"; else print "R" }')
i=$(cat $r | grep -i "evidence for" | sed 's/ /_/g' | sed 's/_among_.*//g' | awk '{if($0~"No_significant") print"NS"; else if($0~"Evidence_for") print$0}' | sed 's/>Evidence_for_//g' | sed 's/_of_selection//g' | tr -d "*" | paste -s -d ",")
fg=$(echo $r | awk -F "/" '{print$NF}' | sed 's/_relax_output.out//g' | sed 's/_as_fg_.*.out//g')
bg=$(echo $r | awk -F "/" '{print$2}' | tr "_" ",")
echo $fg $bg $ka $k1 $pa $p1 $di $i
unset ka k1 pa p1 di i fg bg
done | sed '1i Foreground Background all_k-values final_k-value all_p-values final_p-value Interpretation_based_on_final_k Interpretation_from_test' | column -t > whole_relax_summary

awk '{if($8~"NS") print$1,$2,$4,$6,"Non-significant"; else  print$1,$2,$4,$6,$8}' whole_relax_summary | tr " " "\t" > whole_relax_summary_supplementary_version.txt





