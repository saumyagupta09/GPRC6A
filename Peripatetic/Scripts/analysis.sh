###############################################################################################################################################################################################################################################################################################################
#GPRC6A

###############################################################################################################################################################################################################################################################################################################


###############################################################################################################################################################################################################################################################################################################
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Genome BLAST

#Make BLAST database
cd /media/aswin/Saumya/GPRC6A/Comment/gblast
cp /media/aswin/Saumya/GPRC6A/Comment/fasta/GPRC6A_whole_genic_region.fa .
cp /media/aswin/Saumya/GPRC6A/Comment/fasta/Bos_taurus_ENST00000310357_exon_wise.fa /media/aswin/Saumya/GPRC6A/Comment/fasta/Bos_taurus_ENST00000310357_exon_wise_with_repeat.fa .

makeblastdb -in GPRC6A_whole_genic_region.fa -out GPRC6A_whole_genic_region.fa -dbtype nucl

#Query: GPRC6A exon-4_exon_wise

#Blast in output format 6 : table format
blastn -task blastn -evalue 0.001 -db GPRC6A_whole_genic_region.fa -query Bos_taurus_ENST00000310357_exon_wise.fa -num_threads 4 \
 -outfmt "6 qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand qseq sseq"| sed '1i Query\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\tQuery_sequence\tSubject_sequence\n' > gprc6a.out
#Blast in output format 3 : full report format
blastn -task blastn -evalue 0.001 -db GPRC6A_whole_genic_region.fa -query Bos_taurus_ENST00000310357_exon_wise.fa -num_threads 4 -outfmt 3 -out gprc6a.outfmt3 -line_length 280
grep -v "_Query_covered_per_sub" gprc6a.out | awk '$13>90' | awk '{print$2,$7,$8,$1,1,$18}' | sed 's/plus/+/g' | sed 's/minus/-/g' | sed 's/[ ]\+/\t/g' > gprc6a_unique.bed

#Query: GPRC6A with repeat

#Blast in output format 6 : table format
blastn -task blastn -evalue 0.001 -db GPRC6A_whole_genic_region.fa -query Bos_taurus_ENST00000310357_exon_wise_with_repeat.fa -num_threads 4 \
  -outfmt "6 qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand qseq sseq"| sed '1i Query\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\tQuery_sequence\tSubject_sequence\n' > gprc6a_with_repeat.out
#Blast in output format 3 : full report format
blastn -task blastn -evalue 0.001 -db GPRC6A_whole_genic_region.fa -query Bos_taurus_ENST00000310357_exon_wise_with_repeat.fa -num_threads 4 -outfmt 3 -out gprc6a_with_repeat.outfmt3 -line_length 280
grep -v "_Query_covered_per_sub" gprc6a_with_repeat.out | awk '$13>90' | awk '{print$2,$7,$8,$1,1,$18}' | sed 's/plus/+/g' | sed 's/minus/-/g' | sed 's/[ ]\+/\t/g' > gprc6a_with_repeat_unique.bed

###############################################################################################################################################################################################################################################################################################################
#Identify repeats

#-----------------------------------------------------------------------------------------------------------------------------
#Repeatmasker:

cd /media/aswin/Saumya/GPRC6A/Comment/repeat_masker
grep Exon4 Bos_taurus_repeatmasker.out | awk '{print"ref|NC_037336.1|:33772346-33876376",$6+60552-1,$7+60552,$11"_"$10,1,$9}' | sed 's/C$/-/g' | sed 's/[ ]\+/\t/g' > manuscript_exon_4_repeat.bed

#online repeatmasker
cd /media/aswin/Saumya/GPRC6A/Comment/repeat_masker/repeatmasker_online_rmblast
awk 'NR>3{print$5,$6,$7,$11"_"$10,1,$9}' OFS="\t" RM2_GPRC6A_whole_genic_region.fa_1769438276.out | awk '{if($NF=="C") print$1,$2,$3,$4,$5,"-"; else print$0}' | tr " " "\t" > GPRC6A_whole_genic_region_repeatmasker_rmblast.bed
cd /media/aswin/Saumya/GPRC6A/Comment/repeat_masker/repeatmasker_online_crossmatch
tar -xf RM2_GPRC6A_whole_genic_region.fa_1769439492.tar
awk 'NR>3{print$5,$6,$7,$11"_"$10,1,$9}' OFS="\t" RM2_GPRC6A_whole_genic_region.fa_1769439492.out | awk '{if($NF=="C") print$1,$2,$3,$4,$5,"-"; else print$0}' | tr " " "\t" > GPRC6A_whole_genic_region_repeatmasker_crossmatch.bed

#rerun offline repeatmasker
cd /media/aswin/Saumya/GPRC6A/Comment/repeat_masker/repeatmasker_offline
cp ../../transciption_start_sites/GPRC6A_whole_genic_region.fa .
#try different parameters
time /home/ceglab25/ajs/maker/exe/RepeatMasker/RepeatMasker -e RMBlast -pa 12 -a -s -nolow -gff -u -html -species "Bos taurus" GPRC6A_whole_genic_region.fa -dir rmblast_species_lib_GPRC6A_whole_genic_region -xsmall
time /home/ceglab25/ajs/maker/exe/RepeatMasker/RepeatMasker -e RMBlast -pa 12 -a -s -nolow -gff -u -html -species "Bovidae" GPRC6A_whole_genic_region.fa -dir rmblast_bovidae_lib_GPRC6A_whole_genic_region -xsmall
time /home/ceglab25/ajs/maker/exe/RepeatMasker/RepeatMasker -e nhmmer -pa 12 -a -s -nolow -gff -u -html -species "Bos taurus" GPRC6A_whole_genic_region.fa -dir nhmmer_species_lib_GPRC6A_whole_genic_region -xsmall
time /home/ceglab25/ajs/maker/exe/RepeatMasker/RepeatMasker -e nhmmer -pa 12 -a -s -nolow -gff -u -html -species "mammalia" GPRC6A_whole_genic_region.fa -dir nhmmer_mammalia_lib_GPRC6A_whole_genic_region -xsmall
time /home/ceglab25/ajs/maker/exe/RepeatMasker/RepeatMasker -e nhmmer -pa 12 -a -s -nolow -gff -u -html -species "vertebrates" GPRC6A_whole_genic_region.fa -dir nhmmer_vertebrate_lib_GPRC6A_whole_genic_region -xsmall

cd /media/aswin/Saumya/GPRC6A/Comment/gblast
cp gprc6a_unique.bed gprc6a_with_repeat_unique.bed /media/aswin/Saumya/GPRC6A/Comment/repeat_masker/repeatmasker_offline/

#-----------------------------------------------------------------------------------------------------------------------------
#SINE-base

cd /media/aswin/Saumya/GPRC6A/Comment/SINE_LINE
awk '/bestfit/{flag=1; next} /RepBase/{flag=0} flag' Bos_taurus_ENST00000310357_exon_4_sinebase_sinebank.out | awk '{print"1",$6,$1,"1",$5}' \
| sed 's/-/ /1' | sed 's/\[R\]/-/g' | sed 's/\[F\]/+/g' | sed 's/[ ]\+/\t/g' | awk '{print"ref|NC_037336.1|:33772346-33876376",$2+60553,$3+60553,$4,$5,$6}' OFS="\t" | awk '{if($2>$3) print$1,$3,$2,$4,$5,$6; else print$0}' OFS="\t" > Bos_taurus_ENST00000310357_exon_4_sinebase_sinebank.bed

awk '/bestfit/{flag=1; next} /Vassetzky/{flag=0} flag' Bos_taurus_ENST00000310357_exon_4_sinebase_linebank.out | awk '{print"1",$6,$1,"1",$5}' \
| sed 's/-/ /1' | sed 's/\[R\]/-/g' | sed 's/\[F\]/+/g' | sed 's/[ ]\+/\t/g' | awk '{print"ref|NC_037336.1|:33772346-33876376",$2+60553,$3+60553,$4,$5,$6}' OFS="\t" | awk '{if($2>$3) print$1,$3,$2,$4,$5,$6; else print$0}' OFS="\t" > Bos_taurus_ENST00000310357_exon_4_sinebase_linebank.bed

#-----------------------------------------------------------------------------------------------------------------------------
#Identify tandem repeats

#Using MISA webserver:

cd /media/aswin/Saumya/GPRC6A/Comment/tandem_repeats
cp ../fasta/Bos_taurus_ENST00000310357_exon_wise_with_repeat.fa .

#Run MISA in webserver: (http://misaweb.ipk-gatersleben.de/)
#set thresholds:

#summary
for e in $(ls exon*.gff); do f=$(cat $e | grep -v "#" | awk '/region/{print$5}'); cat $e | grep -v "#" | grep -v "region" | awk -v f="$f" '{print$0,f,$5-$4+1,($4/f)*100,($5/f)*100}' | awk '{$(NF-1)+=0; $NF+=0}1' CONVFMT="%.1f"; unset f; done | column -t > repeat_summary

#View in IGV
cp ../gblast/gprc6a_with_repeat_unique.bed .
while read e
do
ex=$(echo $e | awk '{print$1}')
el=$(awk -v ex="$ex" '$4==ex {print$2-1}' gprc6a_with_repeat_unique.bed)
echo $e | awk -v el="$el" '{print "ref|NC_037336.1|:33772346-33876376",$4+el,$5+el,$3,1,"+"}'
unset ex el
done < repeat_summary | sed 's/[ ]\+/\t/g' > tandem_repeats.bed


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Identify stemloop

cd /media/aswin/Saumya/GPRC6A/Comment/stemloop
cp ../gblast/Bos_taurus_ENST00000310357_exon_wise_with_repeat.fa .
cp /media/aswin/Saumya/GPRC6A/Comment/gblast/gprc6a_with_repeat_unique.bed .

#For exon4
grep exon_4 Bos_taurus_ENST00000310357_exon_wise_with_repeat.fa -A1 > exon_4_Bos_taurus_ENST00000310357_exon_wise_with_repeat.fa
einverted -sequence exon_4_Bos_taurus_ENST00000310357_exon_wise_with_repeat.fa -match 3 -mismatch -4 -gap 6 -threshold 15 -outfile outfile.out -outseq outseq.fa
grep "|" -C1 outfile.out | egrep -v "\||^\-" | awk '{print$1,$NF}' | awk '{print"ref|NC_037336.1|:33772346-33876376",$1+60552,$2+60552,"pair_"int((NR+1)/2),"1"}' | awk '{if($2>$3) print$1,$3,$2,$4,$5,"-"; else print$0,"+"}' | sed 's/[ ]\+/\t/'g > exon_4_stemloop.bed

#For whole gene
einverted -sequence Bos_taurus_ENST00000310357_exon_wise_with_repeat.fa -match 3 -mismatch -4 -gap 6 -threshold 15 -outfile whole_cds_stemloop.out -outseq whole_cds_stemloop.fa
#grep "|" -C1 whole_cds_stemloop.out | egrep -v "\||^\-" | awk '{print$1,$NF}' | awk '{print"ref|NC_037336.1|:33772346-33876376",$1+60552,$2+60552,"pair_"int((NR+1)/2),"1"}' | awk '{if($2>$3) print$1,$3,$2,$4,$5,"-"; else print$0,"+"}' | sed 's/[ ]\+/\t/'g > whole_cds_stemloop.bed

cp ../gblast/gprc6a_with_repeat_unique.bed .

sed 's/exon/>exon/g' whole_cds_stemloop.out | grep -v "|" | sed '/exon/ s/:.*//g' | awk -v RS=">" '{$1=$1}1' | awk '{print$1,$2,$4"\n"$1,$5,$7}' | awk NF > whole_cds_stemloop_formatted.out
while read e
do
ex=$(echo $e | awk '{print$1}')
el=$(awk -v ex="$ex" '$4==ex {print$2-1}' gprc6a_with_repeat_unique.bed)
echo $e | awk -v el="$el" '{print "ref|NC_037336.1|:33772346-33876376",$2+el,$3+el,$1}'
unset ex el
done < whole_cds_stemloop_formatted.out | awk '{print$0"_pair_"int((NR+1)/2),"1"}' | awk '{if($2>$3) print$1,$3,$2,$4,$5,"-"; else print$0,"+"}' | sed 's/[ ]\+/\t/g' > whole_cds_stemloop_formatted.bed

###############################################################################################################################################################################################################################################################################################################
#Identify transcription start sites

bedtools getfasta -fi GPRC6A_whole_genic_region.fa -bed <(echo -e "ref|NC_037336.1|:33772346-33876376\t1\t46426\twhole_upstream\t1\t+") -s -name+ > GPRC6A_whole_upstream_region.fa
grep prediction GPRC6A_whole_upstream_region_promoters.out | sed 's/Marginal prediction/MP/g' | sed 's/Highly likely prediction/HLP/g' | awk '{print"ref|NC_037336.1|:33772346-33876376",$1,$1+1,$3"_"$2,"1","+"}' OFS="\t" > GPRC6A_whole_upstream_region_promoters.bed
awk '$4~"HLP"' GPRC6A_whole_upstream_region_promoters.bed > GPRC6A_whole_upstream_region_promoters_HLP.bed

bedtools getfasta -fi GPRC6A_whole_genic_region.fa -bed <(echo -e "ref|NC_037336.1|:33772346-33876376\t44426\t46426\t2000_upstream\t1\t+") -s -name+ > GPRC6A_2000_upstream_region.fa
grep prediction GPRC6A_2000_upstream_region_promoters.out | sed 's/Marginal prediction/MP/g' | sed 's/Highly likely prediction/HLP/g' | awk '{print"ref|NC_037336.1|:33772346-33876376",44426+$1,44426+$1+1,$3"_"$2,"1","+"}' OFS="\t" > GPRC6A_2000_upstream_region_promoters.bed
awk '$4~"HLP"' GPRC6A_2000_upstream_region_promoters.bed > GPRC6A_2000_upstream_region_promoters_HLP.bed

bedtools getfasta -fi GPRC6A_whole_genic_region.fa -bed <(echo -e "ref|NC_037336.1|:33772346-33876376\t46427\t66333\tcds_region\t1\t+") -s -name+ > GPRC6A_cds_region.fa
grep prediction GPRC6A_cds_region_promoters.out | sed 's/Marginal prediction/MP/g' | sed 's/Highly likely prediction/HLP/g' | awk '{print"ref|NC_037336.1|:33772346-33876376",46427+$1,46427+$1+1,$3"_"$2,"1","+"}' OFS="\t" > GPRC6A_cds_region_promoters.bed
awk '$4~"HLP"' GPRC6A_cds_region_promoters.bed > GPRC6A_cds_region_promoters_HLP.bed

#Intron can have TSS: https://pmc.ncbi.nlm.nih.gov/articles/PMC5435436/?utm_source=chatgpt.com

cd /media/aswin/saumya/GPRC6A/Comment/transciption_start_sites/for_ncbi_track
awk '{print$1,$4,$2,$3,$5,$6}' gblast_unique.bed | sed 's/ref|//g' | tr -d "|" | tr " :-" "\t" | awk '{if($NF=="+") print$1,$2+$5,$2+$6,$4,$7,$8}' OFS="\t" > ncbi_gblast_unique.bed

for bed in gblast_unique.bed GPRC6A_2000_upstream_region_promoters.bed GPRC6A_cds_region_promoters.bed GPRC6A_cds_region_promoters_HLP.bed GPRC6A_whole_upstream_region_promoters.bed GPRC6A_whole_upstream_region_promoters_HLP.bed
do
  awk '{print$1,$4,$2,$3,$5,$6}' $bed | sed 's/ref|//g' | tr -d "|" | tr " :-" "\t" | awk '{if($NF=="+") print$1,$2+$5,$2+$6,$4,$7,$8}' OFS="\t" > ncbi_"$bed"
done

awk '{print$0,33818773-$3}' ncbi_GPRC6A_whole_upstream_region_promoters_HLP.bed

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get TSS identified using CAGE from a paper 

#PAPER: https://www.nature.com/articles/s42003-021-02340-6 (Evolution of tissue and developmental specificity of transcription start sites in Bos taurus indicus)
#Download Supplementary data

#The genome used in that paper
cd /media/aswin/saumya/GPRC6A/Comment/RNA_seq/genomes/GCF_002263795.1
time datasets download genome accession GCF_002263795.1 --include gff3,rna,cds,protein,genome,seq-report

time makeblastdb -in GCF_002263795.1_ARS-UCD1.2_genomic.fna -out GCF_002263795.1_ARS-UCD1.2_genomic.fna -dbtype nucl
sortBed -i GCF_002263795.1_ARS-UCD1.2_genomic.gff | gff2bed > GCF_002263795.1_ARS-UCD1.2_genomic.gff.bed

mkdir Bos_taurus_ENST00000310357_exon_wise
cp ../breed_Holstein/Bos_taurus_ENST00000310357_exon_wise/Bos_taurus_ENST00000310357_exon_wise.fa Bos_taurus_ENST00000310357_exon_wise/
time /media/aswin/programs/gblast_short ../GCF_002263795.1_ARS-UCD1.2_genomic.fna Bos_taurus_ENST00000310357_exon_wise.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=no -iflank ../GCF_002263795.1_ARS-UCD1.2_genomic.gff
awk '{print$2,$7,$8,$1,"1",$18}' blast_output/test.out | sed 's/plus/+/g' | sed 's/minus/-/g' | sed '1,2d' | awk '!a[$4]++' | sed 's/[ \t]\+/\t/g' > all_gblast_uniq_hits.bed

#Convert chromosome IDs in TSS file & convert it to bed
for tss in 42003_2021_2340_MOESM10_ESM.txt  42003_2021_2340_MOESM14_ESM.txt  42003_2021_2340_MOESM2_ESM.txt  42003_2021_2340_MOESM3_ESM.txt  42003_2021_2340_MOESM4_ESM.txt  42003_2021_2340_MOESM9_ESM.txt
do
n=$(head -1 $tss | awk '{print$1}')
echo $tss "-" $n
grep "^chr9:" $tss | awk '{$(NF-1)+=0}1' CONVFMT="%.2f" | awk '{$NF+=0}1' CONVFMT="%.2f" | tr ":" "\t" | awk '{print"NC_037336.1",$2,$4"/"$5,".",$3}' \
 | awk '{sub(/-/," ",$2)}1' | awk '{if($6=="" && $NF=="+") print$1,$2,$2+1,$3,$4,$5; else if($6=="" && $NF=="-") print$1,$2-1,$2,$3,$4,$5; else print$0}' | sed 's/[ \t]\+/\t/g' > "CAGE_TSS_"$n".bed"
done

cd /media/aswin/saumya/GPRC6A/Comment/RNA_seq/genomes/GCF_002263795.1
grep "chr9:" 42003_2021_2340_MOESM2_ESM.txt | awk '{$(NF-1)+=0}1' CONVFMT="%.2f" | awk '{$NF+=0}1' CONVFMT="%.2f" | tr ":" "\t" | awk '{print"NC_037336.1",$2,$4"/"$5,".",$3}' | awk '{sub(/-/," ",$2)}1' \
	| awk '{if($6=="" && $NF=="+") print$1,$2,$2+1,$3,$4,$5; else if($6=="" && $NF=="-") print$1,$2-1,$2,$3,$4,$5; else print$0}' OFS="\t" > Bos_taurus_liver_CAGE_TSS.bed


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Download all bed files of CAGE TSS from ftp

#Paper: https://doi.org/10.1093/g3journal/jkad108 (Improving the annotation of the cattle genome by annotating transcription start sites in a diverse set of tissues and populations using Cap Analysis Gene Expression sequencing )

#Parallel Downloads via aria2
cd /media/aswin/saumya/GPRC6A/Comment/transciption_start_sites/CAGE_TSS/tissues_TSS
URL="https://api.faang.org/files/trackhubs/BOVREG_CAGE_EUROFAANG/ARS-UCD1.2/tissues_TSS/"
wget -qO- "$URL" | grep -oP 'href="\K[^"]+' | grep -v '^/' | sed "s|^|$URL|" | aria2c -i - -j 10 -x 5

cd /media/aswin/saumya/GPRC6A/Comment/transciption_start_sites/CAGE_TSS/tissues_TSS_Enhancers
URL="https://api.faang.org/files/trackhubs/BOVREG_CAGE_EUROFAANG/ARS-UCD1.2/tissues_TSS-Enhancers/"
time wget -qO- "$URL" | grep -oP 'href="\K[^"]+' | grep -v '^/' | sed "s|^|$URL|" | aria2c -i - -j 10 -x 5

#Convert bigbed to bed
cd /media/aswin/saumya/GPRC6A/Comment/transciption_start_sites/CAGE_TSS/tissues_TSS
time for b in $(ls TSS_*.Bigbed)
do
d=$(echo "$b" | sed 's/\.Bigbed/\.bed/g')
echo $d
bigBedToBed $b $d
unset d
done 

#Convert bigbed to bed
cd /media/aswin/saumya/GPRC6A/Comment/transciption_start_sites/CAGE_TSS/tissues_TSS_Enhancers
time for b in $(ls BC_*.Bigbed)
do
d=$(echo "$b" | sed 's/\.Bigbed/\.bed/g')
echo $d
bigBedToBed $b $d
unset d
done 

#Download whole supplementary data from paper: https://api.faang.org/files/trackhubs/BOVREG_CAGE_EUROFAANG/
cd /media/aswin/saumya/GPRC6A/Comment/transciption_start_sites/CAGE_TSS/Supplementary
cd /media/aswin/saumya/GPRC6A/Comment/transciption_start_sites/CAGE_TSS/Supplementary/21769649/Supplementary_file2_CoExp_links_Stretches_superenhencers
awk '{print$1,$2,$3,$6"_"$4}' OFS="\t" Enhancer_stretches_10Kbp_min3_nonpervasive.tsv > Enhancer_stretches_10Kbp_min3_nonpervasive.bed
cd /media/aswin/saumya/GPRC6A/Comment/transciption_start_sites/CAGE_TSS/Supplementary/21769649/Supplementary_file3_Tissue_level_GFF3_All-in-One
gzip -d *

#genome used in paper: ARS-UCD1.2 (GCF_002263795.1 / GCA_002263795.2 , breed : Hereford ; Apr 11, 2018)
#But now this genome is labelled: Warning: contaminated & See latest version: GCF_002263795.3

#link: https://sites.ualberta.ca/~stothard/1000_bull_genomes/
cd /media/aswin/saumya/GPRC6A/Comment/transciption_start_sites/CAGE_TSS/genome
wget https://sites.ualberta.ca/~stothard/1000_bull_genomes/ARS-UCD1.2_Btau5.0.1Y.fa.gz
wget https://sites.ualberta.ca/~stothard/1000_bull_genomes/ARS-UCD1.2_Btau5.0.1Y.ann
gzip -d ARS-UCD1.2_Btau5.0.1Y.fa.gz
samtools faidx ARS-UCD1.2_Btau5.0.1Y.fa

#Run genome blast & find GPRC6A gene
mkdir Bos_taurus_ENST00000310357_exon_wise
cp /media/aswin/saumya/GPRC6A/Comment/RNA_seq/Bos_taurus_ENST00000310357_exon_wise.fa Bos_taurus_ENST00000310357_exon_wise/
time /media/aswin/programs/gblast_short ../ARS-UCD1.2_Btau5.0.1Y.fa Bos_taurus_ENST00000310357_exon_wise.fa -evalue=0.05 -word_size=11 -fix_query -tblastx=yes -iflank
awk '{print$2,$7,$8,$1,"1",$18}' blast_output/test.out | sed 's/plus/+/g' | sed 's/minus/-/g' | sed '1,2d' | awk '!a[$4]++' | sed 's/[ \t]\+/\t/g' > all_gblast_uniq_hits.bed

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#IGV snapshot

cd /media/aswin/saumya/GPRC6A/Comment/transciption_start_sites/CAGE_TSS/genome
#create batch file
cat Bos_taurus_ENST00000310357_exon_wise/all_gblast_uniq_hits.bed | awk 'NR==1{print$1,$2-2500} END{print$3+2500,"GPRC6A_genic_region"}' OFS="\t" | paste -s -d "\t" > GPRC6A_genic_region.bed
/media/aswin/programs/IGV-snapshot-automator-20.11.1/make_IGV_snapshots.py ../tissues_TSS/TSS_adrenal_gland_cortex_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_cerebellum_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_cerebrum_cortex_countAbove10_bb_minusY.bed \
	../tissues_TSS/TSS_colon_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_duodenum_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_heart_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_hypothalamus_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_ileum_countAbove10_bb_minusY.bed \
	../tissues_TSS/TSS_jejunum_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_kidney_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_liver_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_lung_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_lymph_node_countAbove10_bb_minusY.bed \
	../tissues_TSS/TSS_mammary_gland_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_ovary_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_pancreas_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_pituitary_gland_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_rumen_countAbove10_bb_minusY.bed \
	../tissues_TSS/TSS_skeletal_muscle_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_spleen_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_subcutaneous_fat_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_testis_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_thoracic_branch_lymph_countAbove10_bb_minusY.bed \
	../tissues_TSS/TSS_thyroid_gland_countAbove10_bb_minusY.bed ../tissues_TSS/TSS_uterine_horn_countAbove10_bb_minusY.bed \
 Bos_taurus_ENST00000310357_exon_wise/all_gblast_uniq_hits.bed -g ARS-UCD1.2_Btau5.0.1Y.fa -r GPRC6A_genic_region.bed -ht 2000 -o ./ -nosnap -suffix "GPRC6A_Multi_tissue_CAGE_TSS"
#Run batch script
sed '/^snapshot/ s/.png/.svg/g' IGV_snapshots.bat -i
/media/aswin/programs/IGV_Linux_2.19.1/igv_hidpi.sh -b IGV_snapshots.bat

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#IGV reports

# 1. Create a temporary folder for the stripped files
mkdir -p ./tmp_cleaned_bed

# 2. Fix the naming mismatch for ALL bed files matching your wildcard
#for f in ../tissues_TSS/*.bed; do sed 's/^chr//' "$f" > "./tmp_cleaned_bed/$(basename "$f")"; done
#for f in ../tissues_TSS_Enhancers/*.bed; do sed 's/^chr//' "$f" > "./tmp_cleaned_bed/$(basename "$f")"; done

#rename chromosomes to match genome seq IDs & label tissue names inside bed file
for f in ../tissues_TSS/*.bed
do
id=$(echo $f | cut -f3 -d "/" | cut -f2- -d "_" | sed 's/_countAbove10_bb_minusY.bed//g')
sed 's/^chr//' "$f" | awk -v id="$id" '$4=id' OFS="\t"> "./tmp_cleaned_bed/$(basename "$f")"
unset id
done

for f in ../tissues_TSS_Enhancers/*.bed
do
id=$(echo $f | cut -f3 -d "/" | cut -f2- -d "_" | sed 's/_countAbove10_bb_minusY.bed//g')
sed 's/^chr//' "$f" | awk -v id="$id" '$4=id' OFS="\t"> "./tmp_cleaned_bed/$(basename "$f")"
unset id
done

#Combine all tissue data
cat tmp_cleaned_bed/TSS_*.bed > all_TSS_countAbove10_bb_minusY.bed
cat tmp_cleaned_bed/BC_*.bed > all_BS_countAbove10_bb_minusY.bed

#Subset data to GPRC6A region
awk -v chrom="9" -v start=33772607 -v end=33842561 '($1 == chrom) && $2 >= start && $3 <= end' all_TSS_countAbove10_bb_minusY.bed > CAGE_seq_TSS.bed
awk -v chrom="9" -v start=33772607 -v end=33842561 '($1 == chrom) && $2 >= start && $3 <= end' all_BS_countAbove10_bb_minusY.bed > CAGE_seq_Enhancers.bed
#Strand-wise TSS
awk '$NF=="+"' CAGE_seq_TSS.bed > plus_CAGE_seq_TSS.bed 
awk '$NF=="-"' CAGE_seq_TSS.bed > minus_CAGE_seq_TSS.bed 

# 3. Execute create_report using the newly cleaned tracks folder
create_report "$BED1" \
  --fasta "$GENOME" \
  --sequence 1 \
  --begin 2 \
  --end 3 \
  --zero_based true \
  --flanking 2000 \
  --output SRR9136449_GPRC6A_region.html \
  --title "GPRC6A Region (+/- 2000bp Flanking) RNA-Seq" \
  --tracks "$BED2" ./tmp_cleaned_bed/*.bed

cp ../genome/Bos_taurus_ENST00000310357_exon_wise/all_gblast_uniq_hits.bed APOBEC1_exons.bed
  
GENOME="ARS-UCD1.2_Btau5.0.1Y.fa"
BED1="GPRC6A_genic_region.bed"
BED2="APOBEC1_exons.bed"

# 2. Command execution
time create_report "$BED1" --fasta "$GENOME" --sequence 1 --begin 2 --end 3 --zero_based true --flanking 2000 --output CAGE_TSS_Enhancers_Bos_taurus_GPRC6A.html --title "GPRC6A Region (+/- 2000bp Flanking) CAGE TSS Enhancers" --tracks "$BED2" ./tmp_cleaned_bed/*.bed
#create_report "$BED1" --fasta "$GENOME" --sequence 1 --begin 2 --end 3 --zero_based true --flanking 2000 --output test.html --title "GPRC6A Region (+/- 2000bp Flanking) RNA-Seq" --tracks "$BED2" corrected_adrenal_gland.bed

###############################################################################################################################################################################################################################################################################################################
#View all TSS together

#As coordinates remain exactly same between 2 genomes 
#Combine NCBI TSS,Promoter 2.0 TSS & CAGE sequencing TSS frompaper

#Subset TSS data to only GPRC6A region
cd /media/aswin/saumya/GPRC6A/Comment/transciption_start_sites/CAGE_TSS/genome/tmp_cleaned_bed
mkdir -p subsetted

for f in *.bed; do
    if [[ "$f" == subsetted/* ]]; then continue; fi
    base=$(basename "$f" .bed)
    awk -v chrom="9" -v start=33772607 -v end=33842561 '
        ($1 == chrom || $1 == "chr"chrom) && $2 >= start && $3 <= end
    ' "$f" > "subsetted/${base}_subset.bed"
done

###############################################################################################################################################################################################################################################################################################################
#Get Bos taurus GPRC6A TOGA output from github

cd /media/aswin/Saumya/GPRC6A/Comment
git clone --depth 1 --filter=blob:none --sparse https://github.com/saumyagupta09/GPRC6A
cd GPRC6A
git sparse-checkout set TOGA/Bos_taurus_TOGA
mv TOGA/ ../
rm -r GPRC6A

###############################################################################################################################################################################################################################################################################################################
#Identify slippery sequence:

cd /media/aswin/Saumya/GPRC6A/Comment/slippery_sequence
cp ../fasta/Bos_taurus_ENST00000310357_exon_wise_with_repeat.fa .
egrep "AAAAAAA|AAAAAAC|AAAAAAG|AAAAAAT|AAATTTA|AAATTTT|CCCTTTA|CCCTTTT|CCCAAAA|GGGAAAC|GGGAAAT|GGGAAAG|GGGAAAA|GGGTTTA|GGGTTTT|TTTTTTA|TTTAAAC|TTTTTTG|TTTAAAT|and|TTTAAAG" -iB1 Bos_taurus_ENST00000310357_exon_wise_with_repeat.fa

#Tools to identify PRF:
#PKScan: https://www.biorxiv.org/content/10.1101/2021.04.23.441185v1.full
#https://pmc.ncbi.nlm.nih.gov/articles/PMC3705612/



###############################################################################################################################################################################################################################################################################################################
#Check ORFs

#Run ORFFinder in NCBI
#Query: Bos_taurus_ENST00000310357_cds_with_repeat.fa

cd /media/aswin/saumya/GPRC6A/Comment/orfs
makeblastdb -in Bos_taurus_ENST00000310357_cds_with_repeat.fa -out Bos_taurus_ENST00000310357_cds_with_repeat.fa -dbtype nucl
blat Bos_taurus_ENST00000310357_cds_with_repeat.fa Bos_taurus_ENST00000310357_exon_wise_with_repeat.fa > output.psl
cat output.psl  | awk 'NR>5 {print $14,$16,$17,$10,"1",$9}' OFS="\t" > GPRC6A_exons.bed



