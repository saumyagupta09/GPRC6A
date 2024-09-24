#!/bin/bash
# This script generates GC content plots for Caviomorpha species in 100 bp windows.
# 
# Usage:
#   ./gc_content_windows.sh <species_name> <output_file> <genome_file>
#
# Parameters:
#   species_name: The name of the species (e.g., Chinchilla_lanigera).
#   output_file: The blast output file in outfmt 6 (e.g., Heterocephalus_glaber_on_Chinchilla_lanigera_outfmt6_final.out).
#   genome_file: The genome file in FASTA format (e.g., genome.fna).


sps=$1 #give the name of the species with an underscore eg: Chinchilla_lanigera
blast_file=$2 #blast file needs to be in a specific outfmt format. 
genome=$3 #genome file


#1. GET CORRECT EXON POSITION FROM BLAST OUTFMT 6
while read i
 do
     col2=`echo $i | awk '{print $2}'`
     col1=`echo $i | awk '{print $1}'`
     col3=`echo $i | awk '{print $3}'`
     col4=`echo $i | awk '{print $4}'`
     col5=`echo $i | awk '{print $5}'`
     col6=`echo $i | awk '{print $6}'`
     col7=`echo $i | awk '{print $7}'`
     col8=`echo $i | awk '{print $8}'`
     col9=`echo $i | awk '{print $9}'`
     col10=`echo $i | awk '{print $10}'`
     col12=`echo $i | awk '{print $12}'`
 
     if [ $col3 != $col4 ]
     then
         x=`expr $col5 - 1` 
         y=`expr $col3 - $col6` 
         if [ $col9 == plus ]
         then
             ncol7=`expr $col7 - $x`
             ncol8=`expr $col8 + $y`
             echo -e "$col2\t$ncol7\t$ncol8\t$col1\t$col9\t"+"" >> "$sps"_actual_exons.bed
 else
             ncol8=`expr $col7 + $x`  
             ncol7=`expr $col8 - $y`
             echo -e "$col2\t$ncol8\t$ncol7\t$col1\t$col9\t"-"" >> "$sps"_actual_exons.bed
 fi
 
 else 
 echo "$i" | awk -v col2="$col2" -v col1="$col1" -v col7="$col7" -v col8="$col8" -v col9="$col9" \
            '{ if (col9 == "plus") { print col2 "\t" col7 "\t" col8 "\t" col1 "\t" col9 "\t+" }
               else if (col9 == "minus") { print col2 "\t" col7 "\t" col8 "\t" col1 "\t" col9 "\t-" } }' >> "$sps"_actual_exons.bed
 
 fi
 done < $blast_file

#2. ADD INTRON AND EXON INFORMATION

awk '{
    if ($5 == "minus") {
        if ($2 < $3) {
            temp = $2
            $2 = $3
            $3 = temp
        }
    }
    exons[NR] = $0
    if (NR > 1) {
        intron_start = prev[3]
        intron_end = $2 
        if (prev[5] == "minus") {
            if (intron_start > intron_end) {
                temp = intron_start
                intron_start = intron_end
                intron_end = temp
            }
        }
        print prev[1], intron_end, intron_start, "intron_" substr(prev[4], index(prev[4], "_") + 1), prev[5], prev[6]
    }
    split($0, prev)

    print $1, $2, $3, $4, $5, $6
}' OFS='\t' "$sps"_actual_exons.bed > "$sps"_intron_boundary.bed



#3. ADD UPSTREAM AND DOWNSTREAM BASEPAIRS ACCORDINGLY (the division of the absent exon-intron position in the plot is exon 1=0-1kb; exon2=4-5kb ; exon3=7kb position; exon 6= last 2kb bp)
file="$sps"_intron_boundary.bed
output_file="${sps}_adjusted.bed"

first_exon=""
last_exon=""
first_chr=""
first_start=""
last_end=""

while read -r chr start end exon strand rest extra; do
    if [ -z "$first_exon" ]; then
        first_exon=$exon
        first_chr=$chr
        first_start=$start
    fi
    last_exon=$exon
    last_end=$end
    echo -e "$chr\t$start\t$end\t$exon\t$strand\t$rest\t$extra"
done < "$file" > "$output_file"


if [[ $first_exon == "exon_1" ]]; then #recognizes there is no missing exon at the start.
    new_start=$((first_start + 1000)) #adds 1000 bp upstream from exon 1 start position.
    echo -e "$first_chr\t$new_start\t$first_start\tupstream1kb\tminus\t-" | cat - "$output_file" > temp && mv temp "$output_file"
elif [[ $first_exon == "exon_2" ]]; then #recognizes the missing exon is exon 1 
    new_start=$((first_start + 4000)) #adds 4000 bp upstream from exon 2 start position.
    echo -e "$first_chr\t$new_start\t$first_start\texon_1_upstream1kb\tminus\t-" | cat - "$output_file" > temp && mv temp "$output_file"
elif [[ $first_exon == "exon_3" ]]; then #recognizes the missing exons are exon 1 and exon 2
    new_start=$((first_start + 7000)) #adds 7000 bp upstream from exon 2 start position.
    echo -e "$first_chr\t$new_start\t$first_start\texon_1_2_upstream1kb\tminus\t-" | cat - "$output_file" > temp && mv temp "$output_file"
fi

if [[ $last_exon == "exon_5" ]]; then #recognizes the missing exons is exon 6
    new_end=$((last_end - 4000)) #adds 4000 bp downstream from exon 5 end position.
    echo -e "$first_chr\t$last_end\t$new_end\texon_6_downstream1kb\tminus\t-" >> "$output_file"
elif [[ $last_exon == "exon_6" ]]; then #recognizes there is no missing exon at the end 
    new_end=$((last_end - 1000)) #adds 1000 bp downstream from exon 6 end position.
    echo -e "$first_chr\t$last_end\t$new_end\tdownstream1kb\tminus\t-" >> "$output_file"
fi            ##### the results are saved in "${sps}_adjusted.bed"

sed 's/\t*$//' "${sps}_adjusted.bed" -i #removes the tabs at the end of the file

### CALCULATE THE LENGTH OF EXON AND INTRON
awk '{diff = ($2 > $3) ? $2 - $3 : $3 - $2; print $0, diff}' OFS='\t' "${sps}_adjusted.bed" > "${sps}_adjusted_length.bed"


#### MAKE WINDOWS
input_file="${sps}_adjusted_length.bed" #making bedfile for extracting the genic region with coordinates
output_file="${sps}_1kb_flank.bed"

first_column=$(head -n 1 "$input_file" | cut -f1)
min_value=$(awk '{print $2; print $3}' "$input_file" | sort -n | head -n 1)
max_value=$(awk '{print $2; print $3}' "$input_file" | sort -n | tail -n 1)
difference=$((max_value - min_value))
echo -e "$first_column\t$min_value\t$max_value\tGPRC6A_1kb_flank\tminus\t-\t$difference" > "$output_file"


bedtools getfasta -fi $genome -bed "$sps"_1kb_flank.bed -s > ${sps}_1kb_flank.fa

sed 's/:.*//g' ${sps}_1kb_flank.fa -i #correcting the header in fasta file 

samtools faidx ${sps}_1kb_flank.fa #indexes the fasta file for future use

awk '{print$1,$NF}' OFS='\t' "${sps}_1kb_flank.bed" > "${sps}_GPRC6A_gene_size.txt" #saves the length of the gene for makewindows to divide it

bedtools makewindows -g "${sps}_GPRC6A_gene_size.txt" -w 100 > "$sps"_windows.txt #makes windows in 100bp windows

bedtools nuc -fi ${sps}_1kb_flank.fa -bed "$sps"_windows.txt > "$sps"_gccontent.bed #calculates gc and at content in each window

sed '1d' "$sps"_gccontent.bed -i #remove header to avoid complications in plot

awk '{print$1,$2,$3,$5}' OFS='\t' "$sps"_gccontent.bed > "$sps"_gcpcnt.bed # Creates a file showing only gc content

#making a file with lengths of exon 1 and intron in absolute terms to plot using rect in R
start=0
while read line
do
name=`echo $line | awk '{print$4}'`
len=`echo $line | awk '{print$7}'`
end=$(expr $start + $len)

echo -e "$name\t$start\t$end\t$len"

start=$end
done < "${sps}_adjusted_length.bed" > "${sps}_absolute_lengths.txt"  #make sure start becomes 0 again after running this script

unset start
