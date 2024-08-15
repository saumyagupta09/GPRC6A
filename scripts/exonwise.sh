#!/bin/bash

# $1 is the input file, which should be the output from running the following BLAST command:
# time blastn -task blastn -query query.fa -db genome.fa -evalue 0.001 -outfmt "6 qseqid sseqid qlen length qstart qend sstart send sstrand evalue bitscore sseq" -out example_file.out
# $2 is the reference genome file.

# An example input file, named 'example_file.out', is provided for this script. Please refer to it to understand the expected input format.

file=$1
genome=$2

# Reading each line from the input file
while read i
do
    # Extract specific columns from the input line
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

    # Check if column 3 is not equal to column 4
    if [ $col3 != $col4 ]
    then
        # Calculate new values for further processing
        x=`expr $col5 - 1`
        y=`expr $col3 - $col6`

        # If the strand is 'plus'
        if [ $col9 == plus ]
        then
            # Adjust coordinates based on the calculations
            ncol7=`expr $col7 - $x - 1`
            ncol8=`expr $col8 + $y`

            # Create a BED file with the updated coordinates
            echo -e "$col2\t$ncol7\t$ncol8\t$col1" > "$col1".bed

            # Extract the sequence from the genome using BEDTools
            col12=`bedtools getfasta -fi $genome -bed "$col1".bed | grep -v '>'`

        else
            # Adjust coordinates for the reverse strand
            ncol8=`expr $col7 + $x`
            ncol7=`expr $col8 - $y - 1`

            # Create a BED file with the updated coordinates
            echo -e "$col2\t$ncol7\t$ncol8\t$col1" > "$col1".bed

            # Extract the sequence from the genome using BEDTools and save it to a .fa file
            bedtools getfasta -fi $genome -bed "$col1".bed | grep -v '>' > "$col1".fa

            # Generate the reverse complement of the sequence
            revseq -sequence "$col1".fa -outseq outseq_"$col1".fa

            # Extract the reverse complement sequence
            col12=`grep -v ">" outseq_"$col1".fa`
        fi
    fi

    # Print the sequence identifier and the sequence itself
    echo ">"$col1""
    echo "$col12"

# Continue reading lines from the input file
done < $file
