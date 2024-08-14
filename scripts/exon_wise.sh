#!/bin/bash
file=$1
genome=$2

while read i
  do   col2=`echo $i|awk '{print$2}'`
      col1=`echo $i|awk '{print$1}'`
   col3=`echo $i|awk '{print$3}'`
         col4=`echo $i|awk '{print$4}'`
   col5=`echo $i|awk '{print$5}'`
   col6=`echo $i|awk '{print$6}'`
   col7=`echo $i|awk '{print$7}'`
   col8=`echo $i|awk '{print$8}'`
  col9=`echo $i|awk '{print$9}'`
  col10=`echo $i|awk '{print$10}'`
  col12=`echo $i|awk '{print$12}'`
          if [ $col3 != $col4 ]
 then x=`expr $col5 - 1`
 y=`expr $col3 - $col6`
           if [ $col9 == plus ]
  then   ncol7=`expr $col7 - $x - 1`
 ncol8=`expr $col8 + $y`
            echo -e "$col2\t$ncol7\t$ncol8\t$col1" > "$col1".bed
 col12=`bedtools getfasta -fi $genome -bed "$col1".bed | grep -v '>'` 
             else ncol8=`expr $col7 + $x`
   ncol7=`expr $col8 - $y - 1`
     echo -e "$col2\t$ncol7\t$ncol8\t$col1" > "$col1".bed
      bedtools getfasta -fi $genome -bed "$col1".bed | grep -v '>' > "$col1".fa
 revseq -sequence "$col1".fa -outseq outseq_"$col1".fa
 col12=`grep -v ">" outseq_"$col1".fa`
 fi
 fi
        echo ">"$col1""
  echo "$col12"
 done < $file 
