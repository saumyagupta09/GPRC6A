#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load HYPHYMP
wd=`pwd`
cd $1
workdir=`pwd`
#Create a config file to run HYPHY
echo -ne "1\n6\n"$workdir"/Artocarpus_hirsutus_tree_labelled\n" > $1_absrel.config
#Run HYPHYMP
HYPHYMP < $1_absrel.config > "$workdir"/results_"$1"/results.absrel.txt
cp *.ABSREL.json "$workdir"/results_"$1"/
cd $wd
