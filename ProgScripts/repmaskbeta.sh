#!/bin/bash

#executed through Antonios directory in the 1st run
#copied the scripts to my directory after they were fixed

for Assembly in $(ls assembly/spades/*/DSA14_003/ncbi_edits/contigs_fuox14.fasta | grep -v '_2' | grep -v '11055'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism"
echo "$Strain"
OutDir=repeat_masked/$Organism/$Strain/ncbi_edits_repmask
ProgDir=/home/akinya/git_repos/tools/seq_tools/repeat_masking
sbatch $ProgDir/rep_modelingBeta.sh $Assembly $OutDir
done
