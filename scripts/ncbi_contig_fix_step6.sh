#!/bin/bash

for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
NCBI_report=$(ls genome_submission/$Organism/$Strain/initial_sub/Contamination*.txt)
OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
mkdir -p $OutDir
ProgDir=~/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants2.py --keep_mitochondria --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
done
