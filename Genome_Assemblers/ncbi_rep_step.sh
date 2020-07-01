#!/bin/bash 

for Assembly in $(ls assembly/spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -v 'ncbi_edits'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    NCBI_report_dir=genome_submission/$Organism/$Strain/initial_submission
    mkdir -p $NCBI_report_dir
done
