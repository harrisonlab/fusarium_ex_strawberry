#!/bin/bash

# Develop result files  from quast analysis
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/transposed_report.tsv); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    # echo;
    # echo $Organism;
    # echo $Strain;
    cat $Assembly | tail -n +2 | sed "s/contigs_min_500bp/${Organism}_${Strain}/g"
  done > assembly/quast_results.txt
