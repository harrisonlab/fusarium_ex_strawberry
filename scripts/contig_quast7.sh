#!/bin/bash

# quast contig fix step

ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamedA.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
    sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
  done
