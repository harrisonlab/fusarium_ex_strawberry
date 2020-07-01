#!/bin/bash

# quast assembly QC step
ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast

for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta)
do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  # OutDir=assembly/spades/$Organism/$Strain
  OutDir=$(dirname $Assembly)

  sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
done
