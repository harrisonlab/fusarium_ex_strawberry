# Medaka
# Run in medaka env
#A tool to create a consensus sequence from nanopore sequencing data.
# This task is performed using neural networks applied from a pileup of individual sequencing reads against a draft assembly.
# It outperforms graph-based methods operating on basecalled data, and can be competitive with state-of-the-art signal-based methods, whilst being much faster.

for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/racon_10/*_renamed.fasta); do
  ReadsFq=$(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz)
  OutDir=$(dirname $Assembly)/medaka
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
  sbatch $ProgDir/medaka.sh $Assembly $ReadsFq $OutDir
done
