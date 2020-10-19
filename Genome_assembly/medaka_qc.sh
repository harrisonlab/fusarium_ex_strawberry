# Medaka
# Run in medaka env
#A tool to create a consensus sequence from nanopore sequencing data.
# This task is performed using neural networks applied from a pileup of individual sequencing reads against a draft assembly.
# It outperforms graph-based methods operating on basecalled data, and can be competitive with state-of-the-art signal-based methods, whilst being much faster.

# Rename contigs for genome
# If split or remove contigs is needed, provide FCSreport file by NCBI.

# ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
#    touch tmp.txt
#    for Assembly in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/flye_raw/racon_10/assembly_racon_round_4.fasta); do
#        OutDir=$(dirname $Assembly)
#        $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/assembly_r4_renamed_FR.fasta --coord_file tmp.txt > $OutDir/log.txt
#    done
#    rm tmp.txt

for Assembly in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/flye_raw/racon_10/assembly_r4_renamed_FR.fasta); do
  ReadsFq=$(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz)
  OutDir=$(dirname $Assembly)/medaka
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
  sbatch $ProgDir/medaka.sh $Assembly $ReadsFq $OutDir
done
