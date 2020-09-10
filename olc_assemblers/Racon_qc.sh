# Racon
# Need to do for Miniasm, Flye and SMARTdenovo output files

for Assembly in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/assembly.fasta); do
    ReadsFq=$(ls raw_dna/FAL69458.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
