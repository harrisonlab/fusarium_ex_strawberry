# System for rapidly aligning entire genomes
# compare Fofr14 against cepae genome and lycopersici genomes
#SubjectGenome=$(basename $1)
#QueryGenome=$(basename $2)
#Prefix=$3
#OutDir=$4


for SubjectGenome in (ls assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta); do
  Organism=$(echo "$Assembly" | rev | cut -d '/' -f4 | rev)
  Strain=$(echo "$Assembly" | rev | cut -d '/' -f3 | rev)
  echo "$Organism - $Strain"
  Query=path/to/lycopersici/genome
  Prefix=FoFr_14
  OutDir=alignment/miniasm
  mkdir -p $OutDir
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
  sbatch $ProgDir/mummer.sh $SubjectGenome $QueryGenome $Prefix $OutDir
done
