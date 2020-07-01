Read.md

# Assembling DNA utilising Pipelines developed by Andy Armitage - slightly modded to adapt to NewSlurm from old GRDENG

# Spades assembly step
# Assembles illumina sequenced DNA 
for StrainPath in $(ls -d qc_dna/paired/*/*)
do
    ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/spades
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    OutDir=assembly/spades/$Organism/${Strain}
    echo $F_Read
    echo $R_Read
    sbatch $ProgDir/slurm_spades.sh $F_Read $R_Read $OutDir correct 10
done

# Quast assembly QC step
ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast

for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta)
do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  # OutDir=assembly/spades/$Organism/$Strain
  OutDir=$(dirname $Assembly)

  sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
done
