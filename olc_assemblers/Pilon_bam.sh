#!/usr/bin/env bash
#SBATCH -J pilon
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=30

# Align raw reads to a pacbio assembly and then use this alignmeant to correct
# indels and substitutions in the assembly.

# ---------------
# Step 1
# Collect inputs
# ---------------

Assembly=$(basename $1)
BAM=$(basename $2)
OutDir=$3
Iterations=$4
if [ $5 ]; then
  Ploidy=$5
else
  Ploidy="haploid"
fi

#Assembly=$(basename $1)
#Read_F=$(basename $2)
#Read_R=$(basename $3)
#OutDir=$4
#Iterations=$5
#if [ $6 ]; then
  #Ploidy=$6
#else
  #Ploidy="haploid"
#fi

CurDir=$PWD
echo  "Running Pilon with the following inputs:"
echo "Pacbio assembly - $Assembly"
echo "Alignment file - $BAM"
#echo "Forward trimmed reads - $Read_F"
#echo "Reverse trimmed reads - $Read_R"
echo "OutDir - $OutDir"
echo "Running Pilon the following number of times - $Iterations"
echo "Ploidy set to: $Ploidy"

mkdir -p $CurDir/$OutDir

# ---------------
# Step 2
# Copy data
# ---------------

WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir
cp $CurDir/$1 race_1_smartdenovo_racon_round_10_renamed.fasta
#cp $CurDir/$2 $Read_F
#cp $CurDir/$3 $Read_R
cp $CurDir/$2 bwa_mapping.SDEN.sorted.bam

mkdir best_assembly
cp race_1_smartdenovo_racon_round_10_renamed.fasta best_assembly/.

for i in $(seq 1 $Iterations); do
  echo "Running Iteration: $i"
  mkdir $WorkDir/"correction_$i"
  cd $WorkDir/correction_$i
  cp $WorkDir/best_assembly/race_1_smartdenovo_racon_round_10_renamed.fasta .


  # ---------------
  # Step 3.b
  # Run Pilon
  # ---------------
  # Run pilon to polish
  if [ $Ploidy == "haploid" ]; then
    pilon --threads 16 --genome race_1_smartdenovo_racon_round_10_renamed.fasta --changes --bam $BAM --outdir .
    #pilon --threads $Threads --genome assembly.fa --changes --frags assembly.fa_aligned_sorted.bam --outdir .
  elif [ $Ploidy == "diploid" ]; then
    pilon --threads 16 --genome race_1_smartdenovo_racon_round_10_renamed.fasta--changes --diploid --bam $BAM --outdir .
    #pilon --threads $Threads --genome assembly.fa --changes --diploid --frags assembly.fa_aligned_sorted.bam --outdir .
  else
    echo "ploidy not recognised"
  fi
  cp pilon.fasta $WorkDir/best_assembly/race_1_smartdenovo_racon_round_10_renamed.fasta
  # cp pilon.changes $WorkDir/best_assembly/pilon_$i.changes
  cp pilon.fasta $CurDir/$OutDir/pilon_$i.fasta
  cp pilon.changes $CurDir/$OutDir/pilon_$i.changes
  cd $WorkDir
done

mv $WorkDir/best_assembly/race_1_smartdenovo_racon_round_10_renamed.fasta $WorkDir/best_assembly/pilon.fasta

# mkdir -p $CurDir/$OutDir
# cp $WorkDir/best_assembly/* $CurDir/$OutDir/.
rm -r $WorkDir
