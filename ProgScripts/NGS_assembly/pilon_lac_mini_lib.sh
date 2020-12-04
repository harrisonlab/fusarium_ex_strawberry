#!/usr/bin/env bash
#SBATCH -J pilon
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=30

# Align raw reads to a pacbio assembly and then use this alignmeant to correct
# indels and substitutions in the assembly.

Mem="470G"
#Mem="70G"
# Mem="372G"
Threads=8

# ---------------
# Step 1
# Collect inputs
# ---------------

Assembly=$(basename $1)
Read_F=$(basename $2)
Read_R=$(basename $3)
OutDir=$4
Iterations=$5
if [ $6 ]; then
  Ploidy=$6
else
  Ploidy="haploid"
fi

CurDir=$PWD
echo  "Running Pilon with the following inputs:"
echo "Pacbio assembly - $Assembly"
echo "Forward trimmed reads - $Read_F"
echo "Reverse trimmed reads - $Read_R"
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
cp $CurDir/$1 FoLR1_conc_racon_round_6_renamed.fasta
cp $CurDir/$2 $Read_F
cp $CurDir/$3 $Read_R

mkdir best_assembly
cp FoLR1_conc_racon_round_6_renamed.fasta best_assembly/.

for i in $(seq 1 $Iterations); do
  echo "Running Iteration: $i"
  mkdir $WorkDir/"correction_$i"
  cd $WorkDir/correction_$i
  cp $WorkDir/best_assembly/FoLR1_conc_racon_round_6_renamed.fasta .

  # ---------------
  # Step 3.a
  # Align seq reads
  # ---------------
  # Prepare the assembly for alignment
  # Align reads against the assembly
  # Convert the SAM file to BAM in preparation for sorting.
  # Sort the BAM file, in preparation for SNP calling:
  # Index the bam file


  bowtie2-build FoLR1_conc_racon_round_6_renamed.fasta FoLR1_conc_racon_round_6_renamed.fasta.indexed
  bowtie2 -p 12 -x FoLR1_conc_racon_round_6_renamed.fasta.indexed -1 $WorkDir/$Read_F -2 $WorkDir/$Read_R  -S FoLR1_conc_racon_round_6_renamed.fasta_aligned.sam
  samtools view --threads 20 -bS FoLR1_conc_racon_round_6_renamed.fasta_aligned.sam -o FoLR1_conc_racon_round_6_renamed.fasta_aligned.bam
  samtools sort --threads 20 FoLR1_conc_racon_round_6_renamed.fasta_aligned.bam -o FoLR1_conc_racon_round_6_renamed.fasta_aligned_sorted.bam
  samtools index FoLR1_conc_racon_round_6_renamed.fasta_aligned_sorted.bam

  # ---------------
  # Step 3.b
  # Run Pilon
  # ---------------
  # Run pilon to polish
  if [ $Ploidy == "haploid" ]; then
    JavaDir=/scratch/software/pilon-1.23
    java -Xmx$Mem -jar $JavaDir/pilon-1.23.jar --threads 16 --genome FoLR1_conc_racon_round_6_renamed.fasta --changes --frags FoLR1_conc_racon_round_6_renamed.fasta_aligned_sorted.bam --outdir .
  elif [ $Ploidy == "diploid" ]; then
     JavaDir=/scratch/software/pilon-1.23
    java -Xmx$Mem -jar $JavaDir/pilon-1.23.jar --threads 16 --genome FoLR1_conc_racon_round_6_renamed.fasta --changes --diploid --frags FoLR1_conc_racon_round_6_renamed.fasta_aligned_sorted.bam --outdir .
  else
  echo "ploidy not recognised"
  fi
  cp pilon.fasta $WorkDir/best_assembly/FoLR1_conc_racon_round_6_renamed.fasta
  # cp pilon.changes $WorkDir/best_assembly/pilon_$i.changes
  cp pilon.fasta $CurDir/$OutDir/pilon_$i.fasta
  cp pilon.changes $CurDir/$OutDir/pilon_$i.changes
  cd $WorkDir
done

mv $WorkDir/best_assembly/FoLR1_conc_racon_round_6_renamed.fasta $WorkDir/best_assembly/pilon.fasta

# mkdir -p $CurDir/$OutDir
# cp $WorkDir/best_assembly/* $CurDir/$OutDir/.
rm -r $WorkDir
