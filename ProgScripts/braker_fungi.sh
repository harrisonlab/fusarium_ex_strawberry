#!/usr/bin/env bash
#SBATCH -J braker
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=16


WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
Assembly=$1
OutDir=$2
AcceptedHits=$3
GeneModelName=$4
CurDir=$PWD

echo "$Assembly"
echo "$OutDir"
echo "$AcceptedHits"
echo "$GeneModelName"
echo "$CurDir"

mkdir -p $WorkDir
cd $WorkDir

cp $CurDir/$Assembly assembly.fa
cp $CurDir/$AcceptedHits alignedRNA.bam



braker.pl \
  --GENEMARK_PATH=/home/akinya/miniconda3/envs/Repenv/opt/genemark_es/gmes_petap \
  --BAMTOOLS_PATH=/home/akinya/miniconda3/envs/Repenv/bin \
  --overwrite \
  --fungus \
  --gff3 \
  --softmasking on \
  --species=$GeneModelName \
  --genome="assembly.fa" \
  --bam="alignedRNA.bam"

mkdir -p $CurDir/$OutDir
cp -r braker/* $CurDir/$OutDir/.

rm -r $WorkDir
