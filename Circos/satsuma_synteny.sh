#!/usr/bin/env bash
#SBATCH -J satsumasynteny
#SBATCH --partition=long
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30

##########################################################################
#INPUT:
# 1st argument: Genome1
# 2nd argument: Genome2
#OUTPUT:
# File is synteny alignment between two input genomes

Genome1=$1
Genome2=$2
outdir=$3


CurDir=$PWD

WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir

cp $Genome1 $WorkDir
cp $Genome2 $WorkDir
cd $WorkDir


/home/akinya/SatsumaSynteny/satsuma-code-0/SatsumaSynteny -t $Genome1 -q $Genome2 -o $WorkDir


cp $WorkDir/satsuma_summary.chained.out $outdir
rm -r $WorkDir
