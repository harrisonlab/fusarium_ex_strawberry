#!/usr/bin/env bash
#SBATCH -J Circos
#SBATCH --partition=short
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=10

##########################################################################
#INPUT:
# 1st argument: Configuration_file
#OUTPUT:
# Ideogram

Configuration_file=$1
outdir=$2

CurDir=$PWD

WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir

cp $Genome1 $WorkDir
cp $Genome2 $WorkDir
cd $WorkDir


circos=/home/connellj/miniconda2/bin/circos
$circos \
-conf $Configuration_file \
-outputdir $WorkDir



cp $WorkDir/circos.png $outdir
#cp $WorkDir/circos.svg $outdir
rm -r $WorkDir
