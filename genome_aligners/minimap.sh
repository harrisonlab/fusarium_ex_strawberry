#!/usr/bin/env bash
#SBATCH -J minimap
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30




##########################################################################
#INPUT:
# 1st argument: Genome1
# 2nd argument: Genome2
#OUTPUT:
# Commands to make dotplot in R

Genome1=$1
Genome2=$2
Strain=$3
outdir=$4


CurDir=$PWD

WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir

cp $Genome1 $WorkDir
cp $Genome2 $WorkDir
cd $WorkDir

query=$Genome1
target=$Genome2

minimap2=/home/akinya/miniconda3/envs/olc_assemblers/bin/minimap2

$minimap2 -x asm5 -t 36 $target $query > "$Strain"_minimap.paf

cp $WorkDir/"$Strain"_minimap.paf $outdir
rm -r $WorkDir


