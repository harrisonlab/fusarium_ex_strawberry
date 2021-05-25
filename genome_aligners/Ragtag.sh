#!/usr/bin/env bash
#SBATCH -J ragtag
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=10
##########################################################################
# Run ragtag for genome correction
#INPUT:
# 1st argument: refrence genome
# 2nd argument: query genome to be corrected
# 3rd argument: output directory
#OUTPUT:
# corrected genome
reference=$1
query=$2
outdir=$3
CurDir=$PWD
WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cp $reference $WorkDir
cp $query $WorkDir
cd $WorkDir
ragtag=/home/bourquif/miniconda2/pkgs/ragtag-1.1.1-pyh7b7c402_0/python-scripts/ragtag.py
#$ragtag correct \
#-u \
#$reference \
#$query
$ragtag	scaffold \
-u -C \
$reference \
$query
#ragtag_output/*corrected.fasta
#$ragtag	merge \
#-u -C -b hic.bam \
#$query \
#ragtag_output/*.agp
cp ragtag_output/* $outdir
#cp ragtag_output/*.corrected.fasta $outdir
rm -r $WorkDir
