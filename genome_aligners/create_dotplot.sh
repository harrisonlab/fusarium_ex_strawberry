#!/usr/bin/env bash
#SBATCH -J dotplot
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30




##########################################################################
#INPUT:
# 1st argument: paf file
#OUTPUT:
# dorplot

paf=$1
Strain=$(echo $paf_file | rev | cut -f2 -d '/' | rev)
outdir=$2


CurDir=$PWD

WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir

cp $paf $WorkDir
cd $WorkDir


dotplot=/home/akinya/git_repos/fusarium_ex_strawberry/genome_aligners
$dotplot/pafCoordsDotPlotly.R \
-i $paf \
-o paf_plot \
-m 100 \
-q 500 \
-k 10 \
-s -t -l -p 12




cp $WorkDir/paf_plot.png $outdir
rm -r $WorkDir
