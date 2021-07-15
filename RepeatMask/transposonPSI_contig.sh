#!/usr/bin/env bash
#SBATCH -J transposonPSI
#SBATCH --partition=short
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=60

Usage="sbatch transposonPSI.sh <assembled_contigs.fa> [<output_directory>]"

InFile=$1

Organism=F.oxysporum_fsp_cepae
Strain=Fus2_canu_new
contig=contig_21

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

if [ $2 ]; then
  OutDir=$CurPath/$2
else
  OutDir=$CurPath/repeat_masked/$Organism/$Strain/"$contig"_repmask
fi

mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$InFile "$contig"_unmasked.fa

# The full path is needed
/home/gomeza/miniconda3/envs/general_tools/share/transposonPSI/transposonPSI.pl "$contig"_unmasked.fa nuc

mkdir -p $OutDir
cp -r $WorkDir/* $OutDir/.

rm -r $WorkDir
