#!/usr/bin/env bash
#SBATCH -J repeatmasker
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=24

#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking
#BestAssembly=assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta
#OutDir=repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/SMARTdenovo/ncbi_edits_repmask
#sbatch $ProgDir/rep_modeling.sh $BestAssembly $OutDir

InFile=$1
Organism=$(echo $InFile | rev | cut -d "/" -f4 | rev)
Strain=$(echo $InFile | rev | cut -d "/" -f3 | rev)
Assembly=$(echo $InFile | rev | cut -d "/" -f2 | rev)

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

if [ $2 ]; then
  OutDir=$CurPath/$2
else
  OutDir=$CurPath/repeat_masked/$Organism/$Strain/"$Assembly"_repmask
fi

mkdir -p $WorkDir
cd $WorkDir

RepeatModeler -pa 10 -ninja_dir /scratch/software/NINJA-0.97-cluster_only/NINJA -LTRStruct -database “$Strain”_RepMod

mkdir -p $OutDir
rm -r $WorkDir/RM*
cp -r $WorkDir/* $OutDir/.

rm -r $WorkDir
