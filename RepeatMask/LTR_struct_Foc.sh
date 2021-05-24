#!/usr/bin/env bash
#SBATCH -J repeatmasker
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=24

#ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/RepeatMask
#BestAssembly=repeat_masked/F.oxysporum_fsp_fragariae/DSA15_041/ncbi_edits_repmask/DSA15_041_contigs_unmasked.fa
#OutDir=repeat_masked/F.oxysporum_fsp_fragariae/DSA15_041/ncbi_edits_repmask
#sbatch $ProgDir/LTR_struct.sh $BestAssembly $OutDir

InFile=$1
Organism=F.oxysporum_fsp_cepae
Strain=Fus2
Assembly=ncbi_edits

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

if [ $2 ]; then
  OutDir=$CurPath/$2
else
  OutDir=$CurPath/repeat_masked/$Organism/$Strain/"$Assembly"_repmask
fi

mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$InFile "$Strain"_contigs_unmasked.fa
BuildDatabase -name "$Strain"_RepMod "$Strain"_contigs_unmasked.fa
RepeatModeler -pa 10 -ninja_dir /scratch/software/NINJA-0.97-cluster_only/NINJA -LTRStruct -database "$Strain"_RepMod

mkdir -p $OutDir
rm -r $WorkDir/RM*
cp -r $WorkDir/* $OutDir/.

rm -r $WorkDir
