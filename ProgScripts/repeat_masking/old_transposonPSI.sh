#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

Usage="qsub transposonPSI.sh <assembled_contigs.fa> [<output_directory>]"

InFile=$1

Organism=$(echo $InFile | rev | cut -d "/" -f4 | rev)
Strain=$(echo $InFile | rev | cut -d "/" -f3 | rev)
Assembly=$(echo $InFile | rev | cut -d "/" -f2 | rev)

CurPath=$PWD
WorkDir=$TMPDIR/"$Strain"_transposon_psi

if [ $2 ]; then
  OutDir=$CurPath/$2
else
  OutDir=$CurPath/repeat_masked/$Organism/$Strain/"$Assembly"_repmask
fi

mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$InFile "$Strain"_contigs_unmasked.fa

transposonPSI.pl "$Strain"_contigs_unmasked.fa nuc

mkdir -p $OutDir
cp -r $WorkDir/* $OutDir/.

rm -r $TMPDIR
