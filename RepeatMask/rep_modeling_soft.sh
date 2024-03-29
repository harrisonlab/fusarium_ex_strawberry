#!/usr/bin/env bash
#SBATCH -J repeatmasker
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=24G
#SBATCH --cpus-per-task=24

# This script uses repeatmodeler and repeatmasker to mask Interspersed repeats
# and low complexity regions within the genome. Firstly, repeatmodeler identifies
# repeat element boundaries and relationships within repeat families. The repeats
# identified within the genome are provided to repeatmasker, which uses this data
# along with it's own repeat libraries to identify these repetitive regions and
# perform masking. Masking is done at 3 levels:
# Hardmasking = repetitive sequence is replaced with N's.
# Softmasking = repetitive sequence is converted to lower case.
# Ignoring low-complexity regions = only interspersed repetitive elements are masked.

# Edition notes
# Note 1 - The latest version of RepeatModeler (RepeatModeler 2.0.1) copies the final classified results to the families.stk and consensus.fa.
# RepeatClassifier is needed in this version to produced consensi.fa.classified files.
# Note 2 - No nucleotide repeat library is included in RepeatMasker. It is recommended to download a database from RebBase but a licence is needed.
# The RepBase v19.09 was copied to /home/gomeza/RepMaskerDB

Usage="sbatch rep_modeling.sh <assembled_contigs.fa> [<output_directory>]"

InFile=$1
Organism=$2
Strain=$3
Assembly=ncbi_edits

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

if [ $4 ]; then
  OutDir=$CurPath/$4
else
  OutDir=$CurPath/Repeat_masked_2/$Organism/$Strain/miniasm/"$Assembly"_repmask
fi

mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$InFile "$Strain"_contigs_unmasked.fa
BuildDatabase -name "$Strain"_RepMod "$Strain"_contigs_unmasked.fa
RepeatModeler -pa 16 -ninja_dir /scratch/software/NINJA-0.97-cluster_only/NINJA -LTRStruct -database "$Strain"_RepMod


RepeatClassifier -consensi RM_*.*/consensi.fa -stockholm RM_*.*/families.stk

# softmask
RepeatMasker -xsmall -gff -pa 16 -lib RM_*.*/consensi.fa.classified "$Strain"_contigs_unmasked.fa
mv "$Strain"_contigs_unmasked.fa.cat.gz "$Strain"_contigs_softmasked.fa.cat.gz
mv "$Strain"_contigs_unmasked.fa.masked "$Strain"_contigs_softmasked.fa
mv "$Strain"_contigs_unmasked.fa.out "$Strain"_contigs_softmasked.out
mv "$Strain"_contigs_unmasked.fa.out.gff "$Strain"_contigs_softmasked.fa.out.gff
mv "$Strain"_contigs_unmasked.fa.tbl "$Strain"_contigs_softmasked.tbl
grep -v '#' "$Strain"_contigs_softmasked.fa.out.gff > "$Strain"_contigs_softmasked.gff


mkdir -p $OutDir
rm -r $WorkDir/RM*
cp -r $WorkDir/* $OutDir/.

rm -r $WorkDir
