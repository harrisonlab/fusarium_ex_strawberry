#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 12
#$ -l virtual_free=0.9G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

# This script uses repeatmodeler and repeatmasker to mask Interspersed repeats
# and low complexity regions within the genome. Firstly, repeatmodeler identifies
# repeat element boundaries and relationships within repeat families. The repeats
# identified within the genome are provided to repeatmasker, which uses this data
# along with it's own repeat libraries to identify these repetitive regions and
# perform masking. Masking is done at 3 levels:
# Hardmasking = repetitive sequence is replaced with N's.
# Softmasking = repetitive sequence is converted to lower case.
# Ignoring low-complexity regions = only interspersed repetitive elements are masked.

Usage="sbatch rep_modeling.sh <assembled_contigs.fa> [<output_directory>]"

InFile=$1
Orgasnism=$(echo $InFile | rev | cut -d "/" -f4 | rev)
Strain=$(echo $InFile | rev | cut -d "/" -f3 | rev)
Assembly=$(echo $InFile | rev | cut -d "/" -f2 | rev)

CurPath=$PWD
WorkDir=$TMPDIR/"$Strain"_repeatmask

if [ $2 ]; then
  OutDir=$CurPath/$2
else
  OutDir=$CurPath/repeat_masked/$Orgasnism/$Strain/"$Assembly"_repmask
fi

mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$InFile "$Strain"_contigs_unmasked.fa
BuildDatabase -name "$Strain"_RepMod "$Strain"_contigs_unmasked.fa
RepeatModeler -pa 16 -database "$Strain"_RepMod

# hardmask
RepeatMasker -gff -pa 16 -lib RM_*.*/consensi.fa.classified "$Strain"_contigs_unmasked.fa
mv "$Strain"_contigs_unmasked.fa.cat.gz "$Strain"_contigs_hardmasked.fa.cat.gz
mv "$Strain"_contigs_unmasked.fa.masked "$Strain"_contigs_hardmasked.fa
mv "$Strain"_contigs_unmasked.fa.out "$Strain"_contigs_hardmasked.out
mv "$Strain"_contigs_unmasked.fa.out.gff "$Strain"_contigs_hardmasked.fa.out.gff
mv "$Strain"_contigs_unmasked.fa.tbl "$Strain"_contigs_hardmasked.tbl
grep -v '#' "$Strain"_contigs_hardmasked.fa.out.gff > "$Strain"_contigs_hardmasked.gff


# softmask
RepeatMasker -xsmall -gff -pa 16 -lib RM_*.*/consensi.fa.classified "$Strain"_contigs_unmasked.fa
mv "$Strain"_contigs_unmasked.fa.cat.gz "$Strain"_contigs_softmasked.fa.cat.gz
mv "$Strain"_contigs_unmasked.fa.masked "$Strain"_contigs_softmasked.fa
mv "$Strain"_contigs_unmasked.fa.out "$Strain"_contigs_softmasked.out
mv "$Strain"_contigs_unmasked.fa.out.gff "$Strain"_contigs_softmasked.fa.out.gff
mv "$Strain"_contigs_unmasked.fa.tbl "$Strain"_contigs_softmasked.tbl
grep -v '#' "$Strain"_contigs_softmasked.fa.out.gff > "$Strain"_contigs_softmasked.gff

# don't mask low-complexity or simple-repeat sequences, just transposons
RepeatMasker -nolow -gff -pa 16 -lib RM_*.*/consensi.fa.classified "$Strain"_contigs_unmasked.fa
mv "$Strain"_contigs_unmasked.fa.cat.gz "$Strain"_contigs_transposonmasked.fa.cat.gz
mv "$Strain"_contigs_unmasked.fa.masked "$Strain"_contigs_transposonmasked.fa
mv "$Strain"_contigs_unmasked.fa.out "$Strain"_contigs_transposonmasked.out
mv "$Strain"_contigs_unmasked.fa.out.gff "$Strain"_contigs_transposonmasked.fa.out.gff
mv "$Strain"_contigs_unmasked.fa.tbl "$Strain"_contigs_transposonmasked.tbl
grep -v '#' "$Strain"_contigs_transposonmasked.fa.out.gff > "$Strain"_contigs_transposonmasked.gff

mkdir -p $OutDir
rm -r $WorkDir/RM*
cp -r $WorkDir/* $OutDir/.

rm -r $TMPDIR
