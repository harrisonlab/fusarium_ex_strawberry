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
Organism=F.oxysporum_fsp_cepae
Strain=Fus2_canu_new
contig=contig_21

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

if [ $2 ]; then
  OutDir=$CurPath/$2
else
  OutDir=$CurPath/repeat_masked/$contig
fi

mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$InFile "$contig"_unmasked.fa
BuildDatabase -name "$contig"_RepMod "$contig"_unmasked.fa
RepeatModeler -pa 16 -database "$contig"_RepMod

RepeatClassifier -consensi RM_*.*/consensi.fa -stockholm RM_*.*/families.stk

#transposons
RepeatMasker -nolow -gff -pa 16 -lib RM_*.*/consensi.fa.classified "$contig"_unmasked.fa
mv "$contig"_unmasked.fa.cat.gz "$contig"_transposonmasked.fa.cat.gz
mv "$contig"_unmasked.fa.masked "$contig"_transposonmasked.fa
mv "$contig"_unmasked.fa.out "$contig"_transposonmasked.out
mv "$contig"_unmasked.fa.out.gff "$contig"_transposonmasked.fa.out.gff
mv "$contig"_unmasked.fa.tbl "$contig"_transposonmasked.tbl
grep -v '#' "$contig"_transposonmasked.fa.out.gff > "$contig"_transposonmasked.gff

mkdir -p $OutDir
rm -r $WorkDir/RM*
cp -r $WorkDir/* $OutDir/.

rm -r $WorkDir
