#!/usr/bin/env bash
#SBATCH -J swissprot
#SBATCH --partition=long
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=12

 # This script performs blast searches against a database of swissprot genes
 # to determine homology between predicted gene models and previosuly
 # charcaterised gene models.

CurDir=$PWD
WorkDir=$CurDir/${SLURM_JOB_USER}_${SLURM_JOBID}

Proteome=$1
OutDir=$2
SwissDB_Dir=$3
SwissDB_Name=$4

Usage="sub_swissprot.sh <predicted_proteins.fa> <Out_directory> <Swissprot_database_directory> <Swissprot_database_name>"
echo "$Usage"

echo "Query Proteins: $Proteome"
echo "Output directory: $OutDir"
echo "Using Swissprot database from directory: $SwissDB_Dir"
echo "The basename for this database is: $SwissDB_Name"

mkdir -p $WorkDir
cd $WorkDir

cp $CurDir/$Proteome proteins.fa
mkdir -p $WorkDir/db
cd $WorkDir/db
makeblastdb -in ../../../../dbUniprot/swissprot_2020_June/uniprot_sprot.fasta -input_type fasta -dbtype prot -title uniprot_sprot.db -parse_seqids -out uniprot_sprot.db

cd ..

blastp \
  -db db/$SwissDB_Name \
  -query proteins.fa \
  -out swissprot_vJun2020_10_hits.tbl \
  -evalue 1e-100 \
  -outfmt 6 \
  -num_threads 8 \
  -num_alignments 10

ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
$ProgDir/swissprot_parser.py --blast_tbl swissprot_vJun2020_10_hits.tbl --blast_db_fasta db/"$SwissDB_Name".fasta > swissprot_vJun2020_tophit_parsed.tbl

mkdir -p $CurDir/$OutDir
cp -r $WorkDir/*hits.tbl $CurDir/$OutDir/.
cp -r $WorkDir/*parsed.tbl $CurDir/$OutDir/.


