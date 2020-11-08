# fusarium_ex_strawberry
Commands used in the analysis of Fusarium oxysporum isolates ex. strawberry.
test text

# Assembling DNA utilising Pipelines developed by Andy Armitage - slightly modded to adapt to NewSlurm from old GRDENG

# Data qc
# Data quality was visualised using fastqc:

  for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc; # directory has been changed
    sbatch $ProgDir/run_fastqc.sh $RawData;
  done


# Trimming was performed on data to trim adapters from sequences and remove poor quality data.

for StrainPath in $(ls -d raw_dna/paired/*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ReadsF=$(ls $StrainPath/F/*.fastq*)
ReadsR=$(ls $StrainPath/R/*.fastq*)
echo $ReadsF
echo $ReadsR
sbatch $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
done

# Data quality was visualised once again following trimming:

  for RawData in $(ls qc_dna/paired/*/*/*/*.fq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    sbatch $ProgDir/run_fastqc.sh $RawData;
  done

# Sequencing coverage was estimated:

for RawData in $(ls qc_dna/paired/*/*/*/*fq.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData
GenomeSz=35
OutDir=$(dirname $RawData)
sbatch $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done

# Find predicted coverage for these isolates:


for StrainDir in $(ls -d qc_dna/paired/*/*); do
  Strain=$(basename $StrainDir)
  printf "$Strain\t"
  for File in $(ls qc_dna/paired/*/"$Strain"/*/*.txt); do
  echo $(basename $File);
  cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
  done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
done

# Spades assembly step
# Assembles illumina sequenced DNA
for StrainPath in $(ls -d qc_dna/paired/*/*)
do
    ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/spades
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    OutDir=assembly/spades/$Organism/${Strain}
    echo $F_Read
    echo $R_Read
    sbatch $ProgDir/slurm_spades.sh $F_Read $R_Read $OutDir correct 10
done

# Quast assembly QC step
ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast

for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta)
do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  # OutDir=assembly/spades/$Organism/$Strain
  OutDir=$(dirname $Assembly)

  sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
done

# Develop result files  from quast analysis
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/transposed_report.tsv); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    # echo;
    # echo $Organism;
    # echo $Strain;
    cat $Assembly | tail -n +2 | sed "s/contigs_min_500bp/${Organism}_${Strain}/g"
  done > assembly/quast_results.txt

# Remove contaminants from assembled illumina reads and renames them
# Run in a conda env that works with python
ProgDir=$(ls -d /home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants)
  touch tmp.csv
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
    mkdir $OutDir
    $ProgDir/remove_contaminants2.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv

# Generates files for NCBI submission
# Not essential but helps to keep0 a standard for contig names
for Assembly in $(ls assembly/spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -v 'ncbi_edits'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    NCBI_report_dir=genome_submission/$Organism/$Strain/initial_submission
    mkdir -p $NCBI_report_dir
done

# Step to "fix" contigs by removi9ng contaminats from the renamed fasta files
for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  echo "$Organism - $Strain"
  NCBI_report=$(ls genome_submission/$Organism/$Strain/initial_sub/Contamination*.txt)
  OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
  mkdir -p $OutDir
  ProgDir=~/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  $ProgDir/remove_contaminants2.py --keep_mitochondria --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
done

# quast contig fix step
# Qc  check

ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamedA.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
    sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
  done
