# Use for long read assembly programs

# Method was used for Fusarium oxysporum fsp lactucae

####################
# Miniasm assembly
####################
# Step 1
# Login to node srun --partition long --time 0-06:00:00 --mem 40G --cpus-per-task 24 --pty bash
# Concatenate sequence reads first
# Use command below if you are working in the same directory as the raw sequence reads

cat *fastq | gzip -cf > FAL69458.fastq.gz

# Step 2
#Run Porechop before assembly
/scratch/software/Porechop-0.2.3/porechop-runner.py -i FAL69458.fastq.gz -o FAL_trim.fastq.gz --threads 16 > FAL_trim_log.txt

# Step 3
#Need to rename all reads
rename.sh qin=33 in=FAL_trim.fastq.gz out=trimmed_renamed.fasta prefix=FolR1

# Step 4
#Run minimap2
#Need to run read against itself. IT IS NEEDED. it is going to do self mapping
#Fast all-against-all overlap of raw reads

minimap2 -x ava-ont -t8 trimmed_renamed.fasta trimmed_renamed.fasta | gzip -1 > FolR1_fastq_allfiles.paf.gz

# For ONT long read sequences use miniasm to assemble genome
# Run in a screen and a node in the
# Run in olc_assemblers conda shell


# If script doesn't work, see below
for TrimReads in $(ls FAL_trim.fastq.gz); do
    Organism=F.oxysporum_fsp_lactucae
    Strain=race_1
    Prefix="$Strain"_miniasm
    OutDir=assembly/miniasm/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/miniasm.sh $TrimReads $Prefix $OutDir
  done

# Step 5
# Concatenate pieces of read sequences to generate the final sequences
# Can run like this instead: miniasm -f trimmed_renamed.fasta FolR1_fastq_allfiles.paf.gz > reads.gfa
miniasm -f "$Prefix"_rename.fasta $Prefix.paf.gz > reads.gfa

# Step 6
# Convert gfa file to fasta file
  awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > $Prefix.fa


#######################
# Flye assembly
#######################

# Run in olc_assemblers env
# Log into node
# flye assembly method
# size= Expected genome size

for TrimReads in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz); do
       Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev) ;
       Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev) ;
       Prefix="$Strain"_flye;     
       TypeSeq=nanoraw;
       OutDir=assembly/flye/$Organism/$Strain/flye_raw;
       mkdir -p $OutDir;
       Size=60m; # size= Expected genome size
       ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/;
       sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size $TypeSeq;
     done

########################
# SMARTDenovo assembly
########################

# SMartDenovo script
# run in olc_assemblers env
#Use raw inputs unlike miniasm
# Use porechop trimmed output
# Could run like this:
# ~/miniconda3/envs/olc_assemblers/bin/smartdenovo.pl -p race_1_smartdenovo -t 14 -c 1 FolR1_fastq_allfiles.paf.gz > race_1_smartdenovo.mak
# make -f prefix.mak

  for TrimReads in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz); do
    Organism=F.oxysporum_fsp_lactucae
    Strain=race_1
    Prefix="$Strain"_smartdenovo
    OutDir=assembly/SMARTdenovo/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/SMARTdenovo.sh $TrimReads $Prefix $OutDir
  done

# output = race_1_smartdenovo.dmo.lay.utg

#####################
# QC steps
#####################


# quast QC assembly check
# Run in conda env with python 2.7 (betaenv)
# Run on each assembly

ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/assembly.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=assembly/miniasm/$Organism/$Strain/ncbi_edits
    sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
  done
#Updated entire script using https://github.com/harrisonlab/bioinformatics_tools/blob/master/Gene_prediction/README.md
#Look into BuscoDB direc - directory exists
#Run in conda env - BUSCOenv
#Ran on genome(softmasked) and gene models (final_genes_appended_renamed.gene.fasta)

for Assembly in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/assembly.fasta); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
  echo "$Organism - $Strain"
  ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
  BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
  OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
  sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
done
