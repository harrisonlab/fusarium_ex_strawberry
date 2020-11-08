# Use for long read assembly programs

# Method was used for Fusarium oxysporum fsp lactucae

####################
# Miniasm assembly
####################
# Step 1
# Open a screen - scrren -r
# Login to node srun --partition long --time 0-18:00:00 --mem-per-cpu 20G --cpus-per-task 24 --pty bash
# Concatenate sequence reads first
# Use command below if you are working in the same directory as the raw sequence reads

# In /archives/2020_niabemr_nanopore/Fo_fragariae_14_003/Fof14-003-Scha48-1/20201027_1719_X1_FAL69458_56dbaa42/fastq_pass$
# /archives/2020_niabemr_nanopore/Fo_fragariae_14_003/Fo_fragariae_14_003_R2/20201028_1622_X1_FAL69458_e6fdc4e2/fastq_pass do
# Did 2 runs therefore  must join both
cat *fastq | gzip -cf > Fof14R1.fastq.gz
cat *fastq | gzip -cf > Fof14R2.fastq.gz

# Copy concatenated file to raw_dna directory
cp Fof14R1.fastq.gz /projects/fusarium_ex_strawberry/NextGenSeq/raw_dna
cp Fof14R2.fastq.gz /projects/fusarium_ex_strawberry/NextGenSeq/raw_dna

# Merge the contatenated files in directory they were copied into
cat Fof14R1.fastq.gz Fof14R2.fastq.gz > Fof14RT.fastq.gz

# Step 2
#Run Porechop before assembly
/scratch/software/Porechop-0.2.3/porechop-runner.py -i Fof14RT.fastq.gz -o Fof14RT_trim.fastq.gz --threads 16 > Fof14RT_trim_log.txt


## For ONT long read sequences use miniasm to assemble genome
# Step 3 Run in olc_assemblers conda shell
# If script doesn't work, see below
for TrimReads in $(ls raw_dna/Fof14RT.fastq.gz); do
    Organism=F.oxysporum_fsp_fragariae
    Strain=DSA14_003
    Prefix="$Strain"_miniasm
    OutDir=assembly/miniasm/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
    sbatch $ProgDir/miniasm.sh $TrimReads $Prefix $OutDir
  done

# Step 3a
#Need to rename all reads
# run in conda env that has minimap2 and bbmap
rename.sh qin=33 in=Fof14RT_trim.fastq.gz out=Fof14RT_renamed.fasta prefix=Fof14

# Step 4a
#Run minimap2
# Run in olc_assemblers conda shell
#Need to run read against itself. IT IS NEEDED. it is going to do self mapping
#Fast all-against-all overlap of raw reads

minimap2 -x ava-ont -t8 Fof14RT_renamed.fasta Fof14RT_renamed.fasta | gzip -1 > Fof14_fastq_allfiles.paf.gz


# Run in a screen and a node in the
# Step 5
# Concatenate pieces of read sequences to generate the final sequences
# Can run like this instead: miniasm -f trimmed_renamed.fasta FolR1_fastq_allfiles.paf.gz > reads.gfa
miniasm -f raw_dna/Fof14RT_renamed.fasta raw_dna/Fof14_fastq_allfiles.paf.gz > assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/reads.gfa

# Step 6
# Convert gfa file to fasta file
awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > $Prefix.fa

####################
# Flye assembly
####################

for TrimReads in $(ls raw_dna/Fof14RT.fastq.gz); do
      Organism=F.oxysporum_fsp_fragariae
      Strain=DSA14_003
       Prefix="$Strain"_flye;     
       TypeSeq=nanoraw;
       OutDir=assembly/flye/$Organism/$Strain/;
       mkdir -p $OutDir;
       Size=60m; # size= Expected genome size
       ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/;
       sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size $TypeSeq;
     done
