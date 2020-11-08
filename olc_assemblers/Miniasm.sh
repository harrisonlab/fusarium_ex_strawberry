# For ONT long read sequences use miniasm to assemble genome
# Run in a screen and a node in the
#Run in olc_assemblers conda shell

# If script doesn't work, see below
for TrimReads in $(ls FAL_trim.fastq.gz); do
    Organism=F.oxysporum_fsp_lactucae
    Strain=race_1
    Prefix="$Strain"_miniasm
    OutDir=assembly/miniasm/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
    sbatch $ProgDir/miniasm.sh $TrimReads $Prefix $OutDir
  done

# Step 1
# Login to node srun --partition long --time 0-06:00:00 --mem 40G --cpus-per-task 24 --pty bash
# Concatenate sequence reads first
#Use command below if you are working in the same directory as the raw sequence reads

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

home/gomeza/prog/minimap2/minimap2 -x ava-ont -t8 trimmed_renamed.fasta trimmed_renamed.fasta | gzip -1 > FolR1_fastq_allfiles.paf.gz

# Step 5
# Concatenate pieces of read sequences to generate the final sequences
# Can run like this instead: miniasm -f trimmed_renamed.fasta FolR1_fastq_allfiles.paf.gz > reads.gfa
miniasm -f "$Prefix"_rename.fasta $Prefix.paf.gz > reads.gfa

# Step 6
# Convert gfa file to fasta file
awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > $Prefix.fa
