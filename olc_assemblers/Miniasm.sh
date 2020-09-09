#For ONT long read sequences use miniasm
# Run in a screen and a node in the
#Run in olc_assemblers conda shell

#Concatenate sequence reads first
#Use command below if you are working in the same directory as the raw sequence reads

cat *fastq | gzip -cf > FAL69458.fastq.gz

#Run Porechop before assembly
/scratch/software/Porechop-0.2.3/porechop-runner.py -i FAL69458.fastq.gz -o FAL_trim.fastq.gz --threads 16 > FAL_trim_log.txt

#Need to rename all reads
rename.sh qin=33 in=FAL_trim.fastq.gz out=trimmed_renamed.fasta prefix=FolR1

#Run minimap2
#Need to run read against itself. IT IS NEEDED. it is going to do self mapping
#Fast all-against-all overlap of raw reads

home/gomeza/prog/minimap2/minimap2 -x ava-ont -t8 trimmed_renamed.fasta trimmed_renamed.fasta | gzip -1 > FolR1_fastq_allfiles.paf.gz

# Use Miniasm to assemble genome
# Login to node srun --partition long --time 0-06:00:00 --mem 40G --cpus-per-task 24 --pty bash
# Can run like this instead: miniasm -f trimmed_renamed.fasta FolR1_fastq_allfiles.paf.gz > reads.gfa

for TrimReads in $(ls FAL_trim.fastq.gz); do
    Organism=F.oxysporum_fsp_lactucae
    Strain=race_1
    Prefix="$Strain"_miniasm
    OutDir=assembly/miniasm/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/miniasm.sh $TrimReads $Prefix $OutDir
  done
