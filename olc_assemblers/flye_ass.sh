# Run in olc_assemblers
# flye assembly method
# size= Expected genome size
# DO in node and screen

for TrimReads in $(ls raw_dna/FAL69458.fastq.gz) ; do
    Organism=F.oxysporum_fsp_lactucae
    Strain=race_1
    Prefix="$Strain"_flye
    OutDir=assembly/flye/$Organism/$Strain
    mkdir -p $OutDir
    Size=60m
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size
  done

# step 1 - screen -a
# step 2 - srun --partition long --time 0-06:00:00 --mem 10G --cpus-per-task 24 --pty bash (parameters force killed job increase memory)
# step 2.a - srun --partition long --time 0-06:00:00 --mem 20G --cpus-per-task 24 --pty bash (job killed increase memory again)
# step 2.b - srun --partition long --time 0-06:00:00 --mem 40G --cpus-per-task 24 --pty bash
# step 3 - change directory
  flye --nano-raw FAL69458.fastq.gz --out-dir assembly/flye/F.oxysporum_fsp_lactucae/race_1 --genome-size 60m --threads 8
