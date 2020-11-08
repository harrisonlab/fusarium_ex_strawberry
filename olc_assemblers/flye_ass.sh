# Run in olc_assemblers env
#Log into node
# flye assembly method
# size= Expected genome size
# DO in node and screen
# Use porechop trimmed  output
# Need to include qin=33 or qin=64

for TrimReads in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz); do
       Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev) ;
       Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev) ;
       Prefix="$Strain"_flye;
       TypeSeq=nanoraw;
       OutDir=assembly/flye/$Organism/$Strain/flye_raw;
       mkdir -p $OutDir;
       Size=60m;
       ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/;
       sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size $TypeSeq;
     done


for TrimReads in $(ls FAL_trim.fastq.gz) ; do
    Organism=F.oxysporum_fsp_lactucae
    Strain=race_1
    Prefix="$Strain"_flye
    OutDir=assembly/flye/$Organism/$Strain
    mkdir -p $OutDir
    Size=60m
    ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/olc_assemblers
    sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size
  done

  #original flye directory - /home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
  # edited flye directory /home/akinya/git_repos/fusarium_ex_strawberry/olc_assemblers
  # qin=33 is a flag not an variable like "Organism" need to include differently
  # attempted adding qin before TrimReads - caused error: made it read input
  # added qin after Size - caused error as if it wasn't there
  # added after OutDir - caused error: Estimated genome size - qin=33
  # added to flye.sh script line 48 - caused error
  # error was due to redirection to Antonio's data
  # split script
  # rename.sh qin=33 in=raw_dna/FAL_trim.fastq.gz out=flye1_trimmed_renamed.fasta prefix=FolR1
  # then ran
  # flye --nano-raw assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz --out-dir assembly/flye/F.oxysporum_fsp_lactucae/race_1/ --genome-size 60m --threads 8


# step 1 - screen -a
# step 2 - srun --partition long --time 0-06:00:00 --mem 10G --cpus-per-task 24 --pty bash (parameters force killed job increase memory)
# step 2.a - srun --partition long --time 0-06:00:00 --mem 20G --cpus-per-task 24 --pty bash (job killed increase memory again)
# step 2.b - srun --partition long --time 0-06:00:00 --mem 40G --cpus-per-task 24 --pty bash
# step 3 - change directory
  flye --nano-raw FAL69458.fastq.gz --out-dir assembly/flye/F.oxysporum_fsp_lactucae/race_1 --genome-size 60m --threads 8
