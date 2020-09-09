# SMartDenovo script
#Use raw inputs unlike miniasm
# Could run like this:
# ~/miniconda3/envs/olc_assemblers/bin/smartdenovo.pl -p race_1_smartdenovo -t 14 -c 1 FolR1_fastq_allfiles.paf.gz > race_1_smartdenovo.mak
make -f prefix.mak

  for TrimReads in $(ls raw_dna/FAL69458.fastq.gz); do
    Organism=F.oxysporum_fsp_lactucae
    Strain=race_1
    Prefix="$Strain"_smartdenovo
    OutDir=assembly/SMARTdenovo/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
    sbatch $ProgDir/SMARTdenovo.sh $TrimReads $Prefix $OutDir
  done
