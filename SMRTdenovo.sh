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
