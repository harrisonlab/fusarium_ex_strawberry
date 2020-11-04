# Aligning illumina reads against pilon data
# Run in conda env
# Make sure script is executable
# chmod u+x ./sub_pilon.sh

for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/Pilon_SDen/race_1_smartdenovo_racon_round_10_renamed.fasta); do
  FileF=$(ls ../oldhome/groups/harrisonlab/project_files/fusarium/raw_dna/paired/F.oxysporum_fsp_lactucae/AJ520/F/AJ520_S2_L001_R1_001.fastq.gz)
  FileR=$(ls ../oldhome/groups/harrisonlab/project_files/fusarium/raw_dna/paired/F.oxysporum_fsp_lactucae/AJ520/R/AJ520_S2_L001_R2_001.fastq.gz)
OutDir=assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/Pilon_SDen/
Iterations=10
ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
sbatch $ProgDir/sub_pilon.sh $Assembly $FileF $FileR $OutDir $Iterations
done
