# Racon
# Need to do for Miniasm*, Flye* and SMARTdenovo* output files
# *=complete

for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/race_1_smartdenovo.dmo.lay.utg); do
    ReadsFq=$(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz)
    Iterations=10
    OutDir=$(dirname $Assembly)"/racon_$Iterations"
    ProgDir=~/git_repos/assembly_fusarium_ex/ProgScripts
    sbatch $ProgDir/racon.sh $Assembly $ReadsFq $Iterations $OutDir
  done
