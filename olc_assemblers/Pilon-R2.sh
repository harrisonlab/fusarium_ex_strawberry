for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/Pilon_SDen/race_1_smartdenovo_racon_round_10_renamed.fasta); do
Reads=$(pilon/bwa_mapping.SDEN.sorted.bam)
OutDir= assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/Pilon_SDen/
Iterations=10
ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
sbatch $ProgDir/pilon_longread.sh $Assembly $BAM $OutDir $Iterations
done
