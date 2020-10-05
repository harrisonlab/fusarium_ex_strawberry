# Rename contigs for genome
# If split or remove contigs is needed, provide FCSreport file by NCBI.

ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    touch tmp.txt
    for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/racon_10/race_1_smartdenovo_racon_round_3.fasta); do
        OutDir=$(dirname $Assembly)
        $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/race_1_smartdenovo_racon_round_3_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
    done
    rm tmp.txt
