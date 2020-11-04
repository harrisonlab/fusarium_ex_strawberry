#Edited paths from Andy's direc to mine
#Updated entire script using https://github.com/harrisonlab/bioinformatics_tools/blob/master/Gene_prediction/README.md
#Does not have quast segment in script like Andy's pipeline
#Look into BuscoDB direc - directory exists
#Run in conda env - BUSCOenv
#Ran on genome(softmasked) and gene models (final_genes_appended_renamed.gene.fasta)

for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/Pilon_SDen/FolR1_SDen_pilon.fasta); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
ProgDir=~/git_repos/fusarium_ex_strawberry/ProgScripts
BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
done
