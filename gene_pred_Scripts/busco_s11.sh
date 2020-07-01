#Edited paths from Andy's direc to mine
#Updated entire script using https://github.com/harrisonlab/bioinformatics_tools/blob/master/Gene_prediction/README.md
#Does not have quast segment in script like Andy's pipeline
#Look into BuscoDB direc - directory exists
#Run in conda env - BUSCOenv
#Ran on genome(softmasked) and gene models (final_genes_appended_renamed.gene.fasta)

for Assembly in $(ls gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/final/final_genes_appended_renamed.gene.fasta); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
done
