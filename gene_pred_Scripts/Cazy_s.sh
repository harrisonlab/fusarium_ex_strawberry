# CAZY analysis
# Add strain names where Ag02 and co are
# removed line 3 - for Strain in RS305p RS324p; do
# sub_hmmscan.sh needs to be edited and updated

for Strain in DSA14_003 DSA15_041; do
for Proteome in $(ls gene_pred/codingquary/F.oxysporum_fsp_fragariae/$Strain/final/final_genes_combined.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/CAZY/$Organism/$Strain
mkdir -p $OutDir
Prefix="$Strain"_CAZY
CazyHmm=dbCAN/dbCAN-fam-HMMs.txt
ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/
sbatch $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
done
done
