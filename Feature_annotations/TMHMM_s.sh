#TMHMM step
# Identifies transmembrane proteins
# Added strain name

for Strain in DSA14_003 DSA15_041; do
	for Proteome in $(ls gene_pred/codingquary/F.oxysporum_fsp_fragariae/$Strain/final/final_genes_combined.pep.fasta); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
		sbatch $ProgDir/TMHMM.sh $Proteome
	done
done
