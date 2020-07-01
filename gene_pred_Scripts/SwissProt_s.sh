#Copy to progscript dir
#edited SwissDbDir path as it caused an error
#BLAST Database error: No alias or index file found for protein database [db/uniprot_sprot] in search path [/tmp/swissprot::]
#SwissDb has yet to be installed - looking through manual
#/home/akinya/F.oxysporum_Ref_proteomes
#../oldhome/groups/harrisonlab/uniprot/swissprot_2018_March
#../dbUniprot/swissprot_2020_June
# /home/akinya/git_repos/assembly_fusarium_ex/ProgScripts
# /home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation

for Proteome in $(ls gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot20/$Organism/$Strain
SwissDbDir=../dbUniprot/swissprot_2020_June
SwissDbName=uniprot_sprot
ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/ProgScripts
sbatch $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done

#SwissPort part 2
#when finished testing scrips split into 2 separate executables
for SwissTable in $(ls gene_pred/swissprot/*/*/*_hits.tbl); do
		Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_v2017_tophit_parsed.tbl
		ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		$ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta ../../uniprot/swissprot/uniprot_sprot.fasta > $OutTable
	done
