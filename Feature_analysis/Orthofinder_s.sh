# Orthofinder
# RUN in conda env (Repenv)
# Install paths into profile - then update . .profile
# PATH=${PATH}:/scratch/software/orthomclSoftware-v2.0.9/bin
# PATH=${PATH}:/scratch/software/diamond-0.9.27
# PATH=${PATH}:/scratch/software/fastme-2.1.5/bin/bin
# PATH=${PATH}:/scratch/software/orthomclSoftware-v2.0.9/bin

IsolateAbrv=Fof_14-003
WorkDir=analysis/orthology/$IsolateAbrv
mkdir -p $WorkDir
mkdir -p $WorkDir/formatted

Taxon_code=STR1 # Simple and unique taxon identifier
Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/final/final_genes_appended_renamed.pep.fasta)
Id_field=1 # Indicate what field contains the protein ID
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=ST2 # Simple and unique taxon identifier
Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/final/final_genes_appended_renamed.pep.fasta)
Id_field=1 # Indicate what field contains the protein ID
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
