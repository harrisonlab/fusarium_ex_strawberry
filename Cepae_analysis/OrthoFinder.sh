# Orthofinder
# RUN in conda env (Repenv)
# Install paths into profile - then update . .profile
# PATH=${PATH}:/scratch/software/orthomclSoftware-v2.0.9/bin
# PATH=${PATH}:/scratch/software/diamond-0.9.27
# PATH=${PATH}:/scratch/software/fastme-2.1.5/bin/bin
# PATH=${PATH}:/scratch/software/orthomclSoftware-v2.0.9/bin

IsolateAbrv=Foc_F2CN
WorkDir=analysis/orthology/$IsolateAbrv
mkdir -p $WorkDir
mkdir -p $WorkDir/formatted

Taxon_code=STR1 # Simple and unique taxon identifier
Fasta_file=$(ls final/final_genes_combined.pep.fasta)
Id_field=1 # Indicate what field contains the protein ID
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=ST2 # Simple and unique taxon identifier
Fasta_file=$(ls final/final_genes_combined.pep.fasta)
Id_field=1 # Indicate what field contains the protein ID
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

# Main phase Orthofinder
# Run in conda env (Repenv)

for IN_DIR in $(ls -d $WorkDir/formatted) ; do
sbatch /home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis/orthofinder.sh $IN_DIR $IsolateAbrv
done
