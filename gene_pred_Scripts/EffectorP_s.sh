# EffectorP - Effector identification
# Add path to .profile PATH=${PATH}:/scratch/software/EffectorP-2.0/Scripts
# EffectorP did not work first time
# Do not run as script!
# EffectorP.py -o "$BaseName".txt -E "$BaseName".fa -i $Proteome using full paths
# develops .txt and .fa files

mkdir DSA15_041
mv DSA15_041 analysis/effectorP/F.oxysporum_fsp_fragariae/
EffectorP.py -o "$BaseName".txt -E "$BaseName".fa -i $Proteome

#for Proteome in $(ls gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/final/final_genes_appended_renamed.pep.fasta); do
#Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
#Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
#BaseName="$Organism"_"$Strain"_EffectorP
#OutDir=analysis/effectorP/$Organism/$Strain
#mv "$BaseName".txt $OutDir
#mv "$BaseName".fa $OutDir
#ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
#sbatch $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
#done
