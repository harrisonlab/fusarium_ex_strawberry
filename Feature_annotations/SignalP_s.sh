#Signal P script for fungi
# Add your strains name to first line
#Added codingquary to Proteome direc

for Strain in DSA14_003 DSA15_041; do
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  CurPath=$PWD
  for Proteome in $(ls gene_pred/codingquary/*/$Strain/final/final_genes_combined.pep.fasta); do
  Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  SplitDir=gene_pred/final_genes_split/$Organism/$Strain
  mkdir -p $SplitDir
  BaseName="$Organism""_$Strain"_final_preds
  $ProgDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName # Split your input fasta in 500 genes files
    for File in $(ls $SplitDir/*_final_preds_*); do
    #sbatch $ProgDir/pred_signalP.sh $File signalp
    #sbatch $ProgDir/pred_signalP.sh $File signalp-3.0 # Recommended for oomycetes
    sbatch $ProgDir/pred_signalP.sh $File signalp-4.1 # Recommended for fungi
    #sbatch $ProgDir/pred_signalP.sh $File signalp-5.0
    done
  done
done
