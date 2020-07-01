#Secreted proteins merge
# Added strains name
# Ran line by line 
# Use full path for SplitDir

for Strain in DSA14_003; do
   for SplitDir in $(ls -d gene_pred/final_genes_split/F.oxysporum_fsp_fragariae/$Strain); do
   	Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
   	Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
   	InStringAA=''
   	InStringNeg=''
   	InStringTab=''
   	InStringTxt=''
   	SigpDir=final_genes_signalp-4.1
   	for GRP in $(ls -l $SplitDir/*_final_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do
   		InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.aa";
   		InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp_neg.aa";
   		InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.tab";
   		InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.txt";
   	done
   	cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.aa
   	cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_neg_sp.aa
   	tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.tab
   	cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.txt
   done
done
