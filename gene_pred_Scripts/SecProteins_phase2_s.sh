# TMM protein removal

for File in $(ls gene_pred/trans_mem/F.oxysporum_fsp_fragariae/DSA15_041/DSA15_041_TM_genes_neg.txt); do
 Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
 Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
 echo "$Organism - $Strain"
 TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
 cat $File | cut -f1 > $TmHeaders
 SigP=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
 OutDir=$(dirname $SigP)
 ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
 $ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
 cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
done
