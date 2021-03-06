#Run line by line
#Run in conda env (Repenv)
#Remove duplicate and rename genes

  GffAppended=$(ls -d gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3)
  Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
  Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  FinalDir=gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final

  # Remove duplicated genes
  GffFiltered=$FinalDir/filtered_duplicates.gff
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
  $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered

# Rename genes
  GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
  LogFile=$FinalDir/final_genes_appended_renamed.log
  $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
  rm $GffFiltered

# Create renamed fasta files from each gene feature
  Assembly=$(ls repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_softmasked_repeatmasker_TPSI_appended.fa)
  $ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_appended_renamed
  # The proteins fasta file contains * instead of Xs for stop codons, these should be changed
  sed -i 's/\*/X/g' $FinalDir/final_genes_appended_renamed.pep.fasta

# view gene names
 cat $FinalDir/final_genes_appended_renamed.cdna.fasta | grep '>'
