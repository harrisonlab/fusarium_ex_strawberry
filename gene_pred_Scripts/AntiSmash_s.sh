#Antismash
# Run line by line
# RUn like AddTrans.sh

for AntiSmash in $(ls -d gene_pred/antiSMASH/F.oxysporum/DSA14_003/DSA14_003_contigs_unmasked.gbk); do
  Organism=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $AntiSmash | rev | cut -f2 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
  Prefix=$OutDir/"Strain"_antismash_results
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  # Antismash v5 output to gff file
  $ProgDir/antismash2gffv5.py --inp_antismash $AntiSmash --out_prefix $Prefix
  #$ProgDir/antismash2gff.py --inp_antismash $AntiSmash --out_prefix $Prefix # Use only for antismash v4.2 output
  printf "Number of secondary metabolite detected:\t"
  cat "$Prefix"_secmet_clusters.gff | wc -l
  GeneGff=gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/final/final_genes_appended_renamed.gff3
  bedtools intersect -u -a $GeneGff -b "$Prefix"_secmet_clusters.gff > "$Prefix"_secmet_genes.gff
  cat "$Prefix"_secmet_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_antismash_secmet_genes.txt
  bedtools intersect -wo -a $GeneGff -b "$Prefix"_secmet_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -e "s/;Parent=g\w+//g" | perl -p -e "s/;Notes=.*//g" > "$Prefix"_secmet_genes.tsv
  printf "Number of predicted proteins in secondary metabolite clusters:\t"
  cat analysis/secondary_metabolites/antismash/F.oxysporum/DSA14_003/Strain_antismash_results_antismash_secmet_genes.txt | wc -l
  printf "Number of predicted genes in secondary metabolite clusters:\t"
  cat "$Prefix"_secmet_genes.gff | grep -w 'gene' | wc -l
done

# Antismash output correction. Some gene names contain ;. Remove manually with the following command.
# First sed command removes ;. Second and Third remove the cluster kind information (optional)
cat analysis/secondary_metabolites/antismash/F.oxysporum/DSA14_003/Strain_antismash_results_secmet_genes.tsv | sed 's/;//p' | sed 's/;.*//p' | sed 's/Kin.*//p' > analysis/secondary_metabolites/antismash/F.oxysporum/DSA14_003/"Strain"_antismash_results_secmet_genes_corrected.tsv
#Edit output file names from this script after completion
