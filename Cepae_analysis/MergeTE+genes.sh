# Anne helped concatenate TEs and genes
mkdir akin
cd akin
cp /projects/fusarium_ex_strawberry/NextGenSeq/Fo_Cepae/repeat_masked/Fus2_canu_new_contigs_unmasked.fa.TPSI.allHits.chains.bestPerLocus.gff3 .
cp /projects/fusarium_ex_strawberry/NextGenSeq/Fo_Cepae/gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 .
awk 'BEGIN {OFS="\t"}{print $1, $4, $5, $3, $7, $8, $9, $2}' final_genes_appended.gff3 > final_genes_appended.edited.gff3
awk 'BEGIN {OFS="\t"}{print $1, $4, $5, $3, $7, $8, $9, $2}' Fus2_canu_new_contigs_unmasked.fa.TPSI.allHits.chains.bestPerLocus.gff3 >Fus2_canu_new_contigs_unmasked.fa.TPSI.allHits.chains.bestPerLocus.edited.gff3
bedtools intersect -wo -a Fus2_canu_new_contigs_unmasked.fa.TPSI.allHits.chains.bestPerLocus.edited.gff3 -b final_genes_appended.edited.gff3 >akin.gff3
cut -f 1-7,9-15,17 akin.gff3 >transposon_to_genes_all_features.gff3
mkdir /scratch/projects/annew/filedump
mv transposon_to_genes_all_features.gff3 Fus2_canu_new_contigs_unmasked.fa.TPSI.allHits.chains.bestPerLocus.edited.gff3 final_genes_appended.edited.gff3 akin.gff3 /scratch/projects/annew/filedump
