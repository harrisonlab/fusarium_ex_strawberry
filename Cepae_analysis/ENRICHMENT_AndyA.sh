5.4 Enrichment of genes on LS contigs
Fus2 contigs 10, 14, 16, 19, 20, 21 and 22 were identified as lineage specific.

As such interproscan annotation were extracted for these contigs and Functional enrichments were compared to the entire genome.

Enrichment of effector families by contig
GffGenes=gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3
cat $GffGenes | grep 'mRNA' | cut -f1 | sort | uniq -c | sort -n -k2 -t '_'
GffEffP=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted.gff
cat $GffEffP | grep 'mRNA' | cut -f1 | sort | uniq -c | sort -n -k2 -t '_'
GffCAZY=gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted.gff
cat $GffCAZY | grep 'mRNA' | cut -f1 | sort | uniq -c | sort -n -k2 -t '_'
GffSecMet=analysis/antismash/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_secondary_metabolite_regions.gff
cat $GffSecMet | cut -f1 | sort | uniq -c | sort -n -k2 -t '_'
GffMimp=analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp.gff
cat $GffMimp | grep 'mRNA' | cut -f1 | sort | uniq -c | sort -n -k2 -t '_'
GffMimpSec=analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp_secreted.gff
cat $GffMimpSec | grep 'mRNA' | cut -f1 | sort | uniq -c | sort -n -k2 -t '_'

##############################

From this it was clear that the genes within 2kb of a mimp were primarily associated with the pathogen specific contigs. INterproscan information on these genes was extracted

MimpGeneHeaders=analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_prots_in_2kb_mimp.txt
MimpInterpro=analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_prots_in_2kb_mimp_interpro.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/feature_extraction
$ProgDir/filter_by_gene.py --inp_interpro $InterPro --inp_txt $MimpGeneHeaders > $MimpInterpro

###########################

TotalMimps=analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_mimps.gff
# Total in Genome
cat $TotalMimps |wc -l
# In all LS regions
cat $TotalMimps | grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' -e 'contig_14' -e 'contig_20' -e 'contig_22' | wc -l
# In PS regions
cat $TotalMimps | grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' | wc -l
GenesIn2Kb=analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp.gff
# Total in Genome
cat $GenesIn2Kb | grep 'gene'| wc -l
# In all LS regions
cat $GenesIn2Kb | grep 'gene'| grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' -e 'contig_14' -e 'contig_20' -e 'contig_22' | wc -l
# In PS regions
cat $GenesIn2Kb | grep 'gene'| grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' | wc -l
SecretedIn2Kb=analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp_secreted.gff
# Total in Genome
cat $SecretedIn2Kb | grep 'gene'| wc -l
# In all LS regions
cat $SecretedIn2Kb | grep 'gene'| grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' -e 'contig_14' -e 'contig_20' -e 'contig_22' | wc -l
# In PS regions
cat $SecretedIn2Kb | grep 'gene'| grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21'  | wc -l

####################################

GenesIn2Kb=analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp.gff
OutDir=analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new
cat $GenesIn2Kb | grep 'gene'| grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' | cut -f2 -d '=' | tr -d ';' > $OutDir/Fus2_canu_new_genes_in_2kb_mimp_PS_headers.txt
cat $GenesIn2Kb | grep 'gene'| grep -e 'contig_14' -e 'contig_20' -e 'contig_22' | cut -f2 -d '=' | tr -d ';' > $OutDir/Fus2_canu_new_genes_in_2kb_mimp_non-PS_LS_headers.txt
cat $GenesIn2Kb | grep 'gene'| grep -v -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' -e 'contig_14' -e 'contig_20' -e 'contig_22' | cut -f2 -d '=' | tr -d ';' > $OutDir/Fus2_canu_new_genes_in_2kb_mimp_core_headers.txt

EffectorPHeaders=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_headers.txt
cat $EffectorPHeaders | cut -f1 | cut -f1 -d '.' | sort | uniq | grep -f $OutDir/Fus2_canu_new_genes_in_2kb_mimp_PS_headers.txt | wc -l
cat $EffectorPHeaders | cut -f1 | cut -f1 -d '.' | sort | uniq | grep -f $OutDir/Fus2_canu_new_genes_in_2kb_mimp_non-PS_LS_headers.txt | wc -l
cat $EffectorPHeaders | cut -f1 | cut -f1 -d '.' | sort | uniq | grep -f $OutDir/Fus2_canu_new_genes_in_2kb_mimp_core_headers.txt | wc -l

CazyHeaders=gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted_headers.txt
cat $CazyHeaders | cut -f1 | cut -f1 -d '.' | sort | uniq | grep -f $OutDir/Fus2_canu_new_genes_in_2kb_mimp_PS_headers.txt | wc -l
cat $CazyHeaders | cut -f1 | cut -f1 -d '.' | sort | uniq | grep -f $OutDir/Fus2_canu_new_genes_in_2kb_mimp_non-PS_LS_headers.txt | wc -l
cat $CazyHeaders | cut -f1 | cut -f1 -d '.' | sort | uniq | grep -f $OutDir/Fus2_canu_new_genes_in_2kb_mimp_core_headers.txt | wc -l

####################################

AntismashGff=analysis/antismash/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_secondary_metabolite_regions.gff
MimpPlus2KbGff=analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_mimps_exp.gff
bedtools intersect -u -a $AntismashGff -b $MimpPlus2KbGff | less -S

Enrichment was performed for each contig:

for num in $(seq 1 34); do
Contig="contig_"$num"_pilon"
echo $Contig
OutDir=analysis/enrichment/LS_regions/F.oxysporum_fsp_cepae/Fus2_canu_new/by_contig/$Contig
mkdir -p $OutDir
Gff=gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3
GeneList=$OutDir/"$Contig"_interproscan.txt
cat $Gff | grep "$Contig" | grep 'mRNA' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '='  > $GeneList
InterPro=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
InterPro_contig=$OutDir/"$Contig"_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/feature_extraction
$ProgDir/filter_by_gene.py --inp_interpro $InterPro --inp_txt $GeneList > $InterPro_contig
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/feature_extraction
$ProgDir/interpro2wego.py --inp_interpro $InterPro_contig > $OutDir/"$Contig"_GO_WEGO.txt
done

#########################

The WEGO wesite was used to generate summaries of genes on core and LS regions.

A similar process was then repeated for duplicated genes:

OutDir=analysis/enrichment/LS_regions/F.oxysporum_fsp_cepae/Fus2_canu_new/duplicated_genes
mkdir -p $OutDir
DupGenesTab=/home/sobczm/popgen/codon/blast/dagchainer/testing/no_transposon/non-transposon_duplications_summaryf_ann
cat $DupGenesTab | tail -n+2 | cut -f1,3,5 | sed 's/\t/\n/g' | sed 's/,/\n/g' | grep -v "^$" | grep 'Fus2' | sed 's/Fus2_//g' | sort | uniq > $OutDir/duplicated_genes.txt
InterPro=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
InterPro_dups=$OutDir/duplicated_genes_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/feature_extraction
$ProgDir/filter_by_gene.py --inp_interpro $InterPro --inp_txt $OutDir/duplicated_genes.txt > $InterPro_dups
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/feature_extraction
$ProgDir/interpro2wego.py --inp_interpro $InterPro_dups > $OutDir/duplicated_genes_GO_WEGO.txt
