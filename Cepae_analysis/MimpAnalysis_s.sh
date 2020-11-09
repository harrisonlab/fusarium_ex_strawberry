# MIMP analysis step
# Run in Repenv
# used full paths for input files
# Use diff out directory - non writable
# You messed up on last script run you used repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_softmasked_repeatmasker_TPSI_appended.fa
# Should have used repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa

for Assembly in $(ls repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa); do
    Organism=$(echo "$Assembly" | rev | cut -d '/' -f4 | rev)
    Strain=$(echo "$Assembly" | rev | cut -d '/' -f3 | rev)
    GeneGff=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended_renamed.gff3)
    OutDir=../../../../../fusarium_ex_strawberry/F.oxysporum_fsp_cepae/Fus2_canu_new/mimps/V3
    mkdir -p "$OutDir"
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
    $ProgDir/mimp_finder.pl $Assembly $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
    $ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
    echo "The number of mimps identified:"
    cat $OutDir/"$Strain"_mimps.fa | grep '>' | wc -l
    bedtools intersect -u -a $GeneGff -b $OutDir/"$Strain"_mimps_exp.gff > $OutDir/"$Strain"_genes_in_2kb_mimp.gff
    echo "The following transcripts intersect mimps:"
    MimpProtsTxt=$OutDir/"$Strain"_prots_in_2kb_mimp.txt
    MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
    cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | sort | uniq > $MimpProtsTxt
    cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | cut -f1 -d '.'| sort | uniq > $MimpGenesTxt
    cat $MimpProtsTxt | wc -l
    cat $MimpGenesTxt | wc -l
    echo ""
  done
