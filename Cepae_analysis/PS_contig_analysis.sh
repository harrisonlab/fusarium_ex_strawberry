# stat analysis of contig_11
# how long is contig 11
# and is it enriched for mimps
# have a think about how you might test that...
# step 1 - isolate the contig

# create a txt file with contig name
    nano contig_12tmp.txt # type name of contig you want to excise

    faidx -d '|' pilon_10_renamed.fasta $(tr '\n' ' ' < contig_12tmp.txt ) > contig_12.fasta

# opened in word - 2,591,857 characters or bases
# file size in linux is 2591858

# Mimp analysis step
# Run in Repenv
# used full paths for input files

for Assembly in $(ls assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/pilon/contig_12.fasta); do #pilon_10_renamed.fasta
    Organism=$(echo "$Assembly" | rev | cut -d '/' -f4 | rev)
    Strain=$(echo "$Assembly" | rev | cut -d '/' -f3 | rev)
    OutDir=assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/pilon/mimps_contig_12
    mkdir -p "$OutDir"
    echo "$Organism - $Strain"
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
    $ProgDir/mimp_finder.pl $Assembly $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
    $ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
  done

  # There are 45 sequences that contain the consensus mimp motif in miniasm
  # 74 in whole assembly

#tried whole script
for Assembly in $(ls assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon/contig_11.fasta); do #pilon_10_renamed.fasta
    Organism=$(echo "$Assembly" | rev | cut -d '/' -f4 | rev)
    Strain=$(echo "$Assembly" | rev | cut -d '/' -f3 | rev)
    GeneGff=$(ls ../gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/final/final_genes_appended.gff3)
    OutDir=assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon/mimps_contig_11_V2
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
