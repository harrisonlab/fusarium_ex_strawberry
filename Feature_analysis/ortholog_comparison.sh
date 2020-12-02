# Use to extract effectors from a gene.fasta file using names in txt file
# Create text file with known gene names of effectors/mimps you want to remove
# i.e  input gene names in Fof14_genes.txt
# Command uses gene name to copy across fasta seq to Fof genes

# Extract genes from cepae genome NOT FRAGARIAE

# faidx -d '|' final_genes_appended_renamed.cdna.fasta $(tr '\n' ' ' < Fof14_genes.txt ) > Fof14_genes.fasta - example command to test

# /projects/fusarium_ex_strawberry/gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/final/
# Why are you doing this?
# To compare candidate effectors in Fo cepae (which has RNA seq data) against the predicted genes in Fof using cepae RNA seq data

faidx -d '|' final_genes_combined.cdna.fasta $(tr '\n' ' ' < FocFus2_genes.txt ) > eff_ortho_genes.fasta # beta test - don't use

# faidx -d '|' final_genes_combined.cdna.fasta $(tr '\n' ' ' < FoC_cand_mimps.txt ) > Eff_mimp_genes.fasta

# faidx -d '|' final_genes_combined.cdna.fasta $(tr '\n' ' ' < FoC_Six.txt ) > six_ortho_genes.fasta

# Now that you have the cand effector genes, contrast against Fo_fragariae_14_003 long read seq genome
# Run in conda env with perly (Repenv)
# for $Assembly Use files with nucleotides

for Assembly in $(ls assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  Query=../F.oxysporum_fsp_cepae/Fus2_canu_new/final/Eff_mimp_genes.fasta # six_ortho_genes.fasta
  OutDir=assembly/flye/$Organism/$Strain/Orthology/FoFrvsFoCep_mimps
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
  sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly $OutDir
done

# second query../F.oxysporum_fsp_cepae/Fus2_canu_new/final/Eff_mimp_genes.fasta
# Use cut -f1 or * DSA14_003_eff_ortho_genes.fasta_homologs.csv to excise and view column data
# Compare against Andy's known SIX genes
for Assembly in $(ls assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  Query=../../oldhome/groups/harrisonlab/project_files/fusarium/analysis/blast_homology/six_genes/six-appended_parsed.fa
  OutDir=assembly/miniasm/$Organism/$Strain/Orthology
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
  sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly $OutDir
done

#############

# To see if miniasm assembly is chimeric, contrast other assemblies
# They haven't been pilon polished - currently cleaned with medaka
# Use pilon to polish these assemblies
# Use the diff queries ../F.oxysporum_fsp_cepae/Fus2_canu_new/final/Eff_mimp_genes.fasta # six_ortho_genes.fasta
# ../../oldhome/groups/harrisonlab/project_files/fusarium/analysis/blast_homology/six_genes/six-appended_parsed.fa

for Assembly in $(ls assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  Query=../../oldhome/groups/harrisonlab/project_files/fusarium/analysis/blast_homology/six_genes/six-appended_parsed.fa
  OutDir=assembly/flye/$Organism/$Strain/Orthology/FoFrvsFoCep
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
  sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly $OutDir
done

# Then do for flye

for Assembly in $(ls assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/medaka/assembly_racon_round_3_renamed.fasta); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  Query=../../oldhome/groups/harrisonlab/project_files/fusarium/analysis/blast_homology/six_genes/six-appended_parsed.fa
  OutDir=assembly/flye/$Organism/$Strain/Orthology
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
  sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly $OutDir
done

############

# BLAST DOES NOT DO WHOLE GENOME ANALYSIS
# query = assembly that needs checking


##########

# DO a mimp analysis and store data on a spreadsheet

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
