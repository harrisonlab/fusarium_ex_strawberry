#String tie - to be edited
#One for each RNA seq data
#Run in conda env with Python 2.7 (betaenv)
#Codingquarry is another tool for gene prediction that it is able to predict additional genes in fungi
#Merge with Braker to give final gene model set
#Copied stringtie to my dir and made it executable
#/home/akinya/git_repos/assembly_fusarium_ex/scripts
#/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction

for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/stringtie/$Organism/$Strain/concatenated_prelim
    mkdir -p $OutDir
    AcceptedHits=alignment/star/F.oxysporum_fsp_fragariae/DSA1*/concatenated/concatenated.bam
    ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/scripts
    sbatch $ProgDir/stringtie.sh $AcceptedHits $OutDir
   done
