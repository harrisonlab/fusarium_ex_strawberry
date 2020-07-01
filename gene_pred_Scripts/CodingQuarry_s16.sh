#To be edited - completed
#Run in env with Python 2.7 (betaenv)
#After first run, use cquarryV1
#GFT file from stringtie/cufflinks output
#my repo /home/akinya/git_repos/assembly_fusarium_ex/scripts
#Antonio /home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction

for Assembly in $(ls repeat_masked/*/DSA15_041/ncbi_edits_repmask/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary/$Organism/$Strain/
    mkdir -p $OutDir
    GTF=gene_pred/stringtie/F.oxysporum_fsp_fragariae/DSA15_041/concatenated_prelim/*.gtf
    ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/scripts
    sbatch $ProgDir/codingquarry.sh $Assembly $GTF $OutDir
  done
