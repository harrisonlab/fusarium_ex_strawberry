#Braker script - edit complete
#One for each concatenated file-V2
#Run in conda env (Repenv)
#AcceptedHits=alignment/concatenated.bam
#Alternate strain for softmasked and hits - make sure they are same
#Intial run required installation of Hash::Merge and Logger::Simple using cpan
#For DSA15 iso had to add numerical suffix to GeneMN for braker to run for CZAP
#Made mistake of running braker for each RNAseq data instead of concatenated data
#changed prog  dir to mine /home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction

for Assembly in $(ls repeat_masked/*/DSA14_003/ncbi_edits_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/braker/$Organism/$Strain/
    AcceptedHits=alignment/star/F.oxysporum_fsp_fragariae/DSA14_003/concatenated/concatenated.bam
    GeneModelName="$Organism"_"$Strain"_brakerV2
    ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/scripts
    sbatch $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
