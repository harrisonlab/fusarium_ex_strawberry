for Assembly in $(ls busco_phy/4287_chromosomal_contigs_unmasked.fa); do
  Strain=4287
  Organism=F.oxysporum_fsp_lycopersici
  echo "$Organism - $Strain"
  ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/quality_check
  BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
  OutDir=busco_phy/$Strain/busco_sordariomyceta_odb9
  sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
done

15-074_contigs_unmasked.fa
 4287_chromosomal_contigs_unmasked.fa
 DSA14_003_contigs_unmasked.fa
DSA15_041_contigs_unmasked.fa
 Fus2_canu_new_contigs_unmasked.fa
Straw465_contigs_unmasked.fa
