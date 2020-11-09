#Results show this is best STAR script - it worked :)
# Run in Repenv - condaenv
#Fus2_CzapekDox, Fus2_GlucosePeptone, Fus2_PDA and Fus2_PDB are RNAseq data of infected onions
#Only data samples that will map genes of F.oxy accurately
#Need to concatenate data after STAR analysis

for Assembly in $(ls repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa)
  do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    FileF=qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/F/*_trim.fq.gz
    FileR=qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/R/*_trim.fq.gz
    echo $FileF
    echo $FileR
    Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
    echo "$Timepoint"
    Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
    OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
    Preindex=13
    ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
    sbatch $ProgDir/STAR_1.sh $Assembly $FileF $FileR $OutDir
  done
done
