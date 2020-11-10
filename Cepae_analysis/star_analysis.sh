#Results show this is best STAR script - it worked for F.oxysporum_fsp_cepae:)
# Run in Repenv - condaenv
#Fus2_CzapekDox, Fus2_GlucosePeptone, Fus2_PDA and Fus2_PDB are RNAseq data of infected onions
#Only data samples that will map genes of F.oxy accurately
#Need to concatenate data after STAR analysis
#WARNING: --genomeSAindexNbases 13 is too large for the genome size=58458112, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 11

for Assembly in $(ls repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa)
  do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    FileF=qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/F/*_trim.fq.gz
    FileR=qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/R/*_trim.fq.gz
    echo $FileF
    echo $FileR
    Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
    echo "$Timepoint"
    Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
    OutDir=projects/fusarium_ex_strawberry/F.oxysporum_fsp_cepae/Fus2_canu_new/alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
    mkdir -p "$OutDir"
    Preindex=11
    ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
    sbatch $ProgDir/STAR_1.sh $Assembly $FileF $FileR $OutDir $Preindex
  done
