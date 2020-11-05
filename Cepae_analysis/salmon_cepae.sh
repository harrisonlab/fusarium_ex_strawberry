# Run in Salmon envm (salmon)
# This script simply loops through each sample and invokes salmon
# Fus2_CzapekDox, Fus2_GlucosePeptone, Fus2_PDA and Fus2_PDB are RNAseq data of infected onions

for Transcriptome in $(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_combined.cdna.fasta); do
        Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
        Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
        echo "$Organism - $Strain"
        for RNADir in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox); do
        FileF=$(ls $RNADir/F/6_S2_L001_R1_001_trim.fq.gz)
        FileR=$(ls $RNADir/R/6_S2_L001_R2_001_trim.fq.gz)
        echo $FileF
        echo $FileR
        Sample_Name=Fus2_CzapekDox
        echo $Sample_Name
        OutDir=RNAseq_analysis/salmon/$Organism/$Strain/$Sample_Name
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
        sbatch $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
        done
    done

    # qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/F/     7_S3_L001_R1_001_trim.fq.gz
    # qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/R/      7_S3_L001_R2_001_trim.fq.gz
    # qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/F/       4_S1_L001_R1_001_trim.fq.gz
    # qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/R/       4_S1_L001_R2_001_trim.fq.gz
    # qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/F/       9_S4_L001_R1_001_trim.fq.gz
    #  qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/R/      9_S4_L001_R2_001_trim.fq.gz
