# Add path to profile PATH=${PATH}:/projects/oldhome/armita/prog/cufflinks/cufflinks-2.2.1.Linux_x86_64
# update profile - . ~/.profile
# Run in a screen and a node in the
# srun --partition himem --time 0-06:00:00 --mem 20G --cpus-per-task 24 --pty bash
# run for each star alignment
# Alignment = star alignments output files
# OutDir Ref genome gff3
# alignment/star/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_PDB/4_S1_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam
# alignment/star/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_PDA/9_S4_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam
# alignment/star/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_CzapekDox/6_S2_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam
# alignment/star/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_GlucosePeptone/7_S3_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam


    for Alignment in $(ls alignment/star/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_PDA/9_S4_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam); do
    Gff=final/final_genes_appended.gff3
    OutDir=$(dirname $Alignment)
    mkdir -p $OutDir/fpkm
    cufflinks -p 8 -o $OutDir/fpkm -G $Gff $Alignment
  done
