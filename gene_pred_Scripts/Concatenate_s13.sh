#Need to do this step to link RNAseq data into one series
#Concatenate- to link (things) together in a chain or series
#Missed this step therefore did Braker a almost 6 different times

Strain="DSA14_003"
  Organism="F.oxysporum_fsp_fragariae"
  mkdir -p alignment/star/$Organism/$Strain/concatenated
  samtools merge -f alignment/star/$Organism/$Strain/concatenated/concatenated.bam \
  alignment/star/$Organism/$Strain/Fus2_CzapekDox/6_S2_L001_R1_001_trim.fq.gz/star_aligmentAligned.sorted.out.bam \
  alignment/star/$Organism/$Strain/Fus2_GlucosePeptone/7_S3_L001_R1_001_trim.fq.gz/star_aligmentAligned.sorted.out.bam \
  alignment/star/$Organism/$Strain/Fus2_PDA/9_S4_L001_R1_001_trim.fq.gz/star_aligmentAligned.sorted.out.bam \
  alignment/star/$Organism/$Strain/Fus2_PDB/4_S1_L001_R1_001_trim.fq.gz/star_aligmentAligned.sorted.out.bam

#new star otput - star_aligmentAligned.sortedByCoord.out.bam
