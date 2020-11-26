# Aligning illumina reads against pilon data
# Run in conda env
# Make sure script is executable
# chmod u+x ./sub_pilon.sh
#ran in a different folder to test pilon, because the path are different

```bash
for Assembly in $(ls race_1_smartdenovo_racon_round_10_renamed.fasta); do
  Organism=F.oxysporum_fsp_lactucae
  Strain=AJ520
  IlluminaDir=$(ls -d $Strain)
  echo $Strain
  echo $Organism
  TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n2 | tail -n1);
  TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n2 | tail -n1);
  echo $TrimF1_Read
  echo $TrimR1_Read
  OutDir=$(dirname $Assembly)
  Iterations=10
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
  sbatch $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done

  # You might rename your contigs at this point using remove_contaminants.py

  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
      touch tmp.txt
      for Assembly in $(ls assembly/...); do
          OutDir=$(dirname $Assembly)
          $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/race_1_smartdenovo_racon_round_3_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
      done
      rm tmp.txt

#for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/Pilon_SDen/race_1_smartdenovo_racon_round_10_renamed.fasta); do
#  FileF=$(ls ../oldhome/groups/harrisonlab/project_files/fusarium/raw_dna/paired/F.oxysporum_fsp_lactucae/AJ520/F/AJ520_S2_L001_R1_001.fastq.gz)
#  FileR=$(ls ../oldhome/groups/harrisonlab/project_files/fusarium/raw_dna/paired/F.oxysporum_fsp_lactucae/AJ520/R/AJ520_S2_L001_R2_001.fastq.gz)
# OutDir=assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/Pilon_SDen/
# Iterations=10
# ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
# sbatch $ProgDir/sub_pilon.sh $Assembly $FileF $FileR $OutDir $Iterations
# done

# For fragariae
#/projects/fusarium_ex_strawberry/raw_dna/paired/F.oxysporum_fsp_fragariae/DSA14_003/F/14-003_S4_L001_R1_001.fastq.gz
# /projects/fusarium_ex_strawberry/raw_dna/paired/F.oxysporum_fsp_fragariae/DSA14_003/R/14-003_S4_L001_R2_001.fastq.gz
