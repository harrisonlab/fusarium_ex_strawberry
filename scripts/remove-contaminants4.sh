ProgDir=$(ls -d /home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants)
  touch tmp.csv
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
    mkdir $OutDir
    $ProgDir/remove_contaminants2.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
