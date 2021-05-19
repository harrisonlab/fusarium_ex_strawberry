conda create -n geenral_tools

conda activate general_tools
# Gene aligner. Version=v7.475
conda install -c bioconda mafft
# Trim poor aligned sequences. Version=1.4.1
conda install -c bioconda trimal
# Randomized Axelerated Maximum Likelihood. Version=8.2.12
conda install -c bioconda raxml

#Create list of busco IDs
OutDir=analysis/popgen/busco_phylogeny
mkdir -p $OutDir
BuscoDb="sordariomycetes_odb10"
ls -1 /projects/dbBusco/$BuscoDb/hmms/*hmm | rev | cut -f1 -d '/' | rev | sed -e 's/.hmm//' > $OutDir/all_buscos_"$BuscoDb".txt

# For BUSCO version 4
screen -a
srun --partition long --mem-per-cpu 10G --cpus-per-task 10 --pty bash
# Create a folder for each busco gene
mkdir temp_busco
printf "" > analysis/popgen/busco_phylogeny/single_hits.txt
  for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
  echo $Busco
  OutDir=analysis_VP/popgen/busco_phylogeny_2/$Busco
  mkdir -p $OutDir
  # Move all single copy genes to each folder and rename gene headers
    for Fasta in $(ls gene_pred/busco/$Organism/$Strain/*/*/*/*/single_copy_busco_sequences/$Busco*.fna); do
      Strain=$(echo $Fasta | rev | cut -f7 -d '/' | rev)
      Organism=$(echo $Fasta | rev | cut -f8 -d '/' | rev)
      FileName=$(basename $Fasta)
      contig=$(cat $Fasta | grep '>' | sed 's/ <unknown description>//g' | sed 's/>//g')
      echo ">$Busco:$Strain:$contig" > temp_busco/"$Busco"_"$Strain"_new_names.txt
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools
      python $ProgDir/replace_fasta_records.py -i $Fasta -r temp_busco/"$Busco"_"$Strain"_new_names.txt -o $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
      rm temp_busco/"$Busco"_"$Strain"_new_names.txt
    #cat $Fasta | sed "s/:.*.fasta:/:"$Organism"_"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
    done
  # Create fasta file containing all busco for alignment
  cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
  SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | wc -l)
  printf "$Busco\t$SingleBuscoNum\n" >> analysis/popgen/busco_phylogeny/single_hits.txt
  done
rm -r temp_busco

# If all isolates have a single copy of a busco gene, move the appended fasta to a new folder
OutDir=analysis/popgen/busco_phylogeny/alignments
mkdir -p $OutDir
OrganismNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | cut -f2 | sort -nr | head -n1)
for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
echo $Busco
HitNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | grep "$Busco" | cut -f2)
if [ $HitNum == $OrganismNum ]; then
  cp analysis/popgen/busco_phylogeny/$Busco/"$Busco"_appended.fasta $OutDir/.
fi
done
