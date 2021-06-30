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
mkdir temp_busco4
printf "" > analysis/popgen/busco_phylogeny_2/single_hits.txt
  for Busco in $(cat analysis/popgen/busco_phylogeny_2/all_buscos_*.txt); do
  OutDir=analysis/popgen/busco_phylogeny_2/$Busco
  mkdir -p $OutDir
  # Move all single copy genes to each folder and rename gene headers
    for Fasta in $(ls ../fusarium_EX_Lactucae/Assembly2/miniasm/F.oxysporum_fsp_lactucae/race_1/pilon/busco_sordariomycetes_obd10/round_*/pilon_10_renamed/run_sordariomycetes_odb10/busco_sequences/single_copy_busco_sequences/$Busco.fna); do
      Strain=race_1
      Organism=F.oxysporum_fsp_lactucae
      FileName=$(basename $Fasta)
      contig=$(cat $Fasta | grep '>' | sed 's/ <unknown description>//g' | sed 's/>//g')
      echo ">$Busco:$Strain:$contig" > temp_busco3/"$Busco"_"$Strain"_new_names.txt
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools
      python $ProgDir/replace_fasta_records.py -i $Fasta -r temp_busco4/"$Busco"_"$Strain"_new_names.txt -o $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
      rm temp_busco4/"$Busco"_"$Strain"_new_names.txt
    #cat $Fasta | sed "s/:.*.fasta:/:"$Organism"_"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
    done
  # Create fasta file containing all busco for alignment
  cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
  SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | wc -l)
  printf "$Busco\t$SingleBuscoNum\n" >> analysis/popgen/busco_phylogeny_2/single_hits.txt
  done
rm -r temp_busco2

# If all isolates have a single copy of a busco gene, move the appended fasta to a new folder
    OutDir=analysis/popgen/busco_phylogeny_2/alignments
    mkdir -p $OutDir
    OrganismNum=$(cat analysis/popgen/busco_phylogeny_2/single_hits.txt | cut -f2 | sort -nr | head -n1)
    for Busco in $(cat analysis/popgen/busco_phylogeny_2/all_buscos_*.txt); do
    echo $Busco
    HitNum=$(cat analysis/popgen/busco_phylogeny_2/single_hits.txt | grep "$Busco" | cut -f2)
    if [ $HitNum == $OrganismNum ]; then
      cp analysis/popgen/busco_phylogeny_2/$Busco/"$Busco"_appended.fasta $OutDir/.
    fi
    done

## Gene alignemtns
# Since you are logged into a partition, srun will run job there
    AlignDir=analysis/popgen/busco_phylogeny_2/alignments
    CurDir=$PWD
    cd $AlignDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Phylogenetics
    srun $ProgDir/mafft.sh
    cd $CurDir
===>
## Nucleotide diversity (optional)
For closely related organisms, identify genes with high nucleotide diversity (Pi) and average number of pairwise differences, medium number of segregating sites (avoid alignments with low homology and lots of phylogenetically uninformative singletons).
For analyses involving cross-species comparisons involving highly diverged sequences with high nucleotide diversity (e.g. [0.1<Pi<0.4]), looking for genes with the lowest number of segregating sites.

conda install -c bioconda dendropy

    cd analysis/popgen/busco_phylogeny_2/alignments
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Phylogenetics
    python $ProgDir/calculate_nucleotide_diversity.py "*aligned.fasta"

Trim poor alignments
# Edit header name keeping BUSCO name and isolate name using sed. Sample ID should contain at least 1 letter.
cd analysis/popgen/busco_phylogen1y_2/alignments
# e.g.
sed -i 's/:contig.*//g' *_appended_aligned.fasta
sed -i 's/:LD.*//g' *_appended_aligned.fasta
sed -i 's/:NODE.*//g' *_appended_aligned.fasta

Trimming sequence alignments using Trim-Al. Note - automated1 mode is optimised for ML tree reconstruction

  OutDir=analysis/popgen/busco_phylogeny_2/trimmed_alignments
  mkdir -p $OutDir
  for Alignment in $(ls analysis/popgen/busco_phylogeny_2/alignments/*_appended_aligned.fasta); do
    TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
    echo $Alignment
    trimal -in $Alignment -out $OutDir/$TrimmedName -keepheader -automated1
  done

## Randomized Axelerated Maximum Likelihood
  for Alignment in $(ls analysis/popgen/busco_phylogeny_2/trimmed_alignments/*aligned_trimmed.fasta); do
      #sleep option will submit a job every 10s to the short queue
      sleep 10s
      Prefix=$(basename $Alignment | cut -f1 -d '_')
      OutDir=analysis/popgen/busco_phylogeny_2/RAxML/$Prefix
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Phylogenetics
      srun $ProgDir/RAxML.sh $Alignment $Prefix $OutDir
  done

## Astral

Run Astral to build a consensus phylogeny from a collective set of "best phylogenies" from each BUSCO locus.
Note - "Recent versions of ASTRAL output a branch support value even without bootstrapping. Our analyses have revealed that this form of support is more reliable than bootstrapping (under the conditions we explored).
Nevertheless, you may want to run bootstrapping as well." Tutorial tips: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial-template.md
##screen -a
# Edit if necessary. This will run on compute10
##srun --partition long --mem-per-cpu 10G --cpus-per-task 24 --pty bash

OutDir=analysis/popgen/busco_phylogeny_2/ASTRAL
mkdir -p $OutDir

#Concatenate best trees
Name=Name4yourPhylogeny
#cat analysis_VP/popgen/busco_phylogeny_2/RAxML/*/RAxML_bestTree.*  | sed -r "s/CTG.\w+:/:/g" > $OutDir/Nd_phylogeny.appended2.tre
cat analysis/popgen/busco_phylogeny_2/RAxML/*/RAxML_bestTree.*  > $OutDir/"$Name"_phylogeny.appended.tre

# Contract low support brances (below 10% bootsrap support)
nw_ed $OutDir/"$Name"_phylogeny.appended.tre 'i & b<=10' o > $OutDir/"$Name"_phylogeny.appended.trimmed.tre

# Calculate combined tree
ProgDir=/scratch/software/ASTRAL/ASTRAL-5.7.1/Astral
java -jar $ProgDir/astral.5.7.1.jar -i $OutDir/"$Name"_phylogeny.appended.tre -o $OutDir/"$Name"_phylogeny.consensus.tre 2> $OutDir/"$Name"_phylogeny.consensus.log
# Score the resulting tree
java -jar $ProgDir/astral.5.7.1.jar -q $OutDir/"$Name"_phylogeny.consensus.tre -i $OutDir/"$Name"_phylogeny.appended.tre -o $OutDir/"$Name"_phylogeny.consensus.scored.tre 2> $OutDir/"$Name"_phylogeny.consensus.scored.log

Manual edition of the final consensus tree is needed

Step 1: Download consensus tree to local machine

Step 2: Import into geneious and export again in newick format to get around polytomy branches having no branch length.

Step 3 (Probably optional): Terminal branch lengths are meanlingless from ASTRAL and should all be set to an arbitrary value. This will be done by geneious (set to 1), but it also introduces a branch length of 2 for one isolate that needs to be corrected with sed

cat phylogeny.consensus.scored.geneious.tre | sed 's/:2/:1/g' > phylogeny.consensus.scored.geneious2.tre

Work in progress
Plot best scored tree
GGtree was used to make a plot. Tutorial tips: https://bioconnector.org/r-ggtree.html

R version > 4.0

setwd("/data/scratch/gomeza/")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)
library(ggtree)
library(phangorn)
library(treeio)





tree <- read.tree("/Users/armita/Downloads/Aalt/ASTRAL/expanded/Alt_phylogeny.consensus.scored.geneious2.tre")

t<-ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()

mydata <- read.csv("/Users/armita/Downloads/Aalt/ASTRAL/traits.csv", stringsAsFactors=FALSE)
rownames(mydata) <- mydata$Isolate
mydata <- mydata[match(tree$tip.label,rownames(mydata)),]


# Format nodes by values
nodes <- data.frame(t$data)
#nodes <- nodes[!nodes$isTip,]
nodes$label <- as.numeric(nodes[nodes$label,])
as.numeric(nodes$label)
#nodes$label[nodes$label < 0.80] <- ''
nodes$support[nodes$isTip] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label > 0.80)] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label < 0.80)] <- 'unsupported'
nodes$support[(!nodes$isTip) & (nodes$label == '')] <- 'supported'
t <- t + aes(linetype=nodes$support)
nodes$label[nodes$label > 0.80] <- ''
t <- t + geom_nodelab(data=nodes, size=2, hjust=-0.05) # colours as defined by col2rgb



#51 is the node of your outgroup?
tree$edge.length[tree$edge.length == 1] <- 0
tree$edge.length[51] <- 0


t <- ggtree(tree, aes(linetype=nodes$support)) # Core tree


# Adjust terminal branch lengths:
branches <- t$data

branches <- t$data
tree$edge.length[branches$isTip] <- 1.0
# tree$edge.length[tree$edge.length == 1] <- 0
# t <- ggtree(tree, aes(linetype=nodes$support))
#Tree <- branches$branch.length

t <- t + geom_treescale(offset=-1.0, fontsize = 3) # Add scalebar
# t <- t + xlim(0, 0.025) # Add more space for labels


# Colouring labels by values in another df
t <- t %<+% mydata # Allow colouring of nodes by another df
#t <- t + geom_tiplab(aes(color=Source), size=3, hjust=0) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

tips <- data.frame(t$data)
tips$label <- tips$ID
t <- t + geom_tiplab(data=tips, aes(color=Source), size=3, hjust=0, align=T, offset = +0.1) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

# Add in a further set of labels
tree_mod <- data.frame(t$data)
tree_mod$label <- tips$pathotype
t <- t + geom_tiplab(data=tree_mod, aes(label=label, color=Source), align=T, linetype = NULL, size=3, hjust=0, offset = +5.0) +
scale_color_manual(values=c("gray39","black"))

tips$MAT <- factor(tips$MAT)
# t <- t + geom_tippoint(data=tips, aes(shape=MAT), size=2)
t <- t + geom_tiplab(data=tips, aes(label=MAT, color=Source), align=T, linetype = NULL, size=3, hjust=0, offset = +3.5) +
scale_color_manual(values=c("gray39","black"))



# Annotate a clade with a bar line
# t <- t + geom_cladelabel(node=42, label='sect. Alternaria', align=T, colour='black', offset=-1.5)
# t <- t + geom_cladelabel(node=70, label='gaisen clade', align=T, colour='black', offset=-4.5)
# t <- t + geom_cladelabel(node=51, label='tenuissima clade', align=T, colour='black', offset=-4.5)
# t <- t + geom_cladelabel(node=45, label='arborescens clade', align=T, colour='black', offset=-4.5)
t <- t + geom_cladelabel(node=43, label='sect. Alternaria', align=T, colour='black', offset=9.5)
t <- t + geom_cladelabel(node=70, label='gaisen clade', align=T, colour='black', offset=6.5)
t <- t + geom_cladelabel(node=46, label='tenuissima clade', align=T, colour='black', offset=6.5)
t <- t + geom_cladelabel(node=65, label='arborescens clade', align=T, colour='black', offset=6.5)
t <- t + geom_cladelabel(node=65, label='', colour='NA', offset=17.5)

# Save as PDF and force a 'huge' size plot
# t <- ggsave("expanded/Fig3_busco_phylogeny.pdf", width =30, height = 30, units = "cm", limitsize = FALSE)
t <- ggsave("expanded/Fig3_busco_phylogeny.tiff", width =30, height = 30, units = "cm", limitsize = FALSE)
Visually inspect the alignments of selected genes (genes_selected_for_phylogeny.txt) to be used in
constructing the phylogenies and trim them as necessary in MEGA7.
Copy the relevant trimmed alignment FASTA files into

```bash
  # mkdir $CurDir/beast_runs/candidates/select/trimmed

  ##PartitionFinder (nucleotide sequence evolution model)

cd analysis/popgen/busco_phylogeny/phylogeny

config_template=/home/sobczm/bin/PartitionFinder1.1.1/partition_finder.cfg
ct=$(basename "$config_template")

mkdir NEXUS

# prepare directory for PartitionFinder run:
for f in $(ls *fasta); do
sed -i 's/:/_/g' $f
c="$(cat $f | awk 'NR%2==0' | awk '{print length($1)}' | head -1)"
p="${f%.fasta}.phy"
n="${f%.fasta}.NEXUS"
dir="${f%.fasta}"

mkdir $dir
cp $config_template $dir/.

# Substitute the name of the alignment file and the sequence length in the config file to become correct for the current run.
sed -i 's,^\(alignment = \).*,\1'"$p;"',' $dir/$ct
sed -i 's,^\(Gene1_pos1 = \).*,\1'"1-$c\\\3;"',' $dir/$ct
sed -i 's,^\(Gene1_pos2 = \).*,\1'"2-$c\\\3;"',' $dir/$ct
sed -i 's,^\(Gene1_pos3 = \).*,\1'"3-$c\\\3;"',' $dir/$ct

# Convert FASTA to phylip for the Partition Finder run
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
$ProgDir/fasta2phylip.pl $f>$p
mv $p $dir

# Convert FASTA to NEXUS for the BEAST run
$ProgDir/Fasta2Nexus.pl $f>$dir/$n

#Problems running PartitionFinder on the cluster. May have to be run locally on your Mac or Windows machine.
# qsub $ProgDir/sub_partition_finder.sh $dir
done
