#Additional transcripts - to be edited
#Run in perly env (Repenv)
#Type full paths, do not use asterisks
#Run each line indivdually otherwise it will not work
#Do segments one at a time for peace of mind

BrakerGff=$(ls -d gene_pred/braker/F.oxysporum_fsp_fragariae/DSA15_041/F.oxysporum_fsp_fragariae_DSA15_041_brakerV2/augustus.gff3)
	Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev)
	Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
	echo "$Organism - $Strain"
	Assembly=$(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA15_041/ncbi_edits_repmask/DSA15_041_contigs_softmasked_repeatmasker_TPSI_appended.fa)
	CodingQuarryGff=gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA15_041/out/PredictedPass.gff3
	PGNGff=gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA15_041/out/PGN_predictedPass.gff3
	AddDir=gene_pred/codingquary/$Organism/$Strain/additional
	FinalDir=gene_pred/codingquary/$Organism/$Strain/final
	AddGenesList=$AddDir/additional_genes.txt
	AddGenesGff=$AddDir/additional_genes.gff
	FinalGff=$AddDir/combined_genes.gff
	mkdir -p $AddDir
	mkdir -p $FinalDir

# Create a list with the additional transcripts in CondingQuarry gff (and CQPM) vs Braker gene models
#For first line had to put direct paths for -a and -b
	bedtools intersect -v -a $CodingQuarryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
	bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList

# Creat Gff file with the additional transcripts
	ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
	$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuarryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
	$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff

# Create a final Gff file with gene features
	$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3

  # Create fasta files from each gene feature in the CodingQuarry gff3
	$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary

  # Create fasta files from each gene feature in the Braker gff3
	cp $BrakerGff $FinalDir/final_genes_Braker.gff3
  $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker

  # Combine both fasta files
	cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
	cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
	cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
	cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

  # Combine both gff3 files
	GffBraker=$FinalDir/final_genes_CodingQuary.gff3
	GffQuary=$FinalDir/final_genes_Braker.gff3
	GffAppended=$FinalDir/final_genes_appended.gff3
	cat $GffBraker $GffQuary > $GffAppended

  # Check the final number of genes

	for DirPath in $(ls -d $FinalDir); do
    echo $DirPath;
    cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
    echo "";
	done
