# need to get gene predictions for each media
# Check if STAR has been run
# Run Braker then cquary

# Mask it all

conda activate RMask
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking
    BestAssembly=../../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa
    OutDir=Fo_Cepae/repeat_masked
    sbatch $ProgDir/rep_modeling.sh $BestAssembly $OutDir
    sbatch $ProgDir/transposonPSI.sh $BestAssembly $OutDir

## STAR

Run in Repenv - condaenv
Fus2_CzapekDox, Fus2_GlucosePeptone, Fus2_PDA and Fus2_PDB are RNAseq data of infected onions
Only data samples that will map genes of F.oxy accurately
Need to concatenate data after STAR analysis
#--genomeSAindexNbases is unique to each genome and is 11 for FoFR

  for Assembly in $(ls ../../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa);  do
      Strain=Fus2_canu_new
      Organism=F.oxysporum_fsp_cepae
      echo "$Organism - $Strain"
      FileF=../../oldhome/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/F/*_trim.fq.gz
      FileR=../../oldhome/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/R/*_trim.fq.gz
      echo $FileF
      echo $FileR
      Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
      echo "$Timepoint"
      Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
      OutDir=Fo_Cepae/alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/Genome_alignment
      sbatch $ProgDir/STAR_1.sh $Assembly $FileF $FileR $OutDir 11
    done


## Braker

Run in conda env (Repenv)
AcceptedHits=alignment/concatenated.bam
Alternate strain for softmasked
Intial run required installation of Hash::Merge and Logger::Simple using cpan


    #Original prog  dir /home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction

    for Assembly in $(ls /projects/oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
        Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
        echo "$Organism - $Strain"
        OutDir=gene_pred/braker/$Organism/$Strain/flye
        AcceptedHits=alignment/star/F.oxysporum_fsp_fragariae/DSA14_003/concatenated/concatenated.bam
        GeneModelName="$Organism"_"$Strain"_braker_flye_V2
        ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/ProgScripts
        sbatch $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
      done


  ## StringTie

  String tie - to be edited
  Run in conda env with Python 2.7 (betaenv)
  Codingquarry is another tool for gene prediction that it is able to predict additional genes in fungi
  Merge with Braker to give final gene model set
  /home/akinya/git_repos/assembly_fusarium_ex/scripts
  /home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
  # AcceptedHits lists
#Fo_Cepae/alignment/star/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_CzapekDox/6_S2_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam
#Fo_Cepae/alignment/star/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_GlucosePeptone/7_S3_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam
#Fo_Cepae/alignment/star/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_PDA/9_S4_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam
#Fo_Cepae/alignment/star/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_PDB/4_S1_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam

      for Assembly in $(ls ../../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa); do
          Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
          Organism=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
          echo "$Organism - $Strain"
          OutDir=Fo_Cepae/gene_pred/stringtie/$Organism/$Strain/PDB
          mkdir -p $OutDir
          AcceptedHits=Fo_Cepae/alignment/star/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_PDB/4_S1_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam
          ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/ProgScripts/Genome_alignment
          sbatch $ProgDir/stringtie.sh $AcceptedHits $OutDir
         done


## Codingquarry

#To be edited - completed
Run in env with Python 2.7 (betaenv)
After first run, use cquarryV1
GFT file from stringtie/cufflinks output
my repo /home/akinya/git_repos/assembly_fusarium_ex/scripts
Antonio /home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction

  for Assembly in $(ls ../../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
      Organism=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
      echo "$Organism - $Strain"
      OutDir=Fo_Cepae/gene_pred/codingquary/$Organism/$Strain/Gp
      mkdir -p $OutDir
      GTF=Fo_Cepae/gene_pred/stringtie/F.oxysporum_fsp_cepae/Fus2_canu_new/Gp/out.gtf
      ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/ProgScripts
      sbatch $ProgDir/codingquarry2.sh $Assembly $GTF $OutDir
    done

>>>
## Add gene prediction transcripts together

Additional transcripts - to be edited
Run in perly env (Repenv)
Type full paths, do not use asterisks
RUN LINE BY LINE AS IT WILL NOT WORK
Do segments one at a time for peace of mind

    BrakerGff=$(ls -d Fo_Cepae/gene_pred/braker/F.oxysporum_fsp_cepae/Fus2_canu_new/Czap/F.oxysporum_fsp_cepae_Fus2_canu_new_braker/augustus.gff3)
    	Strain=$(echo $BrakerGff| rev | cut -d '/' -f4 | rev)
    	Organism=$(echo $BrakerGff | rev | cut -d '/' -f5 | rev)
    	echo "$Organism - $Strain"
    	Assembly=$(ls /projects/oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    	CodingQuarryGff=Fo_Cepae/gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_canu_new/Czap/out/PredictedPass.gff3
    	PGNGff=Fo_Cepae/gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_canu_new/Czap/out/PredictedPass.gff3
    	AddDir=Fo_Cepae/gene_pred/codingquary/$Organism/$Strain/Czap/additional
    	FinalDir=Fo_Cepae/gene_pred/codingquary/$Organism/$Strain/Czap/final
    	AddGenesList=$AddDir/additional_genes.txt
    	AddGenesGff=$AddDir/additional_genes.gff
    	FinalGff=$AddDir/combined_genes.gff
    	mkdir -p $AddDir
    	mkdir -p $FinalDir

Create a list with the additional transcripts in CondingQuarry gff (and CQPM) vs Braker gene models
For first line had to put direct paths for -a and -b

  	bedtools intersect -v -a $CodingQuarryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
  	bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList

Creat Gff file with the additional transcripts

  	ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
  	$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuarryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
  	$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff

Create a final Gff file with gene features
    $ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3

Got this error:
  Possible precedence issue with control flow operator at /home/gomeza/miniconda3/envs/perly_env/lib/site_perl/5.26.2/Bio/DB/IndexedBase.pm line 805.

Create fasta files from each gene feature in the CodingQuarry gff3
    $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary

Got this error again (got it again after creating fasta files in braker gff3):
      Possible precedence issue with control flow operator at /home/gomeza/miniconda3/envs/perly_env/lib/site_perl/5.26.2/Bio/DB/IndexedBase.pm line 805.

Create fasta files from each gene feature in the Braker gff3
    cp $BrakerGff $FinalDir/final_genes_Braker.gff3
    $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker

Combine both fasta files
    cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
    cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
    cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
    cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

Combine both gff3 files
    GffBraker=$FinalDir/final_genes_CodingQuary.gff3
    GffQuary=$FinalDir/final_genes_Braker.gff3
    GffAppended=$FinalDir/final_genes_appended.gff3
    cat $GffBraker $GffQuary > $GffAppended

Check the final number of genes
  	for DirPath in $(ls -d $FinalDir); do
      echo $DirPath;
      cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
      cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
      cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
      echo "";
  	done

For Czap assembly Braker genes: 17387, CQ: 2030 & combined: 19417
Gp Braker: 17387, Cq : 2184 & comb: 19571
PDA 17438, CQ 1374, COMB 18812
PDB 17271, CQ 1484, COMB 18755

## Gene renaming
Rename genes too
Run line by line
Run in conda env (Repenv)

    #Remove duplicate and rename genes
    GffAppended=$(ls -d gene_pred/codingquary/F.oxysporum_fsp_lactucae/race_1/final/final_genes_appended.gff3)
    Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    FinalDir=gene_pred/codingquary/F.oxysporum_fsp_lactucae/race_1/final

    #Remove duplicated genes
    GffFiltered=$FinalDir/filtered_duplicates.gff
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered

    #Rename genes
    GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
    LogFile=$FinalDir/final_genes_appended_renamed.log
    $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
    rm $GffFiltered

    #Create renamed fasta files from each gene feature
    Assembly=$(ls repeat_masked/F.oxysporum_fsp_lactucae/race_1/flye/ncbi_edits_repmask/race_1_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    $ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_appended_renamed
    #The proteins fasta file contains * instead of Xs for stop codons, these should be changed
    sed -i 's/\*/X/g' $FinalDir/final_genes_appended_renamed.pep.fasta

    view gene names
    cat $FinalDir/final_genes_appended_renamed.cdna.fasta | grep '>'

# Genome annotations

## 1) Interproscan
# DO NOT RUN IN CONDA ENV
# This command will split your gene fasta file and run multiple interproscan jobs.
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
for Genes in $(ls Fo_Cepae/gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_canu_new/Czap/final/final_genes_combined.pep.fasta); do
echo $Genes
$ProgDir/interproscan.sh $Genes
done 2>&1 | tee -a interproscan_submisison_Czap.log


Need to merge interproscan output as follows

#gene_pred/interproscan/Fus2_canu_new/Czap/raw/final_genes_combined.pep.fasta

    ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/Feature_annotation
     for Proteins in $(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_canu_new/Czap/final/final_genes_combined.pep.fasta); do
       Media=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
       Strain=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
       Organism=$(echo $Proteins | rev | cut -d '/' -f5 | rev)
       echo "$Organism - $Strain"
       echo $Strain
       InterProRaw=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Czap/raw/
       $ProgDir/append_interpro_Cep.sh $Proteins $InterProRaw
     done

Use this command to view particular features in interproscan data:
  less path/to/interproscan.tsv | grep 'gene feature' # e.g. transposon
