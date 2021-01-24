# need to get gene predictions for each media
# Check if STAR has been run
# Run Braker then cquary

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

  for Assembly in $(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_contigs_unmasked.fa); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
      Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev)
      echo "$Organism - $Strain"
      OutDir=gene_pred/codingquary/$Organism/$Strain/flye
      mkdir -p $OutDir
      GTF=gene_pred/stringtie/F.oxysporum_fsp_fragariae/DSA14_003/flye/concatenated_prelim/out.gtf
      ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/ProgScripts
      sbatch $ProgDir/codingquarry2.sh $Assembly $GTF $OutDir
    done


## Add gene prediction transcripts together

Additional transcripts - to be edited
Run in perly env (Repenv)
Type full paths, do not use asterisks
RUN LINE BY LINE AS IT WILL NOT WORK
Do segments one at a time for peace of mind

    BrakerGff=$(ls -d gene_pred/braker/F.oxysporum_fsp_fragariae/DSA14_003/flye/F.oxysporum_fsp_fragariae_DSA14_003_braker_flye_V2/augustus.gff3)
    	Strain=$(echo $BrakerGff| rev | cut -d '/' -f4 | rev)
    	Organism=$(echo $BrakerGff | rev | cut -d '/' -f5 | rev)
    	echo "$Organism - $Strain"
    	Assembly=$(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    	CodingQuarryGff=gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/flye/out/PredictedPass.gff3
    	PGNGff=gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/flye/out/PGN_predictedPass.gff3
    	AddDir=gene_pred/codingquary/$Organism/$Strain/flye/additional
    	FinalDir=gene_pred/codingquary/$Organism/$Strain/flye/final
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

For flye assembly Braker genes: 17980, CQ: 1103 & combined: 19083


# Genome annotations

## 1) Interproscan

    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
      for Genes in $(ls gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/flye/final/final_genes_appended_renamed.pep.fasta); do
        echo $Genes
        $ProgDir/interproscan.sh $Genes
      done 2>&1 | tee -a interproscan_submission.log

Interproscan: all jobs failed - couldn't run all jobs simultaneously
ERROR: uk.ac.ebi.interpro.scan.management.model.implementations.RunBinaryStep - Command line failed with exit code: 1
Need to run in batches - Had to run split DNA in sets of 10.
Split gene.pep.fasta like so:
  InFile=gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/flye/final/final_genes_appended_renamed.pep.fasta
  SplitDir=gene_pred/interproscan/$Organism/$Strain/flye
  InterproDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  InName=$(basename $InFile)
  mkdir -p $SplitDir
  $InterproDir/splitfile_500.py --inp_fasta $InFile --out_dir $SplitDir --out_base "$InName"_split

  for file in $(ls gene_pred/interproscan/F.oxysporum_fsp_fragariae/DSA14_003/flye/*_split_9*); do
    sbatch /home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation/run_interproscan.sh $file
    done

Need to merge interproscan output as follows

    ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/Feature_annotation
     for Proteins in $(ls gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/flye/final/final_genes_appended_renamed.pep.fasta); do
       Strain=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
       Organism=$(echo $Proteins | rev | cut -d '/' -f5 | rev)
       echo "$Organism - $Strain"
       echo $Strain
       InterProRaw=gene_pred/interproscan/F.oxysporum_fsp_fragariae/DSA14_003/flye/raw
       $ProgDir/append_interpro.sh $Proteins $InterProRaw
     done

Use this command to view particular features in interproscan data:
  less path/to/interproscan.tsv | grep 'gene feature' # e.g. transposon
