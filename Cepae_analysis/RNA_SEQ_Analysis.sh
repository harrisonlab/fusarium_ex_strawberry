# RNA-Seq analysis
​
## Perform qc on RNA-Seq data
​
```bash
    # Run fastq-mcf
    for RNADir in $(ls -d ../../../oldhome/groups/harrisonlab/project_files/fusarium/raw_rna/paired/F.oxysporum_fsp_cepae/* | grep 'CzapekDox\|PDA\|PDB\|GlucosePeptone'); do
        FileF=$(ls $RNADir/F/*1.fastq.gz | head -n $num | tail -n1)
        FileR=$(ls $RNADir/R/*1.fastq.gz | head -n $num | tail -n1)
        echo $FileF
        echo $FileR
        #Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1.fq.gz//g')
        #echo $Sample_Name
        Media=$(echo $RNADir | rev | cut -d '/' -f1 | rev | sed -r 's/Fus2_//g')
        echo $Media
        OutDir=RNAseq/Fusarium_media/$Media
        IluminaAdapters=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p himem $ProgDir/fastq-mcf_himem.sh $FileF $FileR $IluminaAdapters RNA $OutDir
    done
```
​
```bash
    # Run fastqc
    for RawData in $(ls RNAseq/*/*/*/*.fq.gz); do
        echo $RawData
        OutDir=$(dirname $RawData)
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p short $ProgDir/fastqc.sh $RawData $OutDir
    done
```
​
## Decontamination of rRNA reads in RNAseq data
​
```bash
    for RNADir in $(ls -d RNAseq/*/*); do
        FileF=$(ls $RNADir/F/*trim.fq.gz | head -n $num | tail -n1)
        FileR=$(ls $RNADir/R/*trim.fq.gz | head -n $num | tail -n1)
        echo $FileF
        echo $FileR
        Ref=/data/scratch/gomeza/prog/bbmap/ribokmers.fa.gz
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
        echo $RNADir
        Strain=$(echo $RNADir | rev | cut -d '/' -f1 | rev)
        echo $Strain
        sbatch -p himem $ProgDir/bbduk.sh $Ref "$RNADir"/cleaned $FileF $FileR $ProgDir $Strain
    done
```
​
## Salmon
​
​
```bash
conda activate salmon
​
    for Transcriptome in $(ls ../../../fusarium_ex_strawberry/NextGenSeq/Fo_Cepae/gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_combined.cdna.fasta); do
        Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
        Organism=$(echo $Transcriptome| rev | cut -d '/' -f4 | rev)
        echo "$Organism - $Strain"
        for RNADir in $(ls -d RNAseq/*/*/cleaned); do
            FileF=$(ls $RNADir/F/*cleaned.fq.gz | head -n $num | tail -n1)
            FileR=$(ls $RNADir/R/*cleaned.fq.gz | head -n $num | tail -n1)
            echo $FileF
            echo $FileR
            Media=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
            echo "$Media"
            OutDir=RNAseq/alignment/salmon/Fus2/$Media
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/RNAseq_analysis
            sbatch -p himem $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
        done
    done
```
​
Convert Salmon quasi-quanitifcations to gene counts using an awk script:
​
```bash
mkdir -p RNAseq/alignment/salmon/DeSeq2
#This command creates a two column file with transcript_id and gene_id.
    for File in $(ls RNAseq/alignment/salmon/*/*/quant.sf | head -n1); do
        cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > RNAseq/alignment/salmon/DeSeq2/trans2gene.txt
    done
​
#Put files in a convenient location for DeSeq.
​
    for File in $(ls RNAseq/alignment/salmon/*/*/quant.sf); do
        Prefix=$(echo $File | cut -f5 -d '/' --output-delimiter '_')
        mkdir -p RNAseq/alignment/salmon/DeSeq2/$Prefix
        cp $PWD/$File RNAseq/alignment/salmon/DeSeq2/$Prefix/quant.sf
        # rm RNAseq/alignment/salmon/DeSeq2/$Prefix/quant.sf
    done
```
​
# Import (with tximport) and extract
​
​
```bash
/projects/software/R-3.6.1/bin/R
```
​
```R
setwd("/projects/neonectria_ditissima/gomez_WD/akin")
​
# Load libraries
library(tximport)
​
#Import transcript to gene mapping info
tx2gene <- read.table("RNAseq/alignment/salmon/DeSeq2/trans2gene.txt",header=T,sep="\t")
​
#Import quantification files
txi.reps <- tximport(paste(list.dirs("RNAseq/alignment/salmon/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)
​
#Get the sample names from the folders
mysamples <- list.dirs("RNAseq/alignment/salmon/DeSeq2",full.names=F,recursive=F)
​
#Summarise to gene level. This can be done in the tximport step (txOut=F)
txi.genes <- summarizeToGene(txi.reps,tx2gene)
names(txi.genes)
​
#Set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))
​
write.table(txi.genes,"RNAseq/alignment/salmon/DeSeq2/txigenes.txt",sep="\t",na="",quote=F)
write.table(txi.reps,"RNAseq/alignment/salmon/DeSeq2/txireps.txt",sep="\t",na="",quote=F)
```
