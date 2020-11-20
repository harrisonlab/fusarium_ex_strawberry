# Use for long read assembly programs

## Method was used for Fusarium oxysporum fsp frgariae
#### 1 # , then you have a header, ## subheading , ### small heading - use "tab" key to identify script info

####################
# Miniasm assembly
####################

## Step 1
#Open a screen - scrren -r
Login to node srun --partition long --time 0-18:00:00 --mem-per-cpu 20G --cpus-per-task 24 --pty bash

#Concatenate sequence reads first
Use command below if you are working in the same directory as the raw sequence reads

#Sequnece reads In /archives/2020_niabemr_nanopore/Fo_fragariae_14_003/Fof14-003-Scha48-1/20201027_1719_X1_FAL69458_56dbaa42/fastq_pass$
/archives/2020_niabemr_nanopore/Fo_fragariae_14_003/Fo_fragariae_14_003_R2/20201028_1622_X1_FAL69458_e6fdc4e2/fastq_pass

#Did 2 runs therefore  must join both

    cat *fastq | gzip -cf > Fof14R1.fastq.gz
    cat *fastq | gzip -cf > Fof14R2.fastq.gz

#Copy concatenated file to raw_dna directory

    cp Fof14R1.fastq.gz /projects/fusarium_ex_strawberry/NextGenSeq/raw_dna
    cp Fof14R2.fastq.gz /projects/fusarium_ex_strawberry/NextGenSeq/raw_dna

#Merge the contatenated files in directory they were copied into

    cat Fof14R1.fastq.gz Fof14R2.fastq.gz > Fof14RT.fastq.gz

## Step 2
#Run Porechop before assembly

    /scratch/software/Porechop-0.2.3/porechop-runner.py -i Fof14RT.fastq.gz -o Fof14RT_trim.fastq.gz --threads 16 > Fof14RT_trim_log.txt


#For ONT long read sequences use miniasm to assemble genome
## Step 3 Run in olc_assemblers conda shell
#If script doesn't work, see below

    for TrimReads in $(ls raw_dna/Fof14RT.fastq.gz); do
      Organism=F.oxysporum_fsp_fragariae
      Strain=DSA14_003
      Prefix="$Strain"_miniasm
      OutDir=assembly/miniasm/$Organism/$Strain
      mkdir -p $OutDir
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
      sbatch $ProgDir/miniasm.sh $TrimReads $Prefix $OutDir
  done

## Step 3a
#Need to rename all reads

Run in conda env that has minimap2 and bbmap

    rename.sh qin=33 in=Fof14RT_trim.fastq.gz out=Fof14RT_renamed.fasta prefix=Fof14

## Step 4a
#Run minimap2
Run in olc_assemblers conda shell
#Need to run read against itself. IT IS NEEDED. it is going to do self mapping
Fast all-against-all overlap of raw reads

    minimap2 -x ava-ont -t8 Fof14RT_renamed.fasta Fof14RT_renamed.fasta | gzip -1 > Fof14_fastq_allfiles.paf.gz


#Run in a screen and a node in a conda env with miniasm installed
## Step 5
#Concatenate pieces of read sequences to generate the final sequences
Can run like this instead: miniasm -f trimmed_renamed.fasta FolR1_fastq_allfiles.paf.gz > reads.gfa

    miniasm -f raw_dna/Fof14RT_renamed.fasta raw_dna/Fof14_fastq_allfiles.paf.gz > assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/reads.gfa

## Step 6
#Convert gfa file to fasta file

    awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > Fof14_miniasm.fa

####################
# Flye assembly
####################

    for TrimReads in $(ls raw_dna/Fof14RT.fastq.gz); do
      Organism=F.oxysporum_fsp_fragariae
      Strain=DSA14_003
       Prefix="$Strain"_flye;     
       TypeSeq=nanoraw;
       OutDir=assembly/flye/$Organism/$Strain/;
       mkdir -p $OutDir;
       Size=60m; # size= Expected genome size
       ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/;
       sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size $TypeSeq;
     done

######################
# SMARTDenovo assembly
######################

     for TrimReads in $(ls raw_dna/Fof14RT.fastq.gz); do
       Organism=F.oxysporum_fsp_fragariae
       Strain=DSA14_003
       Prefix="$Strain"_smartdenovo
       OutDir=assembly/SMARTdenovo/$Organism/$Strain
       mkdir -p $OutDir
       ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
       sbatch $ProgDir/SMARTdenovo.sh $TrimReads $Prefix $OutDir
     done

#output = *_smartdenovo.dmo.lay.utg

#####################
# QC steps
#####################


## Quast QC assembly check
#Run in conda env with python 2.7 (betaenv)
#Run on each assembly

    ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
      for Assembly in $(ls assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/assembly.fasta); do
        Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
        Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
        OutDir=assembly/flye/$Organism/$Strain/ncbi_edits
        sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
      done
#Updated entire script using https://github.com/harrisonlab/bioinformatics_tools/blob/master/Gene_prediction/README.md
#Look into BuscoDB direc - directory exists
#Run in conda env - BUSCOenv
#Ran on genome(softmasked) and gene models (final_genes_appended_renamed.gene.fasta)

    for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/DSA14_003_smartdenovo.dmo.lay.utg); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
      Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
      echo "$Organism - $Strain"
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done

#####################
# Racon
#####################

#Racon generates 10 iterations which have polished the genome
Need to do for Miniasm*, Flye* and SMARTdenovo* output files
#Run in condaenv with racon installed (olc_assemblers) - *=complete

    for Assembly in $(ls assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/Fof14_miniasm.fa); do
        ReadsFq=$(ls raw_dna/Fof14RT_trim.fastq.gz)
        Iterations=10
        OutDir=$(dirname $Assembly)"/racon_$Iterations"
        ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
        sbatch $ProgDir/racon.sh $Assembly $ReadsFq $Iterations $OutDir
      done

#Quality check each iteration with Quast and BUSCO
Run in correct conda env (betaenv)
Run for each round of each assembly:
*assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/assembly_racon_round_1.fasta
*assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/DSA14_003_smartdenovo_racon_round_1.fasta
*assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/Fof14_miniasm_racon_round_1.fasta

      ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
        for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/DSA14_003_smartdenovo_racon_round_1.fasta); do
          Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
          Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
          OutDir=assembly/SMARTdenovo/$Organism/$Strain/racon_10/ncbi_edits/round_1
          sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
        done

Run in correct conda env (BUSCOenv)
Run for each round of each assembly:
*assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/assembly_racon_round_1.fasta
*assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/DSA14_003_smartdenovo_racon_round_1.fasta
*assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/Fof14_miniasm_racon_round_1.fasta

    for Assembly in $(ls assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/Fof14_miniasm_racon_round_1.fasta); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
      Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
      echo "$Organism - $Strain"
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/busco_sordariomycetes_obd10/round_1
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done

Select iteration from each assembly with the best BUSCO scores for the next step


#####################
# Medaka
#####################

## Rename contigs for genome
If split or remove contigs is needed, provide FCSreport file by NCBI.

*assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/DSA14_003_smartdenovo_racon_round_4.fasta
*assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/Fof14_miniasm_racon_round_7.fasta
*assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/assembly_racon_round_3.fasta

    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
        touch tmp.txt
        for Assembly in $(ls assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/assembly_racon_round_3.fasta); do
            OutDir=$(dirname $Assembly)
            $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/assembly_racon_round_3_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
        done
        rm tmp.txt


## Run in conda env with medaka installed (medaka)
A tool to create a consensus sequence from nanopore sequencing data.
This task is performed using neural networks applied from a pileup of individual sequencing reads against a draft assembly.
It outperforms graph-based methods operating on basecalled data, and can be competitive with state-of-the-art signal-based methods, whilst being much faster.

*assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/DSA14_003_smartdenovo_racon_round_4_renamed.fasta
*assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/Fof14_miniasm_racon_round_7_renamed.fasta
*assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/assembly_racon_round_3_renamed.fasta

    for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/racon_10/DSA14_003_smartdenovo_racon_round_4_renamed.fasta); do
      ReadsFq=$(ls raw_dna/Fof14RT_trim.fastq.gz)
      OutDir=assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/medaka
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
      sbatch $ProgDir/medaka.sh $Assembly $ReadsFq $OutDir
    done

Run BUSCO and quast on medaka and pick best 1 out of 3 for pilon polishing

    for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/medaka/DSA14_003_smartdenovo_racon_round_4_renamed.fasta); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
      Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
      echo "$Organism - $Strain"
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/medaka/busco_sordariomycetes_obd10
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done

#####################
# Pilon
#####################

Aligning illumina reads against pilon data to polish.
Alternate prog directory /home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
Run in conda env.

INSTALL BOWTIE2 -

    conda install -c bioconda bowtie2

Make sure script is executable

   chmod u+x ./sub_pilon.sh

The best medaka product was from the miniasm assembly - Fof14_miniasm_racon_round_7_renamed.fasta

    for Assembly in $(ls assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/medaka/Fof14_miniasm_racon_round_7_renamed.fasta); do
      Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
      Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
      IlluminaDir=$(ls -d ../raw_dna/paired/F.oxysporum_fsp_fragariae/DSA14_003)
      echo $Strain
      echo $Organism
      TrimF1_Read=$(ls $IlluminaDir/F/*.fastq.gz | head -n2 | tail -n1);
      TrimR1_Read=$(ls $IlluminaDir/R/*.fastq.gz | head -n2 | tail -n1);
      echo $TrimF1_Read
      echo $TrimR1_Read
      OutDir=assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon
      Iterations=10
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
      sbatch $ProgDir/pilon_2_lib.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
    done

Script failed with pilon_1_lib.sh
Edited script inputs in pilon_2_lib.sh script - failed in conda env
failed outside conda env - died after 15s
Installed bowtie2 in olc_assemblers env - worked

Run quast and BUSCO analysis on each iteration
QUAST

    ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
      for Assembly in $(ls assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_*.fasta); do
        Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
        Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
        OutDir=assembly/miniasm/$Organism/$Strain/pilon/ncbi_edits/round_*
        sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
      done

BUSCO

    for Assembly in $(ls assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_*.fasta); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
      Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
      echo "$Organism - $Strain"
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon/busco_sordariomycetes_obd10/round_*
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done


#####################
# Repeat Masking
#####################

Repeat identification and masking is conducted before gene prediction and annotation steps.
The term 'masking' means transforming every nucleotide identified as a repeat to an 'N', 'X' or to a lower case a, t, g or c.

## Repeat mask
Run in conda env (Repenv or seq_tools)

    for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta | grep -v '_2' | grep -v '11055'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism"
    echo "$Strain"
    OutDir=repeat_masked/$Organism/$Strain/ncbi_edits_repmask
    ProgDir=/home/akinya/git_repos/tools/seq_tools/repeat_masking
    sbatch $ProgDir/rep_modeling.sh $Assembly $OutDir
    sbatch $ProgDir/transposonPSI.sh $Assembly $OutDir
    done

## Soft mask
Soft masking means transforming every nucleotide identified as a repeat to a lower case a, t, g or c to be included in later gene prediction stages.

Gives number of masked N's in sequence  - Take physical and digital note of the results.

    for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa); do
      OutDir=$(dirname $File)
      TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
      OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
      echo "$OutFile"
      bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
      echo "Number of masked bases:"
      cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
    done

## Hard Mask
Hard masking  means transforming every nucleotide identified as a repeat to an 'N' or 'X'.

    for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa); do
      OutDir=$(dirname $File)
      TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
      OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
      echo "$OutFile"
      bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
    done

Run BUSCO and Quast qc checks on the softmasked, unmasked and hardmasked assemblies
