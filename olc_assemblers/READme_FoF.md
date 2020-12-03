# Use for long read assembly programs

## Method was used for Fusarium oxysporum fsp frgariae
#### 1 # , then you have a header, ## subheading , ### small heading - use "tab" key to identify script info


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


# QC steps
#####################


## Quast QC assembly check
Run in conda env with python 2.7 (betaenv)

    conda create -n betaenv python=2.7
    # conda create -n betaenv python=2.7 quast=0.1.0 - to install specific version of quast
    conda install -c bioconda quast


Run on each assembly

    ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
      for Assembly in $(ls assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/assembly.fasta); do
        Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
        Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
        OutDir=assembly/flye/$Organism/$Strain/ncbi_edits
        sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
      done
#Updated entire script using https://github.com/harrisonlab/bioinformatics_tools/blob/master/Gene_prediction/README.md
#Look into BuscoDB direc - directory exists
Run in conda env - BUSCOenv

    conda create -n BUSCOenv
    conda install -c bioconda busco

Ran on genome(softmasked) and gene models (final_genes_appended_renamed.gene.fasta)

    for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/DSA14_003_smartdenovo.dmo.lay.utg); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
      Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
      echo "$Organism - $Strain"
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done


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

# Pilon
#####################

Aligning illumina reads against pilon data to polish.
Alternate prog directory /home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
Run in conda env (olc_assemblers).

INSTALL BOWTIE2 -

    conda install -c bioconda bowtie2

Make sure script is executable

   chmod u+x ./sub_pilon.sh

The best medaka product was from the miniasm assembly - Fof14_miniasm_racon_round_7_renamed.fasta
DUE to possible chimeric miniasm assembly, ran pilon on SDEN and flye assemblies


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
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/NGS_assembly
      sbatch $ProgDir/pilon_2_lib.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
    done

Script failed with pilon_1_lib.sh
Edited script inputs in pilon_2_lib.sh script - failed in conda env
failed outside conda env - died after 15s
Installed bowtie2 in olc_assemblers env - worked

Run quast and BUSCO analysis on each iteration
QUAST

    ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
      for Assembly in $(ls assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_*.fasta); do
        Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
        Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
        OutDir=assembly/flye/$Organism/$Strain/pilon/ncbi_edits/round_*
        sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
      done

BUSCO

    for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_*.fasta); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
      Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
      echo "$Organism - $Strain"
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/quality_check/
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/pilon/busco_sordariomycetes_obd10/round_*
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done


# Repeat Masking
#####################

Repeat identification and masking is conducted before gene prediction and annotation steps.
The term 'masking' means transforming every nucleotide identified as a repeat to an 'N', 'X' or to a lower case a, t, g or c.

## Repeat mask
Ensure packages are installed in envs

    conda create -n RMask

    conda install -c bioconda repeatmodeler # repeatmodeler also installs packages below
    #conda install -c bioconda repeatmasker
    #conda install rmblast

Need to manually configure Repeatmasker

    cd /home/USER_ID/miniconda3/envs/general_tools/share/RepeatMasker/ # USER_ID is your user name i.e. akinya

    ./confiure # runs the configuration step

    # Set execution path of tfr, e.g. /home/USER_ID/miniconda3/envs/USER_ENV/bin/trf
    # Add search engine. Option 2 - RMBlast will be used
    # Set path where rmblastn and makeblastdb are found, e.g. /home/USER_ID/miniconda3/envs/USER_ENV/bin
    # 5. Done to exit

### Rename before you run rep mask

ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
    touch tmp.txt
    for Assembly in $(ls assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10.fasta); do
        OutDir=$(dirname $Assembly)
        $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/pilon_10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
    done
    rm tmp.txt

Have 2 paths to choose from to run scripts

### RepeatMask & TPSI path 1

    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking
    BestAssembly=assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta
    OutDir=repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/SMARTdenovo/ncbi_edits_repmask
    sbatch $ProgDir/rep_modeling.sh $BestAssembly $OutDir
    sbatch $ProgDir/transposonPSI.sh $BestAssembly $OutDir

### Rep mask (path 2)

Run in conda env (Repenv) - input for |illumina assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta | grep -v '_2' | grep -v '11055'|

    for Assembly in $(ls assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism"
    echo "$Strain"
    OutDir=repeat_masked/$Organism/$Strain/ncbi_edits_repmask
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking #/home/akinya/git_repos/tools/seq_tools/repeat_masking
    sbatch $ProgDir/rep_modelingBeta.sh $Assembly $OutDir
    done

### TransposonPSI (path 2)
Run in RMask env

    conda install -c bioconda transposonpsi

    for Assembly in $(ls assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism"
    echo "$Strain"
    OutDir=repeat_masked/$Organism/$Strain/ncbi_edits_repmask
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking #/home/akinya/git_repos/tools/seq_tools/repeat_masking
    sbatch $ProgDir/transposonPSI.sh $Assembly $OutDir
    done

## Soft mask
Soft masking means transforming every nucleotide identified as a repeat to a lower case a, t, g or c to be included in later gene prediction stages.

Gives number of masked N's in sequence  - Take physical and digital note of the results.

    for File in $(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_contigs_softmasked.fa); do
      OutDir=$(dirname $File)
      TPSI=$(ls $OutDir/DSA14_003_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
      OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
      echo "$OutFile"
      bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
      echo "Number of masked bases:"
      cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
    done

    # Number of masked bases:
    # miniasm - 6854263     SDEN- 6908635
    # flye - 7069053

## Hard Mask
Hard masking  means transforming every nucleotide identified as a repeat to an 'N' or 'X'.

    for File in $(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/ncbi_edits_repmask/DSA14_003_contigs_hardmasked.fa); do
      OutDir=$(dirname $File)
      TPSI=$(ls $OutDir/DSA14_003_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
      OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
      echo "$OutFile"
      bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
    done

Run BUSCO and Quast qc checks on the softmasked, unmasked and hardmasked assemblies

# Orthology hunt
#####################

Use extracted effectors from a gene.fasta file using names in txt file
Create text file with known gene names of effectors/mimps you want to remove
i.e  input gene names in Fof14_genes.txt
Command uses gene name to copy across fasta seq to Fof genes

Why are you doing this?
To compare candidate effectors in Fo cepae (which has RNA seq data) against the predicted genes in Fof using cepae RNA seq data

Extract genes from reference lycopersici & cepae genomes NOT FRAGARIAE

  faidx -d '|' final_genes_appended_renamed.cdna.fasta $(tr '\n' ' ' < Fof14_genes.txt ) > Fof14_genes.fasta - example command to test

  #faidx -d '|' final_genes_combined.cdna.fasta $(tr '\n' ' ' < FoC_cand_mimps.txt ) > Eff_mimp_genes.fasta
  #faidx -d '|' final_genes_combined.cdna.fasta $(tr '\n' ' ' < FoC_Six.txt ) > six_ortho_genes.fasta

-Now that you have the cand effector genes, contrast against Fo_fragariae_14_003 long read seq genome
-Run in conda env with perly (Repenv)
-For $Assembly Use files with nucleotides

  for Assembly in $(ls assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    Query=../F.oxysporum_fsp_cepae/Fus2_canu_new/final/Eff_mimp_genes.fasta # six_ortho_genes.fasta
    OutDir=assembly/flye/$Organism/$Strain/Orthology/FoFrvsFoCep_mimps
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
    sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly $OutDir
  done

Second query../F.oxysporum_fsp_cepae/Fus2_canu_new/final/Eff_mimp_genes.fasta
Use cut -f1 or * DSA14_003_eff_ortho_genes.fasta_homologs.csv to excise and view column data
Compare against Andy's known SIX genes

  for Assembly in $(ls assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    Query=../../oldhome/groups/harrisonlab/project_files/fusarium/analysis/blast_homology/six_genes/six-appended_parsed.fa
    OutDir=assembly/miniasm/$Organism/$Strain/Orthology
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
    sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly $OutDir
  done

- miniasm assembly was potentially chimeric, therefore did the same for flye and Sden assemblies


# Synteny Check
#####################

# D-genies

Go to http://dgenies.toulouse.inra.fr/ to compare genomes for synteny against FoLy4287 and FoFrvsFoCep.
Using the plots, see which assembly will be best to use for gene_prediction using Fo_cepae data.

  assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta # best assembly


# Gene prediction
#####################

## STAR

Run in Repenv - condaenv
Fus2_CzapekDox, Fus2_GlucosePeptone, Fus2_PDA and Fus2_PDB are RNAseq data of infected onions
Only data samples that will map genes of F.oxy accurately
Need to concatenate data after STAR analysis
#--genomeSAindexNbases is unique to each genome and is 11 for FoFR

  for Assembly in $(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_contigs_unmasked.fa);  do
      Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
      Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
      echo "$Organism - $Strain"
      FileF=../../oldhome/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/F/*_trim.fq.gz
      FileR=../../oldhome/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/R/*_trim.fq.gz
      echo $FileF
      echo $FileR
      Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
      echo "$Timepoint"
      Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
      OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/Genome_alignment
      sbatch $ProgDir/STAR_1.sh $Assembly $FileF $FileR $OutDir 11
    done

Need to concatenate in this step to link RNAseq data into one series
View "star_aligmentLog.final.out" to see uniquely mapped reads %

  Strain=DSA14_003
    Organism=F.oxysporum_fsp_fragariae
    mkdir -p alignment/star/$Organism/$Strain/concatenated
    samtools merge -f alignment/star/$Organism/$Strain/concatenated/concatenated.bam \
    alignment/star/$Organism/$Strain/Fus2_CzapekDox/6_S2_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam \
    alignment/star/$Organism/$Strain/Fus2_GlucosePeptone/7_S3_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam \
    alignment/star/$Organism/$Strain/Fus2_PDA/9_S4_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam \
    alignment/star/$Organism/$Strain/Fus2_PDB/4_S1_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam


## Braker

Run in conda env (Repenv)
AcceptedHits=alignment/concatenated.bam
Alternate strain for softmasked
Intial run required installation of Hash::Merge and Logger::Simple using cpan

    #Original prog  dir /home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction

    for Assembly in $(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
        Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
        echo "$Organism - $Strain"
        OutDir=gene_pred/braker/$Organism/$Strain/flye
        AcceptedHits=alignment/star/F.oxysporum_fsp_fragariae/DSA14_003/concatenated/concatenated.bam
        GeneModelName="$Organism"_"$Strain"_braker_flye
        ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/ProgScripts
        sbatch $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
      done

Got this error:
    failed to execute: perl /home/gomeza/prog/genemark/gmes_linux_64/gmes_petap.pl --sequence=/tmp/akinya_614617/braker/F.oxysporum_fsp_fragariae_DSA14_003_braker_flye/genome.fa --ET=/tmp/akinya_614617/braker/F.oxysporum_fsp_fragariae_DSA14_003_braker_flye/hintsfile.gff --cores=1 --fungus --soft 1000 1>/tmp/akinya_614617/braker/F.oxysporum_fsp_fragariae_DSA14_003_braker_flye/GeneMark-ET.stdout 2>/tmp/akinya_614617/braker/F.oxysporum_fsp_fragariae_DSA14_003_braker_flye/errors/GeneMark-ET.stderr

## StringTie

String tie - to be edited
Run in conda env with Python 2.7 (betaenv)
Codingquarry is another tool for gene prediction that it is able to predict additional genes in fungi
Merge with Braker to give final gene model set
/home/akinya/git_repos/assembly_fusarium_ex/scripts
/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction

    for Assembly in $(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_contigs_unmasked.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
        Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev)
        echo "$Organism - $Strain"
        OutDir=gene_pred/stringtie/$Organism/$Strain/flye/concatenated_prelim
        mkdir -p $OutDir
        AcceptedHits=alignment/star/F.oxysporum_fsp_fragariae/DSA14_003/concatenated/concatenated.bam
        ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/ProgScripts
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

## Add transcripts together

Additional transcripts - to be edited
Run in perly env (Repenv)
Type full paths, do not use asterisks
RUN LINE BY LINE AS IT WILL NOT WORK
Do segments one at a time for peace of mind

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

Create fasta files from each gene feature in the CodingQuarry gff3
    $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary

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

## Gene renaming
Run line by line
Run in conda env (Repenv)
Remove duplicate and rename genes

    GffAppended=$(ls -d gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3)
    Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    FinalDir=gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final

    Remove duplicated genes
    GffFiltered=$FinalDir/filtered_duplicates.gff
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction
    $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered

    Rename genes
    GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
    LogFile=$FinalDir/final_genes_appended_renamed.log
    $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
    rm $GffFiltered

    Create renamed fasta files from each gene feature
    Assembly=$(ls repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    $ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_appended_renamed
    The proteins fasta file contains * instead of Xs for stop codons, these should be changed
    sed -i 's/\*/X/g' $FinalDir/final_genes_appended_renamed.pep.fasta

    view gene names
    cat $FinalDir/final_genes_appended_renamed.cdna.fasta | grep '>'
