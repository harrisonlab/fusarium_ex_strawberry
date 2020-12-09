# Use for long read assembly programs

# Method was used for Fusarium oxysporum fsp lactucae
#### 1 # , then you have a header, ## subheading , ### small heading - use "tab" key to identify script info

####################
# Miniasm assembly
####################

## Step 1
#Login to node srun --partition himem --time 0-06:00:00 --mem-per-cpu 40G --cpus-per-task 24 --pty bash
Concatenate sequence reads first (there is old seq data that was basecalled again include it)
#Use command below if you are working in the same directory as the raw sequence reads

    cat *fastq | gzip -cf > FAL69458.fastq.gz
    /archives/2020_niabemr_nanopore/F.oxyspporum_lactucae_Race1/20180426_AJ520_GA30000$ cat *.fastq.gz | gzip -cf > /projects/fusarium_EX_Lactucae/raw_dna/FoLR12018NBC.fastq.gz
    cat *.fastq.gz | gzip -cf > /projects/fusarium_EX_Lactucae/raw_dna/FoLR1-18nbc.fastq.gz
    cat FAL69458.fastq.gz FoLR1-18nbc.fastq.gz FoLR12018NBC.fastq.gz > FolR1_conc.fastq.gz

#Method above caused error in porechop (died 3 times without running possibly due to concatenating 3 big gz files???)

Tried method 2 - gave same error
    cat 20180426_AJ520_GA30000/*.fastq.gz basecalling-all/*.fastq.gz fastq_pass/*.fastq | gzip -cf > /projects/fusarium_EX_Lactucae/raw_dna/FoLr1cont.fastq.gz

copied gzip files to /projects/fusarium_EX_Lactucae/raw_dna/concatenated
Unzipped the gzip files.
    gzip -d *fastq.gz

Ran pore chop on the unzipped files in the conctenated/ directory
    /scratch/software/Porechop-0.2.3/porechop-runner.py -i concatenated/ -o FoL_CONC_trim.fastq.gz --threads 16 > FAL_CONC_trim_log.txt

## Step 2
#Run Porechop before assembly
Kept the errors to show what can go wrong with trying to run porechop on gzipped re-basecalled runs

    /scratch/software/Porechop-0.2.3/porechop-runner.py -i FAL69458.fastq.gz -o FAL_trim.fastq.gz --threads 16 > FAL_trim_log.txt
    #/scratch/software/Porechop-0.2.3/porechop-runner.py -i FolR1_conc.fastq.gz -o FolR1_conc_trim.fastq.gz --threads 16 > FolR1_conc_trim_log.txt - method failed (multiple times)
    #/scratch/software/Porechop-0.2.3/porechop-runner.py -i FoLr1cont.fastq.gz -o FoLr1cont_trim.fastq.gz --threads 16 > FoLr1cont_trim_log.txt - same error

    Traceback (most recent call last):
      File "/scratch/software/Porechop-0.2.3/porechop-runner.py", line 9, in <module>
        main()
      File "/scratch/software/Porechop-0.2.3/porechop/porechop.py", line 34, in main
        reads, check_reads, read_type = load_reads(args.input, args.verbosity, args.print_dest,
      File "/scratch/software/Porechop-0.2.3/porechop/porechop.py", line 230, in load_reads
        reads, read_type = load_fasta_or_fastq(input_file_or_directory)
      File "/scratch/software/Porechop-0.2.3/porechop/misc.py", line 114, in load_fasta_or_fastq
        file_type = get_sequence_file_type(filename)
      File "/scratch/software/Porechop-0.2.3/porechop/misc.py", line 106, in get_sequence_file_type
        raise ValueError('File is neither FASTA or FASTQ')
    ValueError: File is neither FASTA or FASTQ

Attempted to run porechop on conctenated gzip folders in directory
error below
    /scratch/software/Porechop-0.2.3/porechop-runner.py -i concatenated/ -o FoL_CONC_trim.fastq.gz --threads 16 > FAL_CONC_trim_log.txt
    Traceback (most recent call last):
    File "/scratch/software/Porechop-0.2.3/porechop-runner.py", line 9, in <module>
     main()
    File "/scratch/software/Porechop-0.2.3/porechop/porechop.py", line 35, in main
     args.check_reads)
    File "/scratch/software/Porechop-0.2.3/porechop/porechop.py", line 255, in load_reads
     file_reads, _ = load_fasta_or_fastq(fastq_file)
    File "/scratch/software/Porechop-0.2.3/porechop/misc.py", line 114, in load_fasta_or_fastq
     file_type = get_sequence_file_type(filename)
    File "/scratch/software/Porechop-0.2.3/porechop/misc.py", line 106, in get_sequence_file_type
     raise ValueError('File is neither FASTA or FASTQ')
    ValueError: File is neither FASTA or FASTQ    


## Step 3
#Need to rename all reads
Run in conda env that has minimap2 and bbmap (olc_assemblers)

    rename.sh qin=33 in=FAL_trim.fastq.gz out=trimmed_renamed.fasta prefix=FolR1

For all concatenated reads run:

    rename.sh qin=33 in=FoL_CONC_trim.fastq.gz out=FoL_conc_renamed.fasta prefix=FolR1

#If script doesn't work, see below

        for TrimReads in $(ls FAL_trim.fastq.gz); do #sub in FoL_CONC_trim.fastq.gz
            Organism=F.oxysporum_fsp_lactucae
            Strain=race_1
            Prefix="$Strain"_miniasm
            OutDir=assembly/miniasm/$Organism/$Strain #Assembly2
            mkdir -p $OutDir
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
            sbatch $ProgDir/miniasm.sh $TrimReads $Prefix $OutDir
          done

## Step 4
#Run minimap2
Need to run read against itself. IT IS NEEDED. it is going to do self mapping
#Fast all-against-all overlap of raw reads

    minimap2 -x ava-ont -t8 trimmed_renamed.fasta trimmed_renamed.fasta | gzip -1 > FolR1_fastq_allfiles.paf.gz
    or
    minimap2 -x ava-ont -t8 FoL_conc_renamed.fasta FoL_conc_renamed.fasta | gzip -1 > FolR1C_fastq_all.paf.gz

#For ONT long read sequences use miniasm to assemble genome
Run in a screen and a node in a conda env with miniasm installed
#Run in olc_assemblers conda shell

## Step 5
#Concatenate pieces of read sequences to generate the final sequences
Can run like this instead: miniasm -f trimmed_renamed.fasta FolR1_fastq_allfiles.paf.gz > reads.gfa

    miniasm -f "$Prefix"_rename.fasta $Prefix.paf.gz > reads.gfa
    #or
    miniasm -f FoL_conc_renamed.fasta FolR1C_fastq_all.paf.gz > Assembly2/miniasm/F.oxysporum_fsp_lactucae/race_1/reads.gfa

## Step 6
#Convert gfa file to fasta file

    awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > $Prefix.fa


#######################
# Flye assembly
#######################

#Run in olc_assemblers env
Log into the long node
#flye assembly method
size= Expected genome size

    for TrimReads in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz); do # for concatenated runs raw_dna/FoL_CONC_trim.fastq.gz
           Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev) ;
           Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev) ;
           Prefix="$Strain"_flye;     
           TypeSeq=nanoraw;
           OutDir=assembly/flye/$Organism/$Strain/flye_raw;
           mkdir -p $OutDir;
           Size=60m; # size= Expected genome size
           ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/;
           sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size $TypeSeq;
         done

         or

         for TrimReads in $(ls raw_dna/FoL_CONC_trim.fastq.gz); do # for concatenated runs
                Organism=F.oxysporum_fsp_lactucae;
                Strain=race_1;
                Prefix="$Strain"_flye;     
                TypeSeq=nanoraw;
                OutDir=Assembly2/flye/$Organism/$Strain;
                mkdir -p $OutDir;
                Size=60m;
                ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers;
                sbatch $ProgDir/flye.sh $TrimReads $Prefix $OutDir $Size $TypeSeq;
              done

########################
# SMARTDenovo assembly
########################

## SMartDenovo script
#Run in olc_assemblers env
Use raw inputs unlike miniasm
#Use porechop trimmed output
Could run like this:
#~/miniconda3/envs/olc_assemblers/bin/smartdenovo.pl -p race_1_smartdenovo -t 14 -c 1 FolR1_fastq_allfiles.paf.gz > race_1_smartdenovo.mak
#make -f prefix.mak

    for TrimReads in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz); do
      Organism=F.oxysporum_fsp_lactucae
      Strain=race_1
      Prefix="$Strain"_smartdenovo
      OutDir=assembly/SMARTdenovo/$Organism/$Strain
      mkdir -p $OutDir
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
      sbatch $ProgDir/SMARTdenovo.sh $TrimReads $Prefix $OutDir
    done
     or
     for TrimReads in $(ls raw_dna/FoL_CONC_trim.fastq.gz); do
       Organism=F.oxysporum_fsp_lactucae
       Strain=race_1
       Prefix="$Strain"_smartdenovo
       OutDir=Assembly2/SMARTdenovo/$Organism/$Strain
       mkdir -p $OutDir
       ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
       sbatch $ProgDir/SMARTdenovo.sh $TrimReads $Prefix $OutDir
     done

#output = race_1_smartdenovo.dmo.lay.utg

#####################
# QC steps
#####################


## Quast QC assembly check
#Run in conda env with python 2.7 (betaenv)
Run on each assembly

    ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
      for Assembly in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/assembly.fasta); do
        Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
        Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
        OutDir=assembly/flye/$Organism/$Strain/ncbi_edits
        sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
      done
      or
      ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
      for Assembly in $(ls Assembly2/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/race_1_smartdenovo.dmo.lay.utg); do
        Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
        Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
        OutDir=Assembly2/miniasm/$Organism/$Strain/ncbi_edits
        sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
      done

#Updated entire script using https://github.com/harrisonlab/bioinformatics_tools/blob/master/Gene_prediction/README.md
#Look into BuscoDB direc - directory exists
#Run in conda env - BUSCOenv

    for Assembly in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/assembly.fasta); do # Assembly2/miniasm/F.oxysporum_fsp_lactucae/race_1/FoLR1_conc.fa
      Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
      Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
      echo "$Organism - $Strain"
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done

#####################
# Assembly Polishing
#####################

## Racon
#Racon generates 10 iterations which have polished the genome
Need to do for Miniasm*, Flye* and SMARTdenovo* output files
#Run in condaenv with racon installed (olc_assemblers) - *=complete

    for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/race_1_smartdenovo.dmo.lay.utg); do
        ReadsFq=$(ls raw_dna/FoL_CONC_trim.fastq.gz)
        Iterations=10
        OutDir=$(dirname $Assembly)"/racon_$Iterations"
        ProgDir=~/git_repos/assembly_fusarium_ex/ProgScripts
        sbatch $ProgDir/racon.sh $Assembly $ReadsFq $Iterations $OutDir
      done

      or

      for Assembly in $(ls Assembly2/flye/F.oxysporum_fsp_lactucae/race_1/assembly.fasta); do
          ReadsFq=$(ls raw_dna/FoL_CONC_trim.fastq.gz)
          Iterations=10
          OutDir=$(dirname $Assembly)"/racon_$Iterations"
          ProgDir=~/git_repos/assembly_fusarium_ex/ProgScripts
          sbatch $ProgDir/racon.sh $Assembly $ReadsFq $Iterations $OutDir
        done

#Quality check each iteration with Quast and BUSCO
DO for each iteration

    ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
      for Assembly in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/flye_raw/racon_10/assembly_racon_round_1.fasta); do
        Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
        Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)  
        OutDir=$(dirname $Assembly)/ncbi_edits/round_1
        sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
      done

      or

      ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
        for Assembly in $(ls Assembly2/flye/F.oxysporum_fsp_lactucae/race_1/racon_10/assembly_racon_round_1.fasta); do
          Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
          Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
          OutDir=$(dirname $Assembly)/ncbi_edits/round_1
          sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
        done

#Ended up here for some reason assembly/flye/F.oxysporum_fsp_lactucae/race_1/ncbi_edits/round_*

    for Assembly in $(ls Assembly2/flye/F.oxysporum_fsp_lactucae/race_1/racon_10/assembly_racon_round_1.fasta); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
      Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
      echo "$Organism - $Strain"
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/quality_check
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10/round_1
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done

#Contigs must be renamed before medaka can be run
Rename contigs for genome
If split or remove contigs is needed, provide FCSreport file by NCBI.

    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
        touch tmp.txt
        for Assembly in $(ls Assembly2/flye/F.oxysporum_fsp_lactucae/race_1/racon_10/assembly_racon_round_4.fasta); do
            OutDir=$(dirname $Assembly)
            $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/assembly_racon_round_4.fasta_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
        done
        rm tmp.txt

## Medaka
#Using the best iteration shown by the Quast and BUSCO scores, Medaka was conducted to further polish the genome

#Run in medaka env
#A tool to create a consensus sequence from nanopore sequencing data.
#This task is performed using neural networks applied from a pileup of individual sequencing reads against a draft assembly.
#It outperforms graph-based methods operating on basecalled data, and can be competitive with state-of-the-art signal-based methods, whilst being much faster.

    for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/racon_10/*_renamed.fasta); do
      ReadsFq=$(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz)
      OutDir=$(dirname $Assembly)/medaka
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
      sbatch $ProgDir/medaka.sh $Assembly $ReadsFq $OutDir
    done
    or
    for Assembly in $(ls Assembly2/flye/F.oxysporum_fsp_lactucae/race_1/racon_10/assembly_racon_round_4_renamed.fasta); do
      ReadsFq=raw_dna/FoL_CONC_trim.fastq.gz
      OutDir=Assembly2/flye/F.oxysporum_fsp_lactucae/race_1/medaka
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
      sbatch $ProgDir/medaka.sh $Assembly $ReadsFq $OutDir
    done

Run QC checks with QUAST and BUSCO again to see any changes

# Pilon
#####################

Aligning illumina reads against pilon data to polish.
Alternate prog directory /home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
Run in conda env (olc_assemblers).

INSTALL BOWTIE2 -

    conda install -c bioconda bowtie2

Make sure script is executable

   chmod u+x ./sub_pilon.sh

Raw DNA direc - /projects/oldhome/groups/harrisonlab/project_files/fusarium/raw_dna/paired/F.oxysporum_fsp_lactucae/AJ520/F/AJ520_S2_L001_R1_001.fastq.gz0
Assembly - Assembly2/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/medaka/race_1_smartdenovo_racon_round_2_renamed.fasta
Assembly - Assembly2/miniasm/F.oxysporum_fsp_lactucae/race_1/medaka/FoLR1_conc_racon_round_6_renamed.fasta
Assembly - Assembly2/flye/F.oxysporum_fsp_lactucae/race_1/medaka/assembly_racon_round_4.fasta_renamed.fasta

    for Assembly in $(ls Assembly2/flye/F.oxysporum_fsp_lactucae/race_1/medaka/assembly_racon_round_4_renamed.fasta); do
      Organism=F.oxysporum_fsp_lactucae
      Strain=race_1
      IlluminaDir=$(ls -d ../oldhome/groups/harrisonlab/project_files/fusarium/raw_dna/paired/F.oxysporum_fsp_lactucae/AJ520)
      echo $Strain
      echo $Organism
      TrimF1_Read=$(ls $IlluminaDir/F/AJ520_S2_L001_R1_001.fastq.gz | head -n2 | tail -n1);
      TrimR1_Read=$(ls $IlluminaDir/R/AJ520_S2_L001_R2_001.fastq.gz | head -n2 | tail -n1);
      echo $TrimF1_Read
      echo $TrimR1_Read
      OutDir=Assembly2/flye/F.oxysporum_fsp_lactucae/race_1/pilon
      Iterations=10
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/NGS_assembly
      sbatch $ProgDir/pilon_lac_mini_lib.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
    done

Run BUSCO & QUAST
    for Assembly in $(ls Assembly2/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/pilon/pilon_*.fasta); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
      Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
      echo "$Organism - $Strain"
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/quality_check
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10/round_*
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
    for Assembly in $(ls Assembly2/miniasm/F.oxysporum_fsp_lactucae/race_1/pilon/pilon_10.fasta); do
        OutDir=$(dirname $Assembly)
        $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/pilon_10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
    done
    rm tmp.txt

Have 2 paths to choose from to run scripts

### RepeatMask & TPSI path 1

    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking
    BestAssembly=Assembly2/miniasm/F.oxysporum_fsp_lactucae/race_1/pilon/pilon_10_renamed.fasta
    OutDir=repeat_masked/F.oxysporum_fsp_lactucae/race_1/miniasm/ncbi_edits_repmask
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

    for File in $(ls repeat_masked/F.oxysporum_fsp_lactucae/race_1/SMARTdenovo/ncbi_edits_repmask/race_1_contigs_softmasked.fa); do
      OutDir=$(dirname $File)
      TPSI=$(ls $OutDir/race_1_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
      OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
      echo "$OutFile"
      bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
      echo "Number of masked bases:"
      cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
    done

    # Number of masked bases:
    # miniasm - 11417581     SDEN- 1358959
    # flye -

## Hard Mask
Hard masking  means transforming every nucleotide identified as a repeat to an 'N' or 'X'.

    for File in $(ls repeat_masked/F.oxysporum_fsp_lactucae/race_1/SMARTdenovo/ncbi_edits_repmask/race_1_contigs_hardmasked.fa); do
      OutDir=$(dirname $File)
      TPSI=$(ls $OutDir/race_1_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
      OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
      echo "$OutFile"
      bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
    done

Run BUSCO and Quast qc checks on the softmasked, unmasked and hardmasked assemblies

# Orthology hunt
#####################

Use extracted effectors from a gene.fasta file using names in txt file
Create text file with known gene names of effectors/mimps you want to remove

Why are you doing this?
To compare candidate effectors in Fo cepae (which has RNA seq data) against the predicted genes in Fof using cepae RNA seq data

Extract genes from reference lycopersici & cepae genomes
Command uses gene name to copy across fasta seq to Fo genes
  #faidx -d '|' final_genes_combined.cdna.fasta $(tr '\n' ' ' < FoC_cand_mimps.txt ) > Eff_mimp_genes.fasta
  #faidx -d '|' final_genes_combined.cdna.fasta $(tr '\n' ' ' < FoC_Six.txt ) > six_ortho_genes.fasta

-Now that you have the cand effector genes, contrast against long read seq genome
-Run in conda env with perly (Repenv)
-For $Assembly Use files with nucleotides

  for Assembly in $(ls Assembly2/miniasm/F.oxysporum_fsp_lactucae/race_1/pilon/pilon_10_renamed.fasta); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    Query=../fusarium_ex_strawberry/F.oxysporum_fsp_cepae/Fus2_canu_new/final/Eff_mimp_genes.fasta # six_ortho_genes.fasta
    OutDir=Orthology/miniasm/blastn/$Organism/$Strain/FoFrvsFoCep_mimps
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
    sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly $OutDir
  done

Second query../F.oxysporum_fsp_cepae/Fus2_canu_new/final/Eff_mimp_genes.fasta
Use cut -f1 or * DSA14_003_eff_ortho_genes.fasta_homologs.csv to excise and view column data
Compare against Andy's known SIX genes

  for Assembly in $(ls Assembly2/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/pilon/pilon_10_renamed.fasta); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    Query=../oldhome/groups/harrisonlab/project_files/fusarium/analysis/blast_homology/six_genes/six-appended_parsed.fa
    OutDir=Orthology/SMARTdenovo/blastn/$Organism/$Strain/FoFrvsFoLy
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
    sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly $OutDir
  done


# Synteny Check
#####################

# D-genies

Go to http://dgenies.toulouse.inra.fr/ to compare genomes for synteny against FoLy4287 and FoFrvsFoCep.
Using the plots, see which assembly will be best to use for gene_prediction using Fo_cepae data.


# Gene prediction
#####################

## STAR

Run in Repenv - condaenv
Fus2_CzapekDox, Fus2_GlucosePeptone, Fus2_PDA and Fus2_PDB are RNAseq data of infected onions
Only data samples that will map genes of F.oxy accurately
Need to concatenate data after STAR analysis
#--genomeSAindexNbases is unique to each genome and is 11 for FoFR

  for Assembly in $(ls repeat_masked/F.oxysporum_fsp_lactucae/race_1/SMARTdenovo/ncbi_edits_repmask/race_1_contigs_unmasked.fa);  do
      Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
      Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
      echo "$Organism - $Strain"
      FileF=../oldhome/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/F/*_trim.fq.gz
      FileR=../oldhome/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/R/*_trim.fq.gz
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

  conda install -c thiesgehrmann genemark_es
  find miniconda3/envs/Repenv/ -name genemark_es # find location of program in installed env

Installation instructions for GeneMark* software

a. Copy the content of distribution to desired location.
b. Install the key: copy key "gm_key" into users home directory as:

  cp gm_key ~/.gm_key

Program is ready for execution.

add these paths to your "braker_fungi.sh" program script:
  --GENEMARK_PATH=/home/gomeza/prog/genemark/gmes_linux_64 \
  --BAMTOOLS_PATH=/home/gomeza/miniconda3/envs/gene_pred/bin \
or
  --GENEMARK_PATH=/home/akinya/miniconda3/envs/Repenv/opt/genemark_es/gmes_petap \
  --BAMTOOLS_PATH=/home/akinya/miniconda3/envs/Repenv/bin \

Then copy the .gm_key file like so:
  cp /home/gomeza/.gm_key ~/

    #Original prog  dir /home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction

    for Assembly in $(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
        Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
        echo "$Organism - $Strain"
        OutDir=gene_pred/braker/$Organism/$Strain/flye
        AcceptedHits=alignment/star/F.oxysporum_fsp_fragariae/DSA14_003/concatenated/concatenated.bam
        GeneModelName="$Organism"_"$Strain"_braker_flye_V2
        ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/ProgScripts
        sbatch $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
      done

Got this error:
    failed to execute: perl /home/gomeza/prog/genemark/gmes_linux_64/gmes_petap.pl --sequence=/tmp/akinya_614617/braker/F.oxysporum_fsp_fragariae_DSA14_003_braker_flye/genome.fa --ET=/tmp/akinya_614617/braker/F.oxysporum_fsp_fragariae_DSA14_003_braker_flye/hintsfile.gff --cores=1 --fungus --soft 1000 1>/tmp/akinya_614617/braker/F.oxysporum_fsp_fragariae_DSA14_003_braker_flye/GeneMark-ET.stdout 2>/tmp/akinya_614617/braker/F.oxysporum_fsp_fragariae_DSA14_003_braker_flye/errors/GeneMark-ET.stderr

    --GENEMARK_PATH=/home/gomeza/prog/genemark/gmes_linux_64 \
--BAMTOOLS_PATH=/home/gomeza/miniconda3/envs/gene_pred/bin \
BRAKER CRASHED afte 5 mins of editing paths

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

## Add gene prediction transcripts together

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
