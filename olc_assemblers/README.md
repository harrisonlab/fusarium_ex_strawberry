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

    ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
      for Assembly in $(ls Assembly2/flye/F.oxysporum_fsp_lactucae/race_1/pilon/pilon_1.fasta); do
        Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
        Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
        OutDir=$(dirname $Assembly)/ncbi_edits/round_1
        sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
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
    for Assembly in $(ls Assembly2/flye/F.oxysporum_fsp_lactucae/race_1/pilon/pilon_10.fasta); do
        OutDir=$(dirname $Assembly)
        $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/pilon_10_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
    done
    rm tmp.txt

Have 2 paths to choose from to run scripts if either doesn't work

### RepeatMask & TPSI path 1

    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking
    BestAssembly=Assembly2/flye/F.oxysporum_fsp_lactucae/race_1/pilon/pilon_10_renamed.fasta
    OutDir=repeat_masked/F.oxysporum_fsp_lactucae/race_1/flye/ncbi_edits_repmask
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

    for File in $(ls repeat_masked/F.oxysporum_fsp_lactucae/race_1/flye/ncbi_edits_repmask/race_1_contigs_softmasked.fa); do
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
    # flye - 9530627

## Hard Mask
Hard masking  means transforming every nucleotide identified as a repeat to an 'N' or 'X'.

    for File in $(ls repeat_masked/F.oxysporum_fsp_lactucae/race_1/flye/ncbi_edits_repmask/race_1_contigs_hardmasked.fa); do
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
    OutDir=Orthology/flye/blastn/$Organism/$Strain/FolR1vsFoCep_mimps
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
    sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly $OutDir
  done

Second query../F.oxysporum_fsp_cepae/Fus2_canu_new/final/Eff_mimp_genes.fasta
Use cut -f1 or * resulys_file.fasta_homologs.csv to excise and view column data
Compare against Andy's known SIX genes

  for Assembly in $(ls Assembly2/flye/F.oxysporum_fsp_lactucae/race_1/pilon/pilon_10_renamed.fasta); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    Query=../oldhome/groups/harrisonlab/project_files/fusarium/analysis/blast_homology/six_genes/six-appended_parsed.fa
    OutDir=Orthology/flye/blastn/$Organism/$Strain/FolR!vsFoLy
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
    sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly $OutDir
  done

# Synteny Check
#####################

## D-genies

Go to http://dgenies.toulouse.inra.fr/ to compare genomes for synteny against FoLy4287 and FoFrvsFoCep.
Using the plots, see which assembly will be best to use for gene_prediction using Fo_cepae data.

## Mummer

System for rapidly aligning entire genomes
run in conda env (mummer)
For help go to https://mummer4.github.io/manual/manual.html
Compare Folac against cepae genome and lycopersici genomes (reference genomes)
  SubjectGenome=reference genome
  QueryGenome= my genome i.e FolR1
  Prefix=$3
  OutDir=$4

Do for flye assembly/flye/*/*/pilon/pilon_10_renamed.fasta
assembly/SMARTdenovo/*/*/pilon/pilon_10_renamed.fasta
assembly/miniasm/*/*/pilon/pilon_10_renamed.fasta

Do for 4287_chromosomal also -../../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl_repmask/4287_chromosomal_contigs_unmasked.fa
  for SubjectGenome in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl_repmask/4287_chromosomal_contigs_unmasked.fa); do
    QueryGenome=Assembly2/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/pilon/pilon_10_renamed.fasta
    Organism=$(echo "$QueryGenome" | rev | cut -d '/' -f4 | rev)
    Strain=$(echo "$QueryGenome" | rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    Prefix=Fol_R1
    OutDir=alignment/mummer/SMARTdenovo/FoFrvsFoLyChr
    mkdir -p $OutDir
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
    sbatch $ProgDir/mummer.sh $SubjectGenome $QueryGenome $Prefix $OutDir
  done

To view particular hits on contig, do:
  less FoFr_14_coords.tsv | grep 'contig_11'

do for unmasked cepae genome

for SubjectGenome in $(ls ../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa); do
  QueryGenome=Assembly2/miniasm/F.oxysporum_fsp_lactucae/race_1/pilon/pilon_10_renamed.fasta
  Prefix=Fol_R1
  OutDir=alignment/mummer/miniasm/FoFrvsFoCep
  mkdir -p $OutDir
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
  sbatch $ProgDir/mummer.sh $SubjectGenome $QueryGenome $Prefix $OutDir
done

May produce errors in slurm out file therefore run this after
  #-c	Include percent coverage columns in the output
  #-l	Include sequence length columns in the output
  # -b	Brief output that only displays the non-redundant locations of aligning regions
  # -T	Switch output to tab-delimited format
  /scratch/software/mummer/mummer-4.0.0rc1/show-coords -c -l -b -T yourfiltereddeltafile.delta > coords.tsv


# Gene prediction
#####################

## STAR

Run in Repenv - condaenv
Fus2_CzapekDox, Fus2_GlucosePeptone, Fus2_PDA and Fus2_PDB are RNAseq data of infected onions
Only data samples that will map genes of F.oxy accurately
Need to concatenate data after STAR analysis
#--genomeSAindexNbases is unique to each genome and is 11 for FoFR

  for Assembly in $(ls repeat_masked/F.oxysporum_fsp_lactucae/race_1/flye/ncbi_edits_repmask/race_1_contigs_unmasked.fa);  do
      Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
      Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
      echo "$Organism - $Strain"
      FileF=../oldhome/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/F/*_trim.fq.gz
      FileR=../oldhome/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/R/*_trim.fq.gz
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

  Strain=race_1
    Organism=F.oxysporum_fsp_lactucae
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

    for Assembly in $(ls repeat_masked/F.oxysporum_fsp_lactucae/race_1/flye/ncbi_edits_repmask/race_1_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
        Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
        echo "$Organism - $Strain"
        OutDir=gene_pred/braker/$Organism/$Strain/flye
        AcceptedHits=alignment/star/F.oxysporum_fsp_lactucae/race_1/concatenated/concatenated.bam
        GeneModelName="$Organism"_"$Strain"_braker_flye
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

    for Assembly in $(ls repeat_masked/F.oxysporum_fsp_lactucae/race_1/flye/ncbi_edits_repmask/race_1_contigs_unmasked.fa); do
        Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
        Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev)
        echo "$Organism - $Strain"
        OutDir=gene_pred/stringtie/$Organism/$Strain/flye/concatenated_prelim
        mkdir -p $OutDir
        AcceptedHits=alignment/star/F.oxysporum_fsp_lactucae/race_1/concatenated/concatenated.bam
        ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/ProgScripts
        sbatch $ProgDir/stringtie.sh $AcceptedHits $OutDir
       done

## Codingquarry

Run in env with Python 2.7 (betaenv)
After first run, use cquarryV1
GFT file from stringtie/cufflinks output
my repo /home/akinya/git_repos/assembly_fusarium_ex/scripts
Antonio /home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction

  for Assembly in $(ls repeat_masked/F.oxysporum_fsp_lactucae/race_1/flye/ncbi_edits_repmask/race_1_contigs_unmasked.fa); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
      Organism=$(echo $Assembly| rev | cut -d '/' -f5 | rev)
      echo "$Organism - $Strain"
      OutDir=gene_pred/codingquary/$Organism/$Strain/flye
      mkdir -p $OutDir
      GTF=gene_pred/stringtie/F.oxysporum_fsp_lactucae/race_1/flye/concatenated_prelim/out.gtf
      ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/ProgScripts
      sbatch $ProgDir/codingquarry2.sh $Assembly $GTF $OutDir
    done

Takes about 4-6 hours depends on genome size

## Add gene prediction transcripts together

Additional transcripts - to be edited
Run in perly env (Repenv)
Type full paths, do not use asterisks
RUN LINE BY LINE AS IT WILL NOT WORK
Do segments one at a time for peace of mind

    BrakerGff=$(ls -d gene_pred/braker/F.oxysporum_fsp_lactucae/race_1/flye/F.oxysporum_fsp_lactucae_race_1_braker_flye/augustus.gff3)
    	Strain=$(echo $BrakerGff| rev | cut -d '/' -f4 | rev)
    	Organism=$(echo $BrakerGff | rev | cut -d '/' -f5 | rev)
    	echo "$Organism - $Strain"
    	Assembly=$(ls repeat_masked/F.oxysporum_fsp_lactucae/race_1/flye/ncbi_edits_repmask/race_1_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    	CodingQuarryGff=gene_pred/codingquary/F.oxysporum_fsp_lactucae/race_1/flye/out/PredictedPass.gff3
    	PGNGff=gene_pred/codingquary/F.oxysporum_fsp_lactucae/race_1/flye/out/PGN_predictedPass.gff3
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

For flye:
Braker: 19605 CQ: 1295 combined: 20900

## Gene renaming
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

    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
      for Genes in $(ls gene_pred/codingquary/F.oxysporum_fsp_lactucae/race_1/final/final_genes_appended_renamed.pep.fasta); do
        echo $Genes
        $ProgDir/interproscan.sh $Genes
      done 2>&1 | tee -a interproscan_submission.log

Interproscan: all jobs failed - couldn't run all jobs simultaneously
ERROR: uk.ac.ebi.interpro.scan.management.model.implementations.RunBinaryStep - Command line failed with exit code: 1
Need to run in batches - Had to run split DNA in sets of 10.
Split gene.pep.fasta like so:
  InFile=gene_pred/codingquary/F.oxysporum_fsp_lactucae/race_1/final/final_genes_appended_renamed.pep.fasta
  SplitDir=gene_pred/interproscan/$Organism/$Strain/flye
  InterproDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  InName=$(basename $InFile)
  mkdir -p $SplitDir
  $InterproDir/splitfile_500.py --inp_fasta $InFile --out_dir $SplitDir --out_base "$InName"_split

  for file in $(ls gene_pred/interproscan/F.oxysporum_fsp_lactucae/race_1/flye/*_split_2*); do
    sbatch /home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/Feature_annotation/run_interproscan.sh $file
    done

Need to merge interproscan output as follows

    ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/Feature_annotation
     for Proteins in $(ls gene_pred/codingquary/F.oxysporum_fsp_lactucae/race_1/final/final_genes_appended_renamed.pep.fasta); do
       Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
       Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
       echo "$Organism - $Strain"
       echo $Strain
       InterProRaw=gene_pred/interproscan/F.oxysporum_fsp_lactucae/race_1/raw/
       $ProgDir/append_interpro.sh $Proteins $InterProRaw
     done

Use this command to view particular features in interproscan data:
  less path/to/interproscan.tsv | grep 'gene feature' # e.g. transposon

## 2) SwissProt

SWISS-PROT is a curated protein sequence database which strives to provide a high level of annotation, a minimal level of redundancy and a high level of integration with other databases.
No requirements to run Swissprot
Uniprot databases are downloaded to /projects/dbUniprot

Intructions to create a database- do this first will save you a headache
    dbFasta=$(ls /projects/dbUniprot/swissprot_2020_June/uniprot_sprot.fasta)
    dbType="prot"
    Prefix="uniprot_sprot"
    makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
    #makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out gene_pred/swissprot/F.oxysporum_fsp_fragariae/DSA14_003/$Prefix.db

Now run swissprot

  for Proteome in $(ls gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/flye/final/final_genes_appended_renamed.pep.fasta); do
  Strain=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  Organism=$(echo $Proteome | rev | cut -f5 -d '/' | rev)
  OutDir=gene_pred/swissprot/$Organism/$Strain
  SwissDbDir=../../dbUniprot/swissprot_2020_June
  SwissDbName=uniprot_sprot
  ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/Feature_annotation
  sbatch $ProgDir/sub_swissprot_akin.sh $Proteome $OutDir $SwissDbDir $SwissDbName
  done

## 3) Signal-P
Need to install paths into project_files
  nano .profile # copy paths into profile
    PATH=${PATH}:/home/gomeza/prog/signalp-5.0b
    PATH=${PATH}:/data/scratch/gomeza/prog/signalp/signalp-5.0b/bin
    PATH=/data/scratch/gomeza/prog/java/jdk-11.0.4/bin:${PATH}
    PATH=${PATH}:/data/scratch/gomeza/prog/signalp/signalp-4.1
update your profile
  . ~/.profile

Signal P script for fungi
Add your strains name to first line
Added codingquary to Proteome direc

    for Strain in race_1; do
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/Feature_annotation
      CurPath=$PWD
      for Proteome in $(ls gene_pred/codingquary/F.oxysporum_fsp_lactucae/race_1/final/final_genes_combined.pep.fasta); do
      Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
      Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
      SplitDir=gene_pred/final_genes_split/$Organism/$Strain/flye
      mkdir -p $SplitDir
      BaseName="$Organism""_$Strain"_final_preds
      $ProgDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName # Split your input fasta in 500 genes files
        for File in $(ls $SplitDir/*_final_preds_*); do
        #sbatch $ProgDir/pred_signalP.sh $File signalp
        #sbatch $ProgDir/pred_signalP.sh $File signalp-3.0 # Recommended for oomycetes
        sbatch $ProgDir/pred_signalP.sh $File signalp-4.1 # Recommended for fungi
        #sbatch $ProgDir/pred_signalP.sh $File signalp-5.0
        done
      done
    done

Change output directory name to "final_genes_signalp-4.1"
  mv gene_pred/F.oxysporum_fsp_fragariae_signalp-4.1 gene_pred/final_genes_signalp-4.1

Need to combine the output of the first signal-P run
Make sure path to directories is correct!

  for Strain in race_1; do
   for SplitDir in $(ls -d gene_pred/final_genes_split/F.oxysporum_fsp_lactucae/$Strain); do
    Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
    Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
    InStringAA=''
    InStringNeg=''
    InStringTab=''
    InStringTxt=''
    SigpDir=final_genes_signalp-4.1
    for GRP in $(ls -l $SplitDir/*_final_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do
      InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/flye/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.aa";
      InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/flye/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp_neg.aa";
      InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/flye/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.tab";
      InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/flye/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.txt";
    done
    cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.aa
    cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_neg_sp.aa
    tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.tab
    cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.txt
   done
  done

Having flye in directory path caused small issues therefore I stopped including it from here
Things may be in the wrong directory - use "mv" command to change directory names

## 4) TMHMM step
Identifies transmembrane proteins
Added strain name
Add paths to .profile
  PATH=${PATH}:/data/scratch/gomeza/prog/tmhmm-2.0c/bin

  . ~/.profile # Refresh your profile

for Strain in race_1; do
	for Proteome in $(ls gene_pred/codingquary/F.oxysporum_fsp_lactucae/*/final/final_genes_combined.pep.fasta); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/Feature_annotation
		sbatch $ProgDir/TMHMM.sh $Proteome
	done
done

Proteins with transmembrane domains were removed from lists of Signal peptide containing proteins

    for File in $(ls gene_pred/trans_mem/F.oxysporum_fsp_lactucae/race_1/*_TM_genes_neg.txt); do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
    cat $File | cut -f1 > $TmHeaders
    SigP=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
    OutDir=$(dirname $SigP)
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
    $ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
    cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
    done

1544 proteins had transmembrane domains

## 5) EffectorP - Effector identification

Add path to .profile PATH=${PATH}:/scratch/software/EffectorP-2.0/Scripts
Use full paths to scripts - EffectorP is extremely picky with inputs
  # Make directory first
  mkdir -p analysis/effectorP/$Organism/$Strain/flye
Note down your paths
  Basename="$Organism"_"$Strain"_EffectorP
  Proteome=$(ls -d gene_pred/codingquary/F.oxysporum_fsp_lactucae/race_1/final/final_genes_appended_renamed.pep.fasta)
  OutDir=analysis/effectorP/F.oxysporum_fsp_lactucae/race_1/flye
  EffectorP.py -o analysis/effectorP/F.oxysporum_fsp_lactucae/race_1/flye/F.oxysporum_fsp_lactucae_race_1_EffectorP.txt -E analysis/effectorP/F.oxysporum_fsp_lactucae/race_1/flye/F.oxysporum_fsp_lactucae_race_1_EffectorP.fa -i gene_pred/codingquary/F.oxysporum_fsp_lactucae/race_1/final/final_genes_appended_renamed.pep.fasta

### EffectorP - phase 2
  for File in $(ls analysis/effectorP/F.oxysporum_fsp_lactucae/race_1/flye/F.oxysporum_fsp_lactucae_race_1_EffectorP.txt); do
    Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
    cat $File | grep 'Effector' | cut -f1 > $Headers
    Secretome=$(ls gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_lactucae/race_1/race_1_final_sp_no_trans_mem.aa)
    OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
    $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
    OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
    cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
    cat $OutFileHeaders | wc -l
    Gff=$(ls gene_pred/codingquary/F.oxysporum_fsp_lactucae/race_1/final/final_genes_appended_renamed.gff3)
    EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
    $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
    cat $EffectorP_Gff | grep -w 'gene' | wc -l
  done > tmp.txt

  nano tmp.txt - should state Org, Strain and Numbers
169

## 6) Mimp analysis
Miniature IMPala elements are short autonomous class II  transposable elements and can be used to identify candidate effectors
Run in conda env (Repenv)

  for Assembly in $(ls repeat_masked/F.oxysporum_fsp_lactucae/race_1/flye/ncbi_edits_repmask/race_1_contigs_unmasked.fa); do
      Organism=$(echo "$Assembly" | rev | cut -d '/' -f5 | rev)
      Strain=$(echo "$Assembly" | rev | cut -d '/' -f4 | rev)
      GeneGff=$(ls gene_pred/codingquary/F.oxysporum_fsp_lactucae/race_1/final/final_genes_appended_renamed.gff3)
      OutDir=analysis/mimps/$Organism/$Strain
      mkdir -p "$OutDir"
      echo "$Organism - $Strain"
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
      $ProgDir/mimp_finder.pl $Assembly $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
      $ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
      echo "The number of mimps identified:"
      cat $OutDir/"$Strain"_mimps.fa | grep '>' | wc -l
      bedtools intersect -u -a $GeneGff -b $OutDir/"$Strain"_mimps_exp.gff > $OutDir/"$Strain"_genes_in_2kb_mimp.gff
      echo "The following transcripts intersect mimps:"
      MimpProtsTxt=$OutDir/"$Strain"_prots_in_2kb_mimp.txt
      MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
      cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | sort | uniq > $MimpProtsTxt
      cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | cut -f1 -d '.'| sort | uniq > $MimpGenesTxt
      cat $MimpProtsTxt | wc -l
      cat $MimpGenesTxt | wc -l
      echo ""
    done

The number of mimps identified:
211
The following transcripts intersect mimps:
175
175

## 7) Cazy

CAZY analysis - program for searching and analyzing carbohydrate-active enzymes in a newly sequenced organism using CAZy database
"sub_hmmscan.sh" was edited and updated to "hmmscan.sh"
make sure your directories are all correct

    for Strain in race_1; do
    for Proteome in $(ls gene_pred/codingquary/F.oxysporum_fsp_lactucae/$Strain/final/final_genes_combined.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/CAZY/$Organism/$Strain
    mkdir -p $OutDir
    Prefix="$Strain"_CAZY
    CazyHmm=../dbCAN/dbCAN-fam-HMMs.txt
    ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/Feature_annotation
    sbatch $ProgDir/hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
    done
    done

### CAZy phase 2

Run line by line
Creates a file with CAZy module and gene

    for File in $(ls gene_pred/CAZY/F.oxysporum_fsp_lactucae/race_1/race_1_CAZY.out.dm); do
      Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
      Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
      OutDir=$(dirname $File)
      echo "$Organism - $Strain"
      ProgDir=../dbCAN
      $ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
      CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
      cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders # Extract gene names
      echo "Number of CAZY genes identified:"
      cat $CazyHeaders | wc -l
      Gff=$(ls gene_pred/codingquary/F.oxysporum_fsp_lactucae/race_1/final/final_genes_appended_renamed.gff3)
      CazyGff=$OutDir/"$Strain"_CAZY.gff
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
      $ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff # Creates a gff for all CAZymes
      SecretedProts=$(ls gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_lactucae/race_1/race_1_final_sp_no_trans_mem.aa)
      SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
      cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
      CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
      $ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted # Creates a gff for secreted CAZymes
      echo "Number of Secreted CAZY genes identified:"
      cat $CazyGffSecreted | grep -w 'gene' | cut -f9 | tr -d 'ID=' | wc -l
      done

Number of CAZY genes identified:
954
Number of Secreted CAZY genes identified:
407

## 8) Antismash

Antismash was run to identify clusters of secondary metabolite genes within the genome. Antismash was run using the webserver at: http://antismash.secondarymetabolites.org.
Use you polished unmasked contig genome as the genome input "*_contigs_unmasked.fa" and the complementary "final_genes_appended_renamed.gff3" for the gene input.
Download antiSMASH results

    mkdir -p gene_pred/antiSMASH/F.oxysporum_fsp_lactucae/race_1/
    cd # to new directory
    wget download.link.com/results.zip # not actual result file
    #Then unzip file
    unzip results.zip

Results of web-annotation of gene clusters within the assembly were downloaded to the following directories
Run in conda env ( e.g. Repenv)
Run line by line

    for AntiSmash in $(ls -d gene_pred/antiSMASH/F.oxysporum_fsp_fragariae/DSA14_003/DSA14_003_contigs_unmasked.gbk); do
      Organism=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $AntiSmash | rev | cut -f2 -d '/' | rev)
      echo "$Organism - $Strain"
      OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
      Prefix=$OutDir/"$Strain"_antismash_results
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
      # Antismash v5 output to gff file
      $ProgDir/antismash2gffv5.py --inp_antismash $AntiSmash --out_prefix $Prefix
      #$ProgDir/antismash2gff.py --inp_antismash $AntiSmash --out_prefix $Prefix # Use only for antismash v4.2 output
      printf "Number of secondary metabolite detected:\t"
      cat "$Prefix"_secmet_clusters.gff | wc -l
      GeneGff=gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/flye/final/final_genes_appended_renamed.gff3
      bedtools intersect -u -a $GeneGff -b "$Prefix"_secmet_clusters.gff > "$Prefix"_secmet_genes.gff
      cat "$Prefix"_secmet_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_secmet_genes.txt
      bedtools intersect -wo -a $GeneGff -b "$Prefix"_secmet_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -e "s/;Parent=g\w+//g" | perl -p -e "s/;Notes=.*//g" > "$Prefix"_secmet_genes.tsv
      printf "Number of predicted proteins in secondary metabolite clusters:\t"
      cat analysis/secondary_metabolites/antismash/F.oxysporum_fsp_fragariae/DSA14_003/*_secmet_genes.txt | wc -l     
      printf "Number of predicted genes in secondary metabolite clusters:\t"
      cat "$Prefix"_secmet_genes.gff | grep -w 'gene' | wc -l
    done

Results
  Number of secondary metabolite detected:        57
  Number of predicted proteins in secondary metabolite clusters:  1412
  Number of predicted genes in secondary metabolite clusters:     704

Antismash output correction. Some gene names contain ;. Remove manually with the following command.
First sed command removes ;. Second and Third remove the cluster kind information (optional)
  cat analysis/secondary_metabolites/antismash/F.oxysporum_fsp_fragariae/DSA14_003/DSA14_003_antismash_results_secmet_genes.tsv | sed 's/;//p' | sed 's/;.*//p' | sed 's/Kin.*//p' > analysis/secondary_metabolites/antismash/F.oxysporum_fsp_fragariae/DSA14_003/DSA14_003_antismash_results_secmet_genes_corrected.tsv
Edit output file names from this script after completion
