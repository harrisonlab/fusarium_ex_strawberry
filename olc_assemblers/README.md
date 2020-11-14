# Use for long read assembly programs

# Method was used for Fusarium oxysporum fsp lactucae
#### 1 # , then you have a header, ## subheading , ### small heading - use "tab" key to identify script info

####################
# Miniasm assembly
####################

## Step 1
#Login to node srun --partition himem --time 0-06:00:00 --mem 40G --cpus-per-task 24 --pty bash
Concatenate sequence reads first (there is old seq data that was basecalled again include it)
#Use command below if you are working in the same directory as the raw sequence reads

    cat *fastq | gzip -cf > FAL69458.fastq.gz
    /archives/2020_niabemr_nanopore/F.oxyspporum_lactucae_Race1/20180426_AJ520_GA30000$ cat *.fastq.gz | gzip -cf > /projects/fusarium_EX_Lactucae/raw_dna/FoLR12018NBC.fastq.gz
    cat *.fastq.gz | gzip -cf > /projects/fusarium_EX_Lactucae/raw_dna/FoLR1-18nbc.fastq.gz
    cat FAL69458.fastq.gz FoLR1-18nbc.fastq.gz FoLR12018NBC.fastq.gz > FolR1_conc.fastq.gz

#Method above caused error in porechop (died 3 times without running possibly due to concatenating 3 big gz files???)

Tried method 2 - gave same error
    cat 20180426_AJ520_GA30000/*.fastq.gz basecalling-all/*.fastq.gz fastq_pass/*.fastq | gzip -cf > /projects/fusarium_EX_Lactucae/raw_dna/FoLr1cont.fastq.gz

## Step 2
#Run Porechop before assembly

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


## Step 3
#Need to rename all reads

    rename.sh qin=33 in=FAL_trim.fastq.gz out=trimmed_renamed.fasta prefix=FolR1

#If script doesn't work, see below

        for TrimReads in $(ls FAL_trim.fastq.gz); do
            Organism=F.oxysporum_fsp_lactucae
            Strain=race_1
            Prefix="$Strain"_miniasm
            OutDir=assembly/miniasm/$Organism/$Strain
            mkdir -p $OutDir
            ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers
            sbatch $ProgDir/miniasm.sh $TrimReads $Prefix $OutDir
          done

## Step 4
#Run minimap2
Need to run read against itself. IT IS NEEDED. it is going to do self mapping
#Fast all-against-all overlap of raw reads

    minimap2 -x ava-ont -t8 trimmed_renamed.fasta trimmed_renamed.fasta | gzip -1 > FolR1_fastq_allfiles.paf.gz

#For ONT long read sequences use miniasm to assemble genome
Run in a screen and a node in a conda env with miniasm installed
#Run in olc_assemblers conda shell

## Step 5
#Concatenate pieces of read sequences to generate the final sequences
Can run like this instead: miniasm -f trimmed_renamed.fasta FolR1_fastq_allfiles.paf.gz > reads.gfa

    miniasm -f "$Prefix"_rename.fasta $Prefix.paf.gz > reads.gfa

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

    for TrimReads in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz); do
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

#Updated entire script using https://github.com/harrisonlab/bioinformatics_tools/blob/master/Gene_prediction/README.md
#Look into BuscoDB direc - directory exists
#Run in conda env - BUSCOenv

    for Assembly in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/assembly.fasta); do
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
        ReadsFq=$(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz)
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

#Ended up here for some reason assembly/flye/F.oxysporum_fsp_lactucae/race_1/ncbi_edits/round_*

    for Assembly in $(ls assembly/flye/F.oxysporum_fsp_lactucae/race_1/flye_raw/racon_10/assembly_racon_round_1.fasta); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
      Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
      echo "$Organism - $Strain"
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10/round_1
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done

#Contigs must be renamed before medaka can be run
Rename contigs for genome
If split or remove contigs is needed, provide FCSreport file by NCBI.

    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Assembly_qc
        touch tmp.txt
        for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/racon_10/race_1_smartdenovo_racon_round_1.fasta); do
            OutDir=$(dirname $Assembly)
            $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/race_1_smartdenovo_racon_round_1_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
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

Run QC checks with QUAST and BUSCO again to see any changes

## Pilon
Automatically improves draft assemblies and find variation among strains, including large event detection
Run in conda env (olc_assemblers)

    for Assembly in $(ls race_1_smartdenovo_racon_round_10_renamed.fasta); do
      Organism=F.oxysporum_fsp_lactucae
      Strain=AJ520
      IlluminaDir=$(ls -d $Strain)
      echo $Strain
      echo $Organism
      TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n2 | tail -n1);
      TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n2 | tail -n1);
      echo $TrimF1_Read
      echo $TrimR1_Read
      OutDir=$(dirname $Assembly)
      Iterations=10
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_assemblers/pilon
      sbatch $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
    done
