# Create dotplot synteny analysis from two genomes
First create the minimap.paf file from 2 genomes

## Step 1 install library in command line
    pip install optparse-pretty

## Step 2 install in R
    install.packages("optparse")

## Step 3
Run like so
    Genome1=NextGenSeq/repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_contigs_unmasked.fa        #replace with your genome
    Genome2=repeat_masked/F.oxysporum_fsp_fragariae/DSA15_041/ncbi_edits_repmask/DSA15_041_contigs_unmasked.fa
    Outdir=/projects/fusarium_ex_strawberry/NextGenSeq/Genome_synteny
    #mkdir -p $Outdir
    Progdir=/home/connellj/git_repos/emr_repos/Fv_C-variants/Dotplot_analysis/
    sbatch $Progdir/minimap.sh $Genome1 $Genome2 $Outdir

## Create the dot plot
    paf_file=/projects/fusarium_ex_strawberry/NextGenSeq/Genome_synteny/minimap.paf
    outdir=/projects/fusarium_ex_strawberry/NextGenSeq/Genome_synteny
    Progdir=/home/akinya/git_repos/fusarium_ex_strawberry/genome_aligners
    sbatch $Progdir/create_dotplot.sh $paf_file $outdir

error message:cause still  unknown
    Error in library(optparse) : there is no package called ‘optparse’
    Calls: suppressPackageStartupMessages -> withCallingHandlers -> library
    Execution halted
    cp: cannot stat '/home/akinya/akinya_749626/paf_plot.png': No such file or directory

# Ragtag

Install in conda env
    conda activate Ragtag
    conda install -c bioconda ragtag

Next run it
    reference_genome=/projects/botrytis_cinerea/External_group/renamed/GCF_000143535.2_ASM14353v4_genomic.fasta
    query_genome=/projects/botrytis_cinerea/External_group/sl9/renamed/JACVFN01.1.fasta
    Outdir=/projects/botrytis_cinerea/ragtag_output
    mkdir -p $Outdir
    progdir=/home/bourquif/git_repos/scripts/genome_correction
    sbatch $progdir/ragtag.sh $reference_genome $query_genome $Outdir

Or run like This
    screen -a
    srun --partition himem --time 0-06:00:00 --mem-per-cpu 20G --cpus-per-task 24 --pty bash

    ragtag.py correct ref.fasta query.fasta
    ragtag.py scaffold ref.fasta query.fasta

assembly improvement:
correct         misassembly correction
scaffold        synteny scaffolding
merge           scaffold merging
patch           continuous scaffolding & gap filling

file utilities:
agp2fa          build a FASTA file from an AGP file
agpcheck        check for valid AGP file format
asmstats        assembly statistics
splitasm        split an assembly at gaps
delta2paf       delta to PAF file conversion
paf2delta       PAF to delta file conversion
updategff       update gff intervals
