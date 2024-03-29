# Nanoplot
Steps to take to generate plots for  Nanopore data

Create & activate your environment
    conda create -n Nanoplot
    conda activate Nanoplot

Install Nanoplot
    conda install -c bioconda nanoplot

NanoPlot creates:
A statistical summary
A number of plots
A html summary file

Usage; examples
    NanoPlot -h
    Nanoplot --summary sequencing_summary.txt --loglength -o summary-plots-log-transformed  
    NanoPlot -t 2 --fastq reads1.fastq.gz reads2.fastq.gz --maxlength 40000 --plots dot --legacy hex
    NanoPlot -t 12 --color yellow --bam alignment1.bam alignment2.bam alignment3.bam --downsample 10000 -o bamplots_downsampled

    #--plots uses the plotly package to plot kde and dot plots.
    #--maxlength N         Hide reads longer than length specified.
    #--minlength N         Hide reads shorter than length specified
    #-o, --outdir OUTDIR   Specify directory in which output has to be created.
    # --loglength           Logarithmic scaling of lengths in plots.
    # -c, --color COLOR     Specify a color for the plots, must be a valid matplotlib color
    #--N50                 Show the N50 mark in the read length histogram
    #--title TITLE         Add a title to all plots, requires quoting if using space
    #--summary file        Data is in one or more summary file(s) generated by albacore or guppy
    #--fastq               Data is in one or more default fastq file
    #--legacy plotting of a hex plot currently is only possible using this option,which uses the seaborn and matplotlib package
