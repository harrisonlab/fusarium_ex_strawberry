#Create dotplot synteny analysis from two genomes
# First create the minimap.paf file from 2 genomes

# step 1 install library in command line
    pip install optparse-pretty

# step 2 install in R
#install.packages("optparse")


Genome1=NextGenSeq/repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_contigs_unmasked.fa        #replace with your genome
Genome2=repeat_masked/F.oxysporum_fsp_fragariae/DSA15_041/ncbi_edits_repmask/DSA15_041_contigs_unmasked.fa
Outdir=/projects/fusarium_ex_strawberry/NextGenSeq/Genome_synteny
#mkdir -p $Outdir
Progdir=/home/connellj/git_repos/emr_repos/Fv_C-variants/Dotplot_analysis/
sbatch $Progdir/minimap.sh $Genome1 $Genome2 $Outdir

#Create the dot plot
paf_file=/projects/fusarium_ex_strawberry/NextGenSeq/Genome_synteny/minimap.paf
outdir=/projects/fusarium_ex_strawberry/NextGenSeq/Genome_synteny
Progdir=/home/connellj/git_repos/emr_repos/Fv_C-variants/Dotplot_analysis/
sbatch $Progdir/create_dotplot.sh $paf_file $outdir

Error in library(optparse) : there is no package called ‘optparse’
Calls: suppressPackageStartupMessages -> withCallingHandlers -> library
Execution halted
cp: cannot stat '/home/akinya/akinya_749626/paf_plot.png': No such file or directory
