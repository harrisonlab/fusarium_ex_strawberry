#Create dotplot synteny analysis from two genomes
# First create the minimap.paf file from 2 genomes

Genome1=/path/to/genomes            #replace with your genome
Genome2=/path/to/genomes  #replace with your genome
Outdir=/projects/fusarium_ex_strawberry/*/synteny
#mkdir -p $Outdir
Progdir=/home/connellj/git_repos/emr_repos/Fv_C-variants/Dotplot_analysis/
sbatch $Progdir/minimap.sh $Genome1 $Genome2 $Outdir

#Create the dot plot
paf_file=/projects/botrytis_cinerea/genome_synteny/minimap.paf;
Outdir=/projects/botrytis_cinerea/genome_synteny/
Progdir=/home/connellj/git_repos/emr_repos/Fv_C-variants/Dotplot_analysis/
sbatch $Progdir/create_dotplot.sh $paf_file $Outdir
