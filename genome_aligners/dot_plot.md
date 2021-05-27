# Create dotplot synteny analysis from two genomes
First create the minimap.paf file from 2 genomes

  conda activate olc_assemblers

## Step 3
Run like so using FULL PATHS
    Genome1=/projects/fusarium_ex_strawberry/NextGenSeq/repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_contigs_unmasked.fa       #replace with your genome
    Genome2=/projects/fusarium_ex_strawberry/NextGenSeq/Fo_Cepae/repeat_masked/Fus2_canu_new_contigs_unmasked.fa
    Strain=Fus2
    Outdir=/projects/fusarium_ex_strawberry/NextGenSeq/Genome_synteny/$Strain
    mkdir -p $Outdir
    Progdir=/home/akinya/git_repos/fusarium_ex_strawberry/genome_aligners
    sbatch $Progdir/minimap.sh $Genome1 $Genome2 $Strain $Outdir

    conda activate Ragtag
    conda install -c bioconda Ragtag
    conda install -c conda-forge r-optparse
    conda install -c conda-forge r-ggplot2
    conda install -c conda-forge r-plotly


## Create the dot plot
    for paf_file in $(ls -d /projects/fusarium_ex_strawberry/NextGenSeq/Genome_synteny/*/*_minimap.paf); do
    Strain=$(echo $paf_file | rev | cut -f2 -d '/' | rev)
    outdir=/projects/fusarium_ex_strawberry/NextGenSeq/Genome_synteny/$Strain
    echo $Strain
    Progdir=/home/akinya/git_repos/fusarium_ex_strawberry/genome_aligners
    sbatch $Progdir/create_dotplot.sh $paf_file $outdir
    done
Or
For create_dotplot.sh you need to change options to fit parameters of assemblies
    paf_file=/projects/fusarium_ex_strawberry/NextGenSeq/Genome_synteny/15-074/*_minimap.paf
    Strain=$(echo $paf_file | rev | cut -f2 -d '/' | rev)
    outdir=/projects/fusarium_ex_strawberry/NextGenSeq/Genome_synteny/$Strain
    echo $Strain
    Progdir=/home/akinya/git_repos/fusarium_ex_strawberry/genome_aligners
    sbatch $Progdir/create_dotplot.sh $paf_file $outdir

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

    #ragtag.py correct ref.fasta query.fasta
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
