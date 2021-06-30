# Genome aligners
## STAR
Results show this is best STAR script - it worked :)
Run in Repenv - condaenv
    conda activate Repenv
Fus2_CzapekDox, Fus2_GlucosePeptone, Fus2_PDA and Fus2_PDB are RNAseq data of infected onions
Only data samples that will map genes of F.oxy accurately
Need to concatenate data after STAR analysis.

    for Assembly in $(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA15_041/ncbi_edits_repmask/DSA15_041_contigs_unmasked.fa)
      do
        Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
        Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
        echo "$Organism - $Strain"
        FileF=../oldhome/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/F/*_trim.fq.gz
        FileR=../oldhome/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/R/*_trim.fq.gz
        echo $FileF
        echo $FileR
        Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
        echo "$Timepoint"
        Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
        OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
        ProgDir=/home/akinya/git_repos/assembly_fusarium_ex/ProgScripts
        sbatch $ProgDir/star.sh $Assembly $FileF $FileR $OutDir 11
      done
    done

View "star_aligmentLog.final.out" to see uniquely mapped reads %

      Strain=race_1
        Organism=F.oxysporum_fsp_lactucae
        mkdir -p alignment/star/$Organism/$Strain/concatenated
        samtools merge -f alignment/star/$Organism/$Strain/concatenated/concatenated.bam \
        alignment/star/$Organism/$Strain/Fus2_CzapekDox/6_S2_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam \
        alignment/star/$Organism/$Strain/Fus2_GlucosePeptone/7_S3_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam \
        alignment/star/$Organism/$Strain/Fus2_PDA/9_S4_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam \
        alignment/star/$Organism/$Strain/Fus2_PDB/4_S1_L001_R1_001_trim.fq.gz/star_aligmentAligned.sortedByCoord.out.bam

This following section is somewhat optional.
Dotplot and mummer basically do the same thing, align genomes to generate coordinates which can be used to form a dot plot of the aligned genomes.

So here goes...

## Mummer
System for rapidly aligning entire genomes
Run in conda env (mummer)
    conda activate mummer
For help go to https://mummer4.github.io/manual/manual.html

    SubjectGenome=reference genome
    QueryGenome= my genome i.e FoFr_14
    Prefix=$3
    OutDir=$4

I did this for the assembly of the long read FoFr_14 sequence. Compared it to the highly studied & characterised F.oxysporum_fsp_lycopersici 4827
Location:
    projects/oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl_repmask/4287_chromosomal_contigs_unmasked.fa
The command line entry is as follows:

    for SubjectGenome in $(ls ../../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl_repmask/4287_chromosomal_contigs_unmasked.fa); do
      QueryGenome=assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta
      Organism=$(echo "$QueryGenome" | rev | cut -d '/' -f4 | rev)
      Strain=$(echo "$QueryGenome" | rev | cut -d '/' -f3 | rev)
      echo "$Organism - $Strain"
      Prefix=FoFr_14
      OutDir=alignment/flye/miniasm/FoFrvsFoLy2
      mkdir -p $OutDir
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
      sbatch $ProgDir/mummer.sh $SubjectGenome $QueryGenome $Prefix $OutDir
    done

To view particular hits on contig, do:
    less FoFr_14_coords.tsv | grep 'contig_11'

I repeated the steps for the unmasked cepae genome Fus2

    for SubjectGenome in $(ls ../../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa); do
      QueryGenome=assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta
      Prefix=FoFr_14
      OutDir=alignment/mummer/flye/FoFrvsFoCep
      mkdir -p $OutDir
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
      sbatch $ProgDir/mummer.sh $SubjectGenome $QueryGenome $Prefix $OutDir
    done

May produce errors in slurm out file therefore run with these options after
    -c	Include percent coverage columns in the output
    -l	Include sequence length columns in the output
    -b	Brief output that only displays the non-redundant locations of aligning regions
    -T	Switch output to tab-delimited format
 /scratch/software/mummer/mummer-4.0.0rc1/show-coords -c -l -b -T yourfiltereddeltafile.delta > coords.tsv

## Dot plots

Create dotplot synteny analysis from two genomes
First create the minimap.paf file from 2 genomes

    Genome1=/path/to/genomes            #replace with your genome
    Genome2=/path/to/genomes  #replace with your genome
    Outdir=/projects/fusarium_ex_strawberry/*/synteny
    #mkdir -p $Outdir
    Progdir=/home/connellj/git_repos/emr_repos/Fv_C-variants/Dotplot_analysis/
    sbatch $Progdir/minimap.sh $Genome1 $Genome2 $Outdir

Create the dot plot.
Does this by running R in command line and generates a pretty plot
    paf_file=/path/to/minimap.paf;
    Outdir=/projects/*/genome_synteny/
    Progdir=/home/connellj/git_repos/emr_repos/Fv_C-variants/Dotplot_analysis/
    sbatch $Progdir/create_dotplot.sh $paf_file $Outdir

## ragtag

Firstly create an env.
    conda create -n Ragtag python
    conda install -c bioconda ragtag

Once installed, the path to ragtag is as follows:
    /home/user/miniconda3/pkgs/ragtag-1.1.1-pyh7b7c402_0/python-scripts/ragtag.py

You can run this in a few ways:

    # you can correct your contigs
      ragtag.py correct ref.fasta query.fasta

    # scaffold contigs
      ragtag.py scaffold ref.fa ragtag_output/query.corrected.fasta

    # scaffold with multiple references
    ragtag.py scaffold -o out_1 ref1.fasta query.fasta
    ragtag.py scaffold -o out_2 ref2.fasta query.fasta
    ragtag.py merge query.fasta out_*/*.agp

    # use Hi-C to resolve conflicts
    ragtag.py merge -b hic.bam query.fasta out_*/*.agp

The running script is like this
    reference_genome=/projects/path/to/reference/genome
    query_genome=/projects/path/to/query/genome
    Outdir=/projects/*/ragtag_output
    mkdir -p $Outdir
    progdir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/Genome_alignment
    sbatch $progdir/ragtag.sh $reference_genome $query_genome $Outdir

    or
    ragtag.py scaffold -u -C DSA14_003_contigs_NGS_unmasked.fa DSA14_003_contigs_unmasked.fa
    ragtag.py correct -u -C DSA14_003_contigs_NGS_unmasked.fa DSA14_003_contigs_unmasked.fa

Chr0_RagTag produced by ragtag are reads that could not be aligned to reference. May have to manually map them
  so...
  # create a txt file with contig name
  nano Chr0_RagTag_14.txt
Cut the DNA sequence out into a new file
  faidx -d '|' DSA14_003.scaffolds.fasta $(tr '\n' ' ' < Chr0_RagTag_14.txt) > Chr0_RagTag_14.fasta
