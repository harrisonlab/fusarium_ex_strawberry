# Add path to profile PATH=${PATH}:/projects/oldhome/armita/prog/cufflinks/cufflinks-2.2.1.Linux_x86_64
# update profile - . ~/.profile

    for Alignment in $(ls alignment/star/M.domestica/GD/t1/GD_4A4/star_aligmentAligned.sortedByCoord.out.bam); do # star alignments output files
    Gff=apple_genome/gene_models_20170612.gff3 # Ref genome gff3
    OutDir=$(dirname $Alignment)
    mkdir -p $OutDir/fpkm
    cufflinks -p 8 -o $OutDir/fpkm -G $Gff $Alignment
