# System for rapidly aligning entire genomes
# run in conda env (mummer)
# for help go to http://mummer.sourceforge.net/manual/
# compare Fofr14 against cepae genome and lycopersici genomes (reference genomes)
#SubjectGenome=reference genome
#QueryGenome= my genome i.e FoFr_14
#Prefix=$3
#OutDir=$4


for SubjectGenome in $(ls ../../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lycopersici/4287_v2/fungidb_repmask/4287_v2_contigs_unmasked.fa); do
  QueryGenome=assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta
  Organism=$(echo "$QueryGenome" | rev | cut -d '/' -f4 | rev)
  Strain=$(echo "$QueryGenome" | rev | cut -d '/' -f3 | rev)
  echo "$Organism - $Strain"
  Prefix=FoFr_14
  OutDir=alignment/mummer/miniasm/FoFrvsFoLy
  mkdir -p $OutDir
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
  sbatch $ProgDir/mummer.sh $SubjectGenome $QueryGenome $Prefix $OutDir
done

# to view particular hits on contig, do
# less FoFr_14_coords.tsv | grep 'contig_11'

 # may produce errors in slurm out file therefore run this after
# -c	Include percent coverage columns in the output
# -l	Include sequence length columns in the output
# -b	Brief output that only displays the non-redundant locations of aligning regions
# -T	Switch output to tab-delimited format
 /scratch/software/mummer/mummer-4.0.0rc1/show-coords -c -l -b -T yourfiltereddeltafile.delta > coords.tsv
