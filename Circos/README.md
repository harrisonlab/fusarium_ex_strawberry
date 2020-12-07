README
# Script to create circos plots.
1.) First create txt files for genomes and concatenate
Text file 1 is first genome.Text file 2 is second genome. Text files have contig name followed by the size
2.) Create synteny links between genomes using satsuma synteny
3.) Create circos config file used for building ideogram.


## 1.) Create a .txt file with contig lengths of the genomes to be aligned.
txt file for genome 1
Run in a conda env. Run line by line
  conda create -n Circos python
  conda install -c conda-forge biopython
  conda install -c bioconda perl-bioperl
  conda install -c nanjiang satsuma

  python
  >>> import Bio
  >>> import sets
  >>> press "ctrl + D" to exit python

Make file executable before running
add "import Bio" line to top of python_create_circos_file.py program just under the title if it is not there already
    chmod u+x ./python_create_circos_file.py

Or just run in alphaenv

    Genome1=repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/miniasm/ncbi_edits_repmask/DSA14_003_contigs_unmasked.fa
    #Genome3=repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_contigs_unmasked.fa
    #Genome4=repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/SMARTdenovo/ncbi_edits_repmask/DSA14_003_contigs_unmasked.fa
    OutDir=synteny_analysis/circos/flye
    mkdir -p $OutDir
    ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/Circos
    $ProgDir/python_create_circos_file.py --genome $Genome4 --contig_prefix "FoFr_14_S_" > $OutDir/FoFr_14_Sden_genome.txt
    #FoFr_14_F=flye FoFr_14_S= Sden

### .txt file for genome 2

    Genome2=../../oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl_repmask/4287_chromosomal_contigs_unmasked.fa
    OutDir=synteny_analysis/circos/flye
    ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/Circos
    $ProgDir/python_create_circos_file.py --genome $Genome2 --contig_prefix "Fol_4287_Chr_" > $OutDir/FoLy_illumina_genome.txt

### Concatenate files

Do for all txt files in respective directories: FoFr_14_flye_genome.txt FoFr_14_Sden_genome.txt FoFr_14_ONT_genome.txt

  File=synteny_analysis/circos/SMARTdenovo
  cat $File/FoFr_14_Sden_genome.txt > $File/FoFr_FoLy_genome.txt
  tac $File/FoLy_illumina_genome.txt >> $File/FoFr_FoLy_genome.txt

 #Contigs smaller than 10Kb can be removed if not required.
 eg...
  cat $File/Fv_Fv_genome.txt \
  | grep -v -e "A3_5_contig_87" \
  | grep -v -e "A3_5_contig_88" \
  | grep -v -e "A3_5_contig_89" \
  | grep -v -e "A3_5_contig_90" \
  | grep -v -e "A3_5_contig_91" \
  | grep -v -e "A3_5_contig_92" \
  | grep -v -e "A3_5_contig_93" \
  | grep -v -e "A3_5_contig_94" \
  | grep -v -e "A3_5_contig_95" \
  | grep -v -e "A3_5_contig_96" \
  | grep -v -e "A3_5_contig_97" \
  | grep -v -e "A3_5_contig_98" \
  | grep -v -e "A3_5_contig_99" \
  | grep -v -e "A3_5_contig_100" \
  | grep -v -e "A3_5_contig_101" \
  | grep -v -e "A3_5_contig_102" \
  | grep -v -e "A3_5_contig_103" \
  | grep -v -e "A3_5_contig_104" \
  | grep -v -e "A3_5_contig_105" \
  | grep -v -e "A3_5_MINcontig_5" \
  > $OutDir/Fv_Fv_genome_edited.txt


## 2.) Create synteny links between genomes.

Type "SatsumaSynteny" to check if it is installed in conda env.

Note that SatsumaSynteny calls other executables (FilterGridSeeds, HomologyByXCorr, HomologyByXCorrSlave, MergeXCorrMatches), and thus has to be invoked by either supplying the full path of the executable, or
“./SatsumSynteny”

USE FULL PATHS AS YOU RUN PROG FROM HOME DIRECTORY!!!!!!!!!!!!!!

  mkdir SatsumaSynteny
  cd SatsumaSynteny/
  wget ftp://ftp.broadinstitute.org/distribution/software/spines/satsuma-3.0.tar.gz
  gunzip satsuma-3.0.tar.gz
  tar -xf satsuma-3.0.tar
  cd satsuma-code-0
  # make a files executable
  chmod u+x ./M*
  chmod u+x ./S*
  chmod u+x ./t*
  chmod u+x ./H*
  chmod u+x ./F*
  chmod u+x ./C*
  chmod u+x ./B*
  chmod u+x ./c*
  chmod u+x ./T*
  # copy directory into "satsuma_synteny.sh" place where command is

    /home/*/SatsumaSynteny/satsuma-code-0/


  Genome1=/projects/fusarium_ex_strawberry/NextGenSeq/repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/miniasm/ncbi_edits_repmask/DSA14_003_contigs_unmasked.fa
  #Genome3=/projects/fusarium_ex_strawberry/NextGenSeq/repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_contigs_unmasked.fa
  #Genome4=/projects/fusarium_ex_strawberry/NextGenSeq/repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/SMARTdenovo/ncbi_edits_repmask/DSA14_003_contigs_unmasked.fa
  Genome2=/projects/oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lycopersici/4287_v2/fungidb_repmask/4287_v2_contigs_unmasked.fa
  #Genome5=/projects/oldhome/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl_repmask/4287_chromosomal_contigs_unmasked.fa
  OutDir=/projects/fusarium_ex_strawberry/NextGenSeq/synteny_analysis/circos/miniasm/satsuma_alignment/FoFr_14VSFoLy
  mkdir -p $OutDir
  ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/Circos
  sbatch $ProgDir/satsuma_synteny.sh $Genome1 $Genome2 $OutDir



## 3.) Edit contig prefix to match step 1, this is required as contig prefix must remain uniform.
This script assumes your contigs are called "contig_". Change pre1 and pre2 to new contig prefix identifier.
Changes output of satsuma contig names - major key as it distinguishes information


    Synteny_file=synteny_analysis/circos/flye/satsuma_alignment/FoFr_14VSFoLy/satsuma_summary.chained.out           #Synteny file location and file name here
    OutDir=synteny_analysis/circos/flye/satsuma_alignment/FoFr_14VSFoLy/satsuma_summary_editedforcircos.chained.out  #Outfile location and file name here
    pre1=FoFr_14_contig_  #change to your contig prefix
    pre2=Fol_4287_Chr_ #change to your contig prefix
    #ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/Circos/
    sbatch $ProgDir/Python_edit_contig_names.py --synteny $Synteny_file --contig_prefix_1 $pre1 --contig_prefix_2 $pre2 --outfile $OutDir

Had to edit line 34 in the "Python_edit_contig_names.py" script
  # From "print out_line" to "print (out_line)"


#3.) Run circos
# Circos relies on a configuration file which has other branching files involved in the ideogram configuration.
# The 2D plot must link to satsuma_summary_chained_out file from previous step


Conf=/home/connellj/git_repos/emr_repos/Fv_C-variants/Circos/circos_configuration.sh
OutDir=/home/connellj/Circos/ideogram
mkdir -p $OutDir
ProgDir=/home/connellj/git_repos/emr_repos/Fv_C-variants/Circos
sbatch $ProgDir/circos.sh $Conf $OutDir
