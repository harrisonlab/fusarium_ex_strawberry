### TE Hunt
# If you haven't already, you need to:
## Repeat mask
Ensure packages are installed in envs

    conda create -n RMask
    conda install -c bioconda repeatmodeler # repeatmodeler also installs packages below
    #conda install -c bioconda repeatmasker
    #conda install rmblast

Need to manually configure Repeatmasker

    cd /home/USER_ID/miniconda3/envs/general_tools/share/RepeatMasker/ # USER_ID is your user name i.e. akinya

    ./confiure # runs the configuration step

Set execution path of tfr, e.g. /home/USER_ID/miniconda3/envs/USER_ENV/bin/trf
Add search engine. Option 2 - RMBlast will be used
Set path where rmblastn and makeblastdb are found, e.g. /home/USER_ID/miniconda3/envs/USER_ENV/bin
5. Done to exit

### RepeatMask & TPSI path 1

    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking
    BestAssembly=assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta
    OutDir=repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/SMARTdenovo/ncbi_edits_repmask
    sbatch $ProgDir/rep_modeling.sh $BestAssembly $OutDir
    sbatch $ProgDir/transposonPSI.sh $BestAssembly $OutDir

output location - /analysis/blast_homology/Organism/strain
change blast output to xyz.tsv and open in excel
wc -l command gives number of what ou're looking for

    less *_unmasked.fa.TPSI.allHits.chains.bestPerLocus.gff3 | grep 'TY1_Copia' | grep 'E=0.0' | wc -l

To grep 2 lines

    less repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_RepMod-families.stk | grep -E 'ID|hAT'

Running repeat modeler to find LTR structure
    screen -a
    srun --partition himem --time 0-06:00:00 --mem-per-cpu 40G --cpus-per-task 24 --pty bash

    RepeatModeler -pa 10 -ninja_dir /scratch/software/NINJA-0.97-cluster_only/NINJA -LTRStruct -database “$Strain”_RepMod

Onecodetofindthemall/README

Make files executable
    chmod u+x build_dictionary.pl
    chmod u+x one_code_to_find_them_all.pl

Usage
    ./build_dictionary.pl --rm filename [--unknown] [--fuzzy] > output

--rm filename
Indicates the code to run on RepeatMasker output file filename. If filename is a directory, all .out files inside this directory and all sub-directories (recursively) will be scanned.
--unknown
Indicates to the code that the RepeatMasker output file passed in the --rm option may contain transposons of "Unkown" class/family, which is the case in particular when a local library was used in running RepeatMasker. If this option is passed, these "Unknown" elements will be included
--fuzzy
Indicates to the code to be less stringent in the criterions used to match names between corresponding subparts. May be useful to reconstruct all internal-LTR pairs in the data, but with a much higher proportion of false positives

Usage:
    ./one_code_to_find_them_all.pl --rm filename --ltr file_dictionary [--length file_length] [--strict] [--choice] [--dry-run]  [--fasta file_fasta [--flanking X]] [--insert Y]

ltr.csv and .transposons.csv) report the copies found
Score                   the score reported by RepeatMasker in the input file.
%_Div                   the percentage of divergence of the copy from the reference, computed by averaging the fragments divergences
%_Del                   similar to %_Div for deletion
%_Ins                   similar to %_Div for insertions
Query                   scaffold on which the copy was found
Beg.                            start position, on the scaffold, of the copy
End                             end position, on the scaffold, of the copy
Length                  the element length in the genomic sequence
Sense                   + if the transposon is inserted in 5'->3', C for a 3'->5' insertion
Element                 the transposable element name
Family                  the transposable element family or class
Pos_Repeat_Beg  the start of the actual sequence relative to the reference
Pos_Repeat_End  the end of the actual sequence relative to the reference
Pos_Repeat_End  the count of bases in the reference located after the end of the matching actual sequence
ID                              the RepeatMasker fragment ID
Num_Assembled   the number of fragments merged or assembled into this copy
