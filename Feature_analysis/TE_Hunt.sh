###
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

# Set execution path of tfr, e.g. /home/USER_ID/miniconda3/envs/USER_ENV/bin/trf
# Add search engine. Option 2 - RMBlast will be used
# Set path where rmblastn and makeblastdb are found, e.g. /home/USER_ID/miniconda3/envs/USER_ENV/bin
# 5. Done to exit

### RepeatMask & TPSI path 1

    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking
    BestAssembly=assembly/SMARTdenovo/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta
    OutDir=repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/SMARTdenovo/ncbi_edits_repmask
    sbatch $ProgDir/rep_modeling.sh $BestAssembly $OutDir
    sbatch $ProgDir/transposonPSI.sh $BestAssembly $OutDir

# output location - /analysis/blast_homology/Organism/strain
# change blast output to xyz.tsv and open in excel
# wc -l command gives number of what ou're looking for

less *_unmasked.fa.TPSI.allHits.chains.bestPerLocus.gff3 | grep 'TY1_Copia' | grep 'E=0.0' | wc -l

# to grep 2 lines
less repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_RepMod-families.stk | grep -E 'ID|hAT'

# Running repeat modeler to find LTR structure
screen -a
srun --partition himem --time 0-06:00:00 --mem-per-cpu 40G --cpus-per-task 24 --pty bash

RepeatModeler -pa 10 -ninja_dir /scratch/software/NINJA-0.97-cluster_only/NINJA -LTRStruct -database “$Strain”_RepMod
