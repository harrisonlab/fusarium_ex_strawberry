# Phylogeny analysis
Based on BUSCO genes
If you have yet to run it already, run Busco
Make a BUSCOenv to correctly analyse genome

    conda create -n BUSCOenv
    conda install -c bioconda busco

Running script
    for Assembly in $(ls path/to/genome); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
      Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
      echo "$Organism - $Strain"
      ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/quality_check
      BuscoDB=$(ls -d /projects/dbBusco/sordariomycetes_odb10)
      OutDir=$(dirname $Assembly)/busco_sordariomycetes_obd10
      sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir
    done

Create a list of all BUSCO IDs

    OutDir=analysis/popgen/busco_phylogeny
    mkdir -p $OutDir
    BuscoDb="sordariomycetes_odb10"
    ls -1 /projects/dbBusco/$BuscoDb/hmms/*hmm | rev | cut -f1 -d '/' | rev | sed -e 's/.hmm//' > $OutDir/all_buscos_"$BuscoDb".txt

For each busco gene create a folder and move all single copy busco hits from each assembly to the folder. Then create a fasta file containing all the aligned reads for each busco gene for alignment later.
repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/ncbi_edits_repmask/busco_sordariomycetes_obd10/DSA14_003_contigs_unmasked/run_sordariomycetes_odb10/busco_sequences/single_copy_busco_sequences/100017at147550.fna
assembly/F.oxysporum_fsp_fragariae/15-074/busco_sordariomycetes_obd10/15-074_contigs_unmasked/run_sordariomycetes_odb10/busco_sequences/single_copy_busco_sequences/100017at147550.fna
repeat_masked/busco_sordariomycetes_obd10/4287_chromosomal_contigs_unmasked/run_sordariomycetes_odb10/busco_sequences/single_copy_busco_sequences/
repeat_masked/busco_sordariomycetes_obd10/Fus2_canu_new_contigs_unmasked/run_sordariomycetes_odb10/busco_sequences/single_copy_busco_sequences/

    printf "" > analysis/popgen/busco_phylogeny/single_hits.txt
    for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
    echo $Busco
    OutDir=analysis/popgen/busco_phylogeny/$Busco
    mkdir -p $OutDir
    for Fasta in $(ls NextGenSeq/Fo_lycopersici/repeat_masked/busco_sordariomycetes_obd10/4287_chromosomal_contigs_unmasked/run_sordariomycetes_odb10/busco_sequences/single_copy_busco_sequences/$Busco.fna); do
    Strain=4827
    Organism=F.oxysporum_fsp_lycopersici
    FileName=$(basename $Fasta)
    cat $Fasta | sed "s/:.*.fasta:/:"$Organism"_"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
    done
    cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
    SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | wc -l)
    printf "$Busco\t$SingleBuscoNum\n" >> analysis/popgen/busco_phylogeny/single_hits.txt
    done

Repeat for next gen seq
    for Busco in $(cat NextGenSeq/analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
    echo $Busco
    OutDir=analysis/popgen/busco_phylogeny/$Busco
    mkdir -p $OutDir
    for Fasta in $(ls NextGenSeq/assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/pilon/busco_sordariomycetes_obd10/rounds/pilon_10/run_sordariomycetes_odb10/busco_sequences/single_copy_busco_sequences/$Busco*.fna); do
    Strain=$(echo $Fasta | rev | cut -f9 -d '/' | rev)
    Organism=$(echo $Fasta | rev | cut -f10 -d '/' | rev)
    FileName=$(basename $Fasta)
    cat $Fasta | sed "s/:.*.fasta:/:"$Organism"_"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco"_NGS.fasta
    done
    cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended_NGS.fasta
    SingleBuscoNum=$(cat $OutDir/"$Busco"_appended_NGS.fasta | grep '>' | wc -l)
    printf "$Busco\t$SingleBuscoNum\n" >> NextGenSeq/analysis/popgen/busco_phylogeny/single_hits_NGS.txt
    done

Check for multiple hits
    less analysis/popgen/busco_phylogeny/single_hits.txt | sort -k2 -n
 Antonio's fasta file, N.ditissima_Ag02_EOG0933000J.fasta , looks like this
  #>EOG0933000J:Ag02_contigs_unmasked.fa:NODE_29_length_270755_cov_66.213485:16680-31581
  atggccatggttgacgtctcacggcagaggaggtcgcttctccaggacgcgactgttttggaacacctgccatcagaattactggccattatca...
Need to rename headers in fasta files so it has a downstream effect
    cd path/to/FASTAfiles
    sed -i 's/Organism_strain/whats is present in fastaheader apart from busconumber/g' *
    # e.g.
    sed -i 's/$Busco:F.oxysporum_fsp_fragariae_15-074:NODE_/NODE_/g' F.oxysporum_fsp_fragariae_15-074_*.fasta
test it first
    sed -i 's/$Busco:F.oxysporum_fsp_fragariae_15-074:NODE_/NODE_/g' F.oxysporum_fsp_fragariae_15-074_100017at147550.fasta

Remove busco genes that give unexpected multiple hits in the previous analysis.
    rm path/to/*/*/*/single_copy_busco_sequences/

If all isolates have a single copy of a busco gene, move the appended fasta to a new folder (edited for Fo_cepae and lycopersici directories)

    OutDir=../../analysis/popgen/busco_phylogeny/alignments
    mkdir -p $OutDir
    OrganismNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | cut -f2 | sort -nr | head -n1)
    for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
    echo $Busco
    HitNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | grep "$Busco" | cut -f2)
    if [ $HitNum == $OrganismNum ]; then
      cp ../../analysis/popgen/busco_phylogeny/$Busco/"$Busco"_appended.fasta $OutDir/.
    fi
    done

Submit alignment for single copy busco genes with a hit in each organism
Activtate any conda env with Python V3 +
    conda install -c bioconda mafft
Raw usage
    mafft-linsi input > output
    #Convert to single line FASTA for easy parsing
    awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $output >temp && mv temp $output

Command line usage
    AlignDir=analysis/popgen/busco_phylogeny/alignments
    CurDir=$PWD
    cd $AlignDir
    ProgDir=/home/akinya/git_repos/fusarium_ex_strawberry/ProgScripts/Genome_alignment
    sbatch $ProgDir/Mafft.sh
    cd $CurDir

For closely related organisms (same species etc.): identify genes with high nucleotide diversity (Pi) and average number of pairwise differences, medium number of segregating sites (avoid alignments with low homology and lots of phylogenetically uninformative singletons). For analyses involving cross-species comparisons involving highly diverged sequences with high nucleotide diversity (e.g. 0.1<Pi<0.4), looking for genes with the lowest number of segregating sites.

Run in env with perl
    conda activate Repenv
Need to install dendropy
    python3 -m pip install -U dendropy
Run analysis
    cd analysis/popgen/busco_phylogeny/alignments
    scripts=/home/gomeza/git_repos/scripts/bioinformatics_tools/Phylogenetics
    python $scripts/calculate_nucleotide_diversity.py "*aligned.fasta"
    #### you made by accident 100017at147550_appended_aligned.fasta 100017at147550_appended_aligned_aligned.fasta 100017at147550_appended_aligned_aligned_aligned.fasta & 100017at147550_appended_aligned_aligned_aligned_aligned_aligned.fasta

Got message but not an error:
Traceback (most recent call last):
  File "/home/gomeza/git_repos/scripts/bioinformatics_tools/Phylogenetics//calculate_nucleotide_diversity.py", line 26, in <module>
    pi = dendropy.calculate.popgenstat.nucleotide_diversity(seqs, ignore_uncertain=True)
  File "/home/akinya/miniconda3/envs/Repenv/lib/python3.8/site-packages/dendropy/calculate/popgenstat.py", line 184, in nucleotide_diversity
    return _nucleotide_diversity(char_matrix.sequences(), char_matrix.default_state_alphabet, ignore_uncertain)
  File "/home/akinya/miniconda3/envs/Repenv/lib/python3.8/site-packages/dendropy/calculate/popgenstat.py", line 90, in _nucleotide_diversity
    return _count_differences(char_sequences, state_alphabet, ignore_uncertain)[1]
  File "/home/akinya/miniconda3/envs/Repenv/lib/python3.8/site-packages/dendropy/calculate/popgenstat.py", line 45, in _count_differences
    raise Exception("sequences of unequal length")
Exception: sequences of unequal length

## Trimming sequence alignments using Trim-Al
Trimming sequence alignments using Trim-Al. Note - automated1 mode is optimised for ML tree reconstruction
install trimal
    conda install -c bioconda trimal

Run script
    OutDir=analysis/popgen/busco_phylogeny/trimmed_alignments
    mkdir -p $OutDir
    for Alignment in $(ls analysis/popgen/busco_phylogeny/alignments/*_appended_aligned.fasta); do
      TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
      echo $Alignment
      trimal -in $Alignment -out $OutDir/$TrimmedName -keepheader -automated1
    done

Randomized Axelerated Maximum Likelihood
    #Edit header name keeping BUSCO name and isolate name
    cd analysis/popgen/busco_phylogeny/trimmed_alignments
    sed -i 's/_contigs_.*//g' *_appended_aligned_trimmed.fasta                                                                              
    sed -i 's/:LD.*//g' *_appended_aligned_trimmed.fasta
    sed -i 's/N.ditissima_//g' *_appended_aligned_trimmed.fasta

    screen -a
    #This will run in the short queue for testing
    for Alignment in $(ls analysis/popgen/busco_phylogeny/trimmed_alignments/*aligned_trimmed.fasta); do
        Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
        while [ $Jobs -gt 5 ]; do
        sleep 5m
        printf "."
        Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
        done
        printf "\n"
        Prefix=$(basename $Alignment | cut -f1 -d '_')
        OutDir=analysis/popgen/busco_phylogeny/RAxML/$Prefix
        ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Phylogenetics
        squeue $ProgDir/RAxML.sh $Alignment $Prefix $OutDir
    done
