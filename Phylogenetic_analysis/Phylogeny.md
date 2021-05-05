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

    printf "" > analysis/popgen/busco_phylogeny/single_hits_15xStraw.txt
    for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
    echo $Busco
    OutDir=analysis/popgen/busco_phylogeny/$Busco
    mkdir -p $OutDir
    for Fasta in $(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/ncbi_edits_repmask/busco_sordariomycetes_obd10/DSA14_003_contigs_unmasked/run_sordariomycetes_odb10/busco_sequences/single_copy_busco_sequences/$Busco*.fna); do
    Strain=$(echo $Fasta | rev | cut -f8 -d '/' | rev)
    Organism=$(echo $Fasta | rev | cut -f9 -d '/' | rev)
    FileName=$(basename $Fasta)
    cat $Fasta | sed "s/:.*.fasta:/:"$Organism"_"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
    done
    cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
    SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | wc -l)
    printf "$Busco\t$SingleBuscoNum\n" >> analysis/popgen/busco_phylogeny/single_hits_15xStraw.txt
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

Remove busco genes that give unexpected multiple hits in the previous analysis.
    rm path/to/*/*/*/single_copy_busco_sequences/

If all isolates have a single copy of a busco gene, move the appended fasta to a new folder

    OutDir=analysis/popgen/busco_phylogeny/alignments
    mkdir -p $OutDir
    OrganismNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | cut -f2 | sort -nr | head -n1)
    for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
    echo $Busco
    HitNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | grep "$Busco" | cut -f2)
    if [ $HitNum == $OrganismNum ]; then
      cp analysis/popgen/busco_phylogeny/$Busco/"$Busco"_appended.fasta $OutDir/.
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

    cd analysis/popgen/busco_phylogeny/alignments
    scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/phylogenetics
    python $scripts/calculate_nucleotide_diversity.py "*aligned.fasta"
