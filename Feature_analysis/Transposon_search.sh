# Finding transposon sequences for qPCR primers

# 1) make a text file for the genes you think are TE's from interproscan data
  nano DNA_H_cand.txt # or transposon_cand.txt

# 2) copy gene names listed with this command into "DNA_H_cand.txt" or "transposon_cand.txt"
# Use cut -f1 or * DSA14_003_eff_ortho_genes.fasta_homologs.csv to excise and view column data
  less gene_pred/interproscan/F.oxysporum_fsp_fragariae/DSA14_003/DSA14_003_interproscan.tsv | grep 'DNA helicase Pif1-like'
g10299.t1
g13188.t1
g13188.t1
g13188.t1
g15913.t1
g15705.t1
g16426.t2
g16584.t1
g16535.t1
g16906.t1
g16906.t1
g17176.t1
g17610.t1
g17702.t1
g17762.t1
g18007.t1

# 3) Extract genesequences from ILLUMINA FRAGARIAE genome
faidx -d '|' final_genes_appended_renamed.cdna.fasta $(tr '\n' ' ' < DNA_H_cand.txt ) > DNA_H_Pif1.fasta

# 4) Now you have gene sequence find intron and exon regions in gene in DSA14_003_interproscan

# look for reverse complement too -  cat dummyfile.fa | tr -d "\n" > rev_comp.fa

cat repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_RepMod-families.stk | grep 'DNA' | cut -f1 | cut -f2 | sort | uniq > DNA_Tsons_RepMod.txt
# not tsv file so it won't remove things in "column1 & 2"
# count no of each element
cat DNA_Tsons_RepMod.txt | grep 'Helitron-1' | wc -l

####
# Manual blast for imapala sequences
# found sequences on NCBI

# Blast pipe search
# Run in conda env with perly (Repenv)
# for $Assembly Use files with nucleotides
# Assembly can be final_genes_appended_renamed.fasta or DSA14_003_contigs_unmasked.fasta
# Do both
# Had to cp and edit blast_pipe - slight error in directories

    for Assembly in $(ls assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta); do
      Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
      Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
      echo "$Organism - $Strain"
      Query=Impala_seq.fa
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
      sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly
    done

# output location - /analysis/blast_homology/Organism/strain
