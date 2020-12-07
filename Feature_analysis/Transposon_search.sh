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
