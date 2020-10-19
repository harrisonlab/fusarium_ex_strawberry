# Run in Salmon env
# This script simply loops through each sample and invokes salmon
# Fus2_CzapekDox, Fus2_GlucosePeptone, Fus2_PDA and Fus2_PDB are RNAseq data of infected onions
# -i - where to find the index
# -l A - automatically determines the library type of the sequencing reads
# -1 and -2 - where to find the left and right reads for this sample
# -p 8 - threads count
# -o output dir

#Make index first - done for both isolates
salmon index -t codingquary/F.oxysporum_fsp_fragariae/DSA14_003/final/final_genes_appended_renamed.cdna.fasta -i salmon/DSA14_003_index

# srun --partition long --time 0-06:00:00 --mem 40G --cpus-per-task 24 --pty bash
#!/bin/bash
for fn in ../oldhome/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/ ;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i gene_pred/salmon/DSA14_003_index -l A \
         -1 ${fn}/F/6_S2_L001_R1_001_trim.fq.gz \
         -2 ${fn}/R/6_S2_L001_R2_001_trim.fq.gz \
         -p 8 --validateMappings -o gene_pred/salmon/quants/${samp}_quant
done

# gave same error 3 TrimReads#
# Detected a *potential* strand bias > 1% in an unstranded protocol check the file: gene_pred/salmon/quants/Fus2_CzapekDox_quant/lib_format_counts.json for details
# > ./bin/salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fq -2 reads2.fq --validateMappings -o transcripts_quant
