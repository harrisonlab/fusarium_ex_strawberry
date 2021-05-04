Fusarium-Oxy-fragariae
Adapted scripts from Andy's Verticillium nonalfalfae ex. hop pipeline

These scripts allow steps from the pipeline to be followed in the NewSlurm as opposed to the OldGridengine format that was previously used.

These scripts follow from the quast assembly step as this is where the edits in the pipeline had to be made.
The scripts are run in this order:
1. spades_step1.sh
2. quast_step2.sh
3. result_quast_step3.sh
4. remove-contaminants4.sh - This step should be run in an Anaconda shell with the latest version of Python (v3.8) the python file to run

The sequence reads and contigs need to be submitted to NCBI as Bioproject and a Biosample submission of genomes. After genome submission 
you are sent a "contamination.txt" to correct the contigs. They may send you a corrected .fasta file also

The steps continue as follows:
5. ncbi_rep_step.sh
6. ncbi_contig_fix_step6.sh
7. contig_quast7.sh
8. repeatmask_step8.sh

Still editing script and trouble shooting -- To be continued
