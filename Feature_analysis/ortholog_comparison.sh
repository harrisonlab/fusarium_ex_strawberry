# Use to extract effectors from a gene.fasta file using names in txt file
# Create text file with known gene names of effectors/mimps you want to remove
# i.e  input gene names in Fof14_genes.txt
# Command uses gene name to copy across fasta seq to Fof genes

# Extract genes from cepae genome NOT FRAGARIAE

# faidx -d '|' final_genes_appended_renamed.cdna.fasta $(tr '\n' ' ' < Fof14_genes.txt ) > Fof14_genes.fasta - example command to test

# /projects/fusarium_ex_strawberry/gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA14_003/final/
# Why are you doing this?
# To compare candidate effectors in Fo cepae (which has RNA seq data) against the predicted genes in Fof using cepae RNA seq data

faidx -d '|' final_genes_combined.cdna.fasta $(tr '\n' ' ' < FocFus2_genes.txt ) > eff_ortho_genes.fasta

# Now that you have the cand effector genes, contrast against Fo_fragariae_14_003 long read seq genome
# Run in conda env with perly (Repenv)
# for $Assembly Use files with nucleotides

for Assembly in $(ls assembly/miniasm/F.oxysporum_fsp_fragariae/DSA14_003/pilon/pilon_10_renamed.fasta); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  Query=../F.oxysporum_fsp_cepae/Fus2_canu_new/final/eff_ortho_genes.fasta
  OutDir=assembly/miniasm/$Organism/$Strain/Orthology
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
  sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly $OutDir
done

# Use cut -f1 or * DSA14_003_eff_ortho_genes.fasta_homologs.csv to excise and view column data
