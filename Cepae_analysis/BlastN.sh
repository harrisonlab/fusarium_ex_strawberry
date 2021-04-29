# Blast pipe search
# Run in conda env with perly (Repenv)
# Query=path/to/query/fasta - path to six genes oldhome/groups/harrisonlab/project_files/fusarium/analysis/blast_homology/six_genes/six-appended_parsed.fa
# for $Assembly Use files with nucleotides
# Assembly can be final_genes_appended_renamed.fasta or DSA14_003_contigs_unmasked.fasta
# Do both
# Had to cp and edit blast_pipe - slight error in directories

for Assembly in $(ls repeat_masked/F.oxysporum_fsp_fragariae/15-074/15-074_contigs_unmasked.fa); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
  echo "$Organism - $Strain"
  Query=../oldhome/groups/harrisonlab/project_files/fusarium/analysis/blast_homology/six_genes/six-appended_parsed.fa
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis
  OutDir=analysis/blast_homology/$Organism/$Strain
  sbatch $ProgDir/blast_pipe.sh $Query dna $Assembly $OutDir
done

# To turn dna from DNA to proteins use this but use full paths
# sbatch $ProgDir/blast_pipe.sh $Query DNA $Assembly
# output location - /analysis/blast_homology/Organism/strain
