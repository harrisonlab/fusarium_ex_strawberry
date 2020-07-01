#DNA was successfully split but interpro script didn't run for 33/35 of the split DNA
#Had to run split DNA in sets of 6
#Putting data in queue in short partition forced an error
#Everthing in queue would run for 2 seconds and crash

for file in $(ls gene_pred/interproscan/F.oxysporum_fsp_fragariae/DSA14_003/*_split_11*); do
  sbatch /home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation/run_interproscan.sh $file
  done
