# Main phase Orthofinder
# Run in conda env (Repenv)

for IN_DIR in $(ls -d $WorkDir/formatted) ; do
sbatch /home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_analysis/orthofinder.sh $IN_DIR $IsolateAbrv
done
