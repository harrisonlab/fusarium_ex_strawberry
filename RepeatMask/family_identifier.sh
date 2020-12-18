#Use to find families in Rep modeler output
#Identifying TE families on cand LS/PS contig

Sort no. TEs by family and function
    for TEfamily in $(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_RepMod-families.stk); do
    	Organism=$(echo $Antismash | rev | cut -f4 -d '/' | rev)
    	Strain=$(echo $Antismash | rev | cut -f3 -d '/' | rev)
    	echo "$Organism - $Strain"
    	cat $TEfamily | grep 'DNA_Transposon' | sort | uniq -c
    done

    # grep -w -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139' | grep -e 'A8' -e 'PG' -e 'A23'

    for TEfamily in $(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_RepMod-families.stk); do
      Organism=$(echo $Antismash | rev | cut -f4 -d '/' | rev)
      Strain=$(echo $Antismash | rev | cut -f3 -d '/' | rev)
      echo "$Organism - $Strain"
      cat $TEfamily | grep -E 'rnd-1\|DNA_Transposon\|contig_12' | sort | uniq -c
    done

>>>

#Finding things with grep
#grep 'extra\|value\|service' sample.txt
#grep -E ‘extra|value|service’ sample.txt
#egrep ‘extra|value|service’ sample.txt
#grep -e extra -e value -e service sample.txt

# Ran repeatmasker on contig 12 fasta
# repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/contig_12/DSA14_003_RepMod-families.stk
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Repeat_masking
  BestAssembly=assembly/flye/F.oxysporum_fsp_fragariae/DSA14_003/pilon/contig_12.fasta
  OutDir=repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/contig_12
  sbatch $ProgDir/rep_modeling.sh $BestAssembly $OutDir
  sbatch $ProgDir/transposonPSI.sh $BestAssembly $OutDir

# pulling families just on contig 12 repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/contig_12/DSA14_003_RepMod-families.stk
# search for diff values: Retrotransposon , Retrotransposed_Element , DNA_Transposon

for TEfamily in $(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/contig_12/DSA14_003_RepMod-families.stk); do
  Organism=$(echo $Antismash | rev | cut -f5 -d '/' | rev)
  Strain=$(echo $Antismash | rev | cut -f4 -d '/' | rev)
  echo "$Organism - $Strain"
  cat $TEfamily | grep 'Retrotransposed_Element' | sort | uniq -c
done
