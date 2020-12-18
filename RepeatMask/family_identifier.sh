#Use to find families in Rep modeler output

Sort no. TEs by family and function
    for TEfamily in $(ls repeat_masked/F.oxysporum_fsp_fragariae/DSA14_003/flye/ncbi_edits_repmask/DSA14_003_RepMod-families.stk); do
    	Organism=$(echo $Antismash | rev | cut -f4 -d '/' | rev)
    	Strain=$(echo $Antismash | rev | cut -f3 -d '/' | rev)
    	echo "$Organism - $Strain"
    	cat $TEfamily | grep 'DNA_Transposon' | sort | uniq -c
    done

    # grep -w -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139' | grep -e 'A8' -e 'PG' -e 'A23'
