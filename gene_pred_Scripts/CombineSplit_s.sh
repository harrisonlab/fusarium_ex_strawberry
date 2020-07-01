#Output files: DSA1*_interpro.gff3 , DSA*_interproscan.tsv , DSA1*_interproscan.xml
#outdir gene_pred/interproscan/F.oxysporum_fsp_fragariae/DSA*

ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  for Proteins in $(ls gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA15_041/final/final_genes_appended_renamed.pep.fasta); do
    Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    echo $Strain
    InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
    $ProgDir/append_interpro.sh $Proteins $InterProRaw
  done

#tsv Output file gives is ordered like this:
#1-Protein Accession (e.g. P51587)
#Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
#Sequence Length (e.g. 3418)
#Analysis (e.g. Pfam / PRINTS / Gene3D)
#Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
#Signature Description (e.g. BRCA2 repeat profile)
#Start location
#Stop location
#Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
#Status - is the status of the match (T: true)
#Date - is the date of the run
#(InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
#(InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
#(GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
#15-(Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
