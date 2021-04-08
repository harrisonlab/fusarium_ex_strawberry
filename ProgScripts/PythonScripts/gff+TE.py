#!/usr/bin/python

'''
Script to take coordinates of genes and predicted TEs. Genes & TEs are identified and labelled in new gff3 file.
'''

import sys,argparse
from collections import defaultdict
from sets import Set


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument('--inp_gff',required=True,type=str,help='input genes file')
ap.add_argument('--te_gff', required =True, type=str, default = False, help = 'input TE file')
ap.add_argument('--Output', required=True,type=str,help='output file required')
conf = ap.parse_args() #sys.argv


#-----------------------------------------------------
# Step 2
# Order gff features by input contigs and gene start
#-----------------------------------------------------

open(conf.inp_gff)
conf.inp_gff = f
inp_lines = f.readlines()
contig_list=[] # This will be the final_genes_appended.gff which is finished annotation of genes in genome
gene_hash = {} # Will giive the gene names i.e. g1.t1 a value in the contig it is located in to sort
for line in inp_lines:
    line = line.strip("\n")
    if line.startswith("contig"):
        print line
        continue
    split_line = line.split()
    if split_line[2] == 'gene':
        contig = split_line[0]
        gene_start = split_line[3]
        gene_end = split_line[4]
        gene_id = split_line[8]
        if contig not in contig_list:
            contig_list.append(contig)
        gene_start_dict[contig].append(int(gene_start))
        key = "_".join([contig, gene_start])
    features_dict[key].append(line)
    if not 'gene':
        continue # only interested in genes that overlap TEs
    def overlap(gene_start,gene_end,TE_start,TE_end):
         if gene_start<=TE_end and TE_start<=gene_end:
            return true
    else:
            return false if not

open(conf.te_gff)
conf.te_gff = TE
INP_ROWS = TE.readlines()
CONTIG_LIST = []
transp_hash = {}
for each line in INP_ROWS:
    ROW = ROW.strip("\n")
    if LINE.startswith("contig"):
        print ROW
        continue
    split_ROW = ROW.split()
    if split_ROW[2] == 'translated_nucleotide_match':
        CONTIG = split_ROW[0]
        TE_start = split_ROW3]
        TE_end = split_ROW[4]
        TE_ID = split_ROW[8]

if Output:
    z = open(Output, "w")
for each contig in gene_hash:
    for each gene in list:
        for each transp item in transp_hash for this contig_list
            if overlap(gene_start,gene_end,TE_start,TE_end):
                z.write(contig,gene_id,gene_start,gene_end,TE_start,TE_end,TE_id)
            z.close()
