#!/usr/bin/python

'''
Script to take coordinates of genes and predicted TEs. Genes & TEs are identified and labelled in new gff3 file.
'''
# $ProgDir/TE_to_gff.py --inp_gff $Gff --te_gff $TE_gff > $Gff+TE
import sys,argparse
from collections import defaultdict
from sets import Set


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument('--inp_gff',required=True,type=str,help='input gff file')
ap.add_argument('--te_gff', required =True, type=str, default = False, help = 'TE gff file')
conf = ap.parse_args() #sys.argv

open(conf.inp_gff)
conf.inp_gff = f
inp_lines = f.readlines()

with open(conf.te_gff) as g:
    te_lines= g.readlines()


#-----------------------------------------------------
# Step 2
# Order gff features by input contigs and gene start
#-----------------------------------------------------

CONTIG_LIST = []
TE_start_dict = defaultdict(list)
TE_end_dict = defaultdict(list)
Tson_dict = defaultdict(list)
contig_list = []
gene_start_dict = defaultdict(list)
gene_end_dict = defaultdict(list)
features_dict = defaultdict(list)

for line in inp_lines:
    line = line.strip("\n")
    if line.startswith("contig"):
        print line
        continue
    split_line = line.split()
    # Col9 = split_line[8]
    if split_line[2] == 'gene':
        contig = split_line[0]
        gene_start = split_line[3]
        gene_end = split_line[4]
        gene_ID = split_line[8]
        if contig not in contig_list:
            contig_list.append(contig)
        gene_start_dict[contig].append(int(gene_start))
        gene_end_dict[contig].append(int(gene_end))
        key = "_".join([contig, gene_start, gene_end])
    features_dict[key].append(line)

gff_lines = []
# for contig in contig_list:
for contig in sorted(contig_list, key = lambda x: (int(x.split('_')[1]))):
    gene_start_list = gene_start_dict[contig]
    for gene_start in sorted(gene_start_list):
        # print "\t".join([contig, str(gene_start)])
        key = "_".join([contig, str(gene_start, gene_end)])
        gff_lines.extend(features_dict[key])
        # print("\n".join(features_dict[key]))

for ROW in te_lines:
    ROW = ROW.strip("\n")
    if ROW.startswith("contig"):
        print ROW
        continue
    split_ROW = ROW.split()
    # Col9 = split_line[8]
    if split_ROW[2] == 'translated_nucleotide_match':
        CONTIG = split_ROW[0]
        TE_start = split_ROW[3]
        TE_end = split_ROW[4]
        if CONTIG not in CONTIG_LIST:
            CONTIG_LIST.append(CONTIG)
        TE_start_dict[CONTIG].append(int(TE_start))
        TE_end_dict[CONTIG].append(int(TE_end))
        key = "_".join([CONTIG, TE_start, TE_end])
    Tson_dict[key].append(ROW)

TE_lines = []
# for CONTIG in contig_list:
for CONTIG in sorted(CONTIG_LIST, key = lambda x: (int(x.split('_')[1]))):
    TE_start_list = TE_start_dict[CONTIG]
    for TE_start in sorted(TE_start_list):
        # print "\t".join([contig, str(gene_start)])
        key = "_".join([CONTIG, str(TE_start, TE_end)])
        TE_lines.extend(Tson_dict[key])
        # print("\n".join(features_dict[key]))

#-----------------------------------------------------
# Step 3
# MERGE OVERLAPS
#-----------------------------------------------------

for each line in TE_lines:
    def function has_overlap(gene_start,gene_end,TE_start,TE_end):
         if gene_start<=TE_end and TE_start<=gene_end
            return true
        else
            return false if not

#open output file
#for each contigid in genehash:
#    for each gene item in list:
#        for each transp item in transphash for this contig_list
#            if has_overlap(start1,end1,start2,end2):
#                write contigid,geneid,transpid etc # print array of line add 10th column for TE ID and hit no from bestPerLocus file
