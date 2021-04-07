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
ap.add_argument('--inp_gff',required=True,type=str,help='input gff file')
ap.add_argument('--conversion_log', required = False, type=str, default = False, help = 'File name to store converted gene if required')
conf = ap.parse_args() #sys.argv

if conf.conversion_log:
    logfile = conf.conversion_log

with open(conf.inp_gff) as f:
    inp_lines = f.readlines()

#-----------------------------------------------------
# Step 2
# Order gff features by input contigs and gene start
#-----------------------------------------------------


def function has_overlap(start1,end1,start2,end2):
    # This function checks to determine if two ranges overlap.
    # The ranges are defined by two co-ordinates - start and end
    # note this function assumes that start <= end and that the contig ids match
    # Why?
    # The function returns true if?
    # and False if not.
     if start1<=end2 and start2<=end1 (or something like that)
        return true
    else
        return false if not
open gene gff file # This will be the final_genes_appended.gff which is finished annotation of genes in genome
gene_hash = {} # Will giive the gene names i.e. g1.t1 a value in the contig it is located in to sort
for each line:
    split line into list of tokens
    if not 'gene' item line continue to next line # only interested in genes that overlap TEs
    if contig_list not yet in genehash create blank list under this key
    contig_list=[] # create blank list using square brackets - positioning may be incorrect
    store [geneid,genestart,geneend] in genehash under contig_list
open transp gff file
transp_hash = {}
for each line:
    split line into list of tokens
    if not relevant line continue to next line # Excises
    if transpid not yet in transphash create blank list under this key
    store [transpid,start,end] in transphash under contigid # if the transcript ID matches an entry in the hash then print the previous gene as the parent feature (in gff3 format)
    transpid=[] # - positioning may be incorrect
open output file
for each contigid in genehash:
    for each gene item in list:
        for each transp item in transphash for this contig_list
            if has_overlap(start1,end1,start2,end2):
                write contigid,geneid,transpid etc # print array of line add 10th column for TE ID and hit no from bestPerLocus file
