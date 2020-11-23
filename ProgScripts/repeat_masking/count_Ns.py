#!/usr/bin/python


'''
This toolcounts the number of N's present in a fasta file and returns the number
of N's in total and in each accession
'''

import sys,argparse
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument('--inp_fasta',required=True,type=str,help='.fasta file containing query contigs')
ap.add_argument('--out_txt',required=True,type=str,help='output text file')

conf = ap.parse_args()

out_file = open(conf.out_txt, "w")

for record in SeqIO.parse(conf.inp_fasta, "fasta"):
    N_count = record.seq.count("N") + record.seq.count("n")
    seq_length = len(record)
    perc_N = float(N_count*100.0/seq_length)
    perc_N = "%.2f" % round(perc_N, 2)
    out_file.write("\t".join([str(record.id), str(seq_length), str(N_count), str(perc_N)]) + "\n")
out_file.close()
