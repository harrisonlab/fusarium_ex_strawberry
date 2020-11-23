#!/usr/bin/python

import os
import sys
import optparse
import os.path
import re

'''
Add ID and Name fields to repeat masker gff3

Usage: fix_rm_gff3.py -i <repeat_masker.gff> -o <corrected_repeat_masker.gff>

Written by Rob Vickerstaff, East Malling Research

Modified by Andrew Armitage, East Malling Research
'''

from optparse import OptionParser

parser=OptionParser()
parser.add_option("-i",   dest="inp",   default="NONE",   help="-i <repeat_masker.gff>"   )
parser.add_option("-o",   dest="out",   default="NONE",   help="-o <corrected_repeat_masker.gff>"   )
(options, args) = parser.parse_args()

inp	= options.inp
out = options.out

#inp='../../zipper_003/scr/ref2/020_nocontam/atlan20131014_030_masked.gff3'
#out='../../zipper_003/scr/ref2/020_nocontam/atlan20131014_030_masked_fixed.gff3'

f = open(inp)
fout = open(out,'wb')

ct = 0
for line in f:
    if line.startswith('#'):
        fout.write(line)
        continue
      
    ct += 1
    line = line.strip()
    tok = line.split('\t')[8]
    name = re.split(':|\"',tok)[2]
    line += ';ID=rep%d;Name=%s\n'%(ct,name)
    fout.write(line)

fout.close()
f.close()
