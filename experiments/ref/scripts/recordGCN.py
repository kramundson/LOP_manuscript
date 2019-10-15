#!/usr/env/python
# RecordGCN.py
# Kirk Amundson
# 25 October 2018

"""
Consume a multifasta file. Record GC and N content in nonoverlapping bins of a
size that is both fixed and specified by the user.

USAGE: python script.py -s binsize -f fastafile -o outfile
"""

import argparse, re
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Find GC and N content in windows')
parser.add_argument('-s', '--step', dest='step' , type=int, help='Nonoverlapping step size')
parser.add_argument('-f', '--file', dest='fasta', type=str, help='Input file')
parser.add_argument('-o', '--out' , dest='ofh'  , type=str, help='Output file')
args = parser.parse_args()

step = args.step
o = open(args.ofh, 'w')
o.write('\t'.join(["chrom","start","end","binsize","GC_content","N_content\n"]))

with open(args.fasta, 'r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        print("Now working on sequence:", record.name)
        start = 0
        while start < len(record.seq):
            test = record.seq[start:start+step].upper()
            out = [str(i) for i in [record.name,
                                    start,
                                    start+step,
                                    step,
                                    len(re.findall("[GC]", str(test))) / len(test),
                                    len(re.findall("[N]", str(test))) / len(test)]
                                   ]
            o.write('\t'.join(out)+'\n')
            start += step
o.close()
