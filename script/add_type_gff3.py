#!/usr/bin/env python
# coding: utf-8

import re
import gzip
import argparse
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Adding amino acid type to gff3 file.", add_help=False, 
                                     epilog='date:2025/09/16 author:guisen chen email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('gff3', metavar='gff3', help='A input file of gff3 format.')
    required.add_argument('pep', metavar='pep', help='A input file of fasta format, including protein sequences.')
    required.add_argument('-o', '--output', metavar='str', help='A output file of gff3 format.', required=True)
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00', help="Show program's version number and exit.")
    args = parser.parse_args()
    gff3 = args.gff3
    pep = args.pep
    outfile = args.output
    
mRNA2Type = {}
statType = {"complete":0, "3prime_partial":0, "5prime_partial":0, "internal":0}
for record in SeqIO.parse(pep, 'fasta'):
    if record.seq[0].upper() == "M" and record.seq[-1].upper() == "*":
        mRNA2Type.setdefault(record.id, 'complete')
        statType['complete'] += 1
    if record.seq[0].upper() == "M" and record.seq[-1].upper() != "*":
        mRNA2Type.setdefault(record.id, '3prime_partial')
        statType['3prime_partial'] += 1
    if record.seq[0].upper() != "M" and record.seq[-1].upper() == "*":
        mRNA2Type.setdefault(record.id, '5prime_partial')
        statType['5prime_partial'] += 1
    if record.seq[0].upper() != "M" and record.seq[-1].upper() != "*":
        mRNA2Type.setdefault(record.id, 'internal')
        statType['internal'] += 1

print("Type\tGene Number")
for i in statType:
    print(i,statType[i], sep='\t')
    
out = open(outfile, 'w')
with open(gff3, 'r') as f:
    for l in f:
        if l.strip() != "" and not l.startswith('#'):
            if l.strip().split('\t')[2] == 'mRNA':
                ID = re.search("ID=(.*?);",l.strip().split('\t')[8]).group(1)
                if l.strip()[-1] == ";":
                    print(l.strip() + f"Type={mRNA2Type[ID]};", file=out)
                else:
                    print(l.strip() + f";Type={mRNA2Type[ID]};", file=out)
            else:
                print(l.strip(), file=out)
        else:
            print(l.strip(), file=out)
out.close()
