#!/usr/bin/env python
# coding: utf-8

import sys
import gzip
import argparse
from Bio import SeqIO

def main(infile, outfile):
    if outfile == None:
        out = sys.stdout
    elif outfile.endswith('.gz'):
        out = gzip.open(outfile, 'rt')
    else:
        out = open(outfile, 'w')

    for record in SeqIO.parse(infile, 'fasta'):
        print(f'>{record.id}', file=out)
        seqence = ''
        for i in record.seq:
            if i.islower():
                seqence += 'N'
            else:
                seqence += i
        print(seqence, file=out)
    out.close()
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Soft-masked (dna_sm) genome convert to masked(dna_rm) genome.',  add_help=False,
                                     epilog='Date:2024/03/14 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i',  '--input',  metavar='str',  help='Input file of fasta format',  required=True)
    optional.add_argument('-o',  '--output',  metavar='str',  help='Output file of fasta format.')
    optional.add_argument('-h',  '--help',  action='help',  help="Show program's help message and exit.")
    optional.add_argument('-v',  '--version',  action='version',  version='v1.00',  help="Show program's version number and exit.")
    args = parser.parse_args()
    main(infile=args.input, outfile=args.output)
