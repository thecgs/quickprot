#!/usr/bin/env python
# coding: utf-8

import re
import sys
import argparse
from natsort import natsorted

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sort gff3 file.',  add_help=False,
                                     epilog='date:2024/11/18 author:guisen chen email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('input', metavar='str',
                          help='A input file of gff3 format.')
    optional.add_argument('-s', '--software', help='Change the source of the second line of gff3.', default=None)
    optional.add_argument('-o', '--output', metavar='str', default=None,  
                          help=f'A output file of gff3 format. defualt=None')
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00',  
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    inputfile = args.input
    outputfile = args.output
    software = args.software

def get_geneID_sorted(genes):
    geneID_list = [(geneID,  genes[geneID][0], genes[geneID][3]) for geneID in genes]
    geneID_list = natsorted(geneID_list, key=lambda x: (x[1], x[2]))
    return [i[0] for i in geneID_list]

genes = {}
mRNAs = {}
features = {}
gene2mRNA = {}

if inputfile == '-':
    f = sys.stdin
else:
    f = open(inputfile, 'r')

for l in f:
    if not l.startswith('#') and l.strip() != '':
        l = l.split('\t')
        if l[2] == "gene":
            geneID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
            l[8] = l[8].strip()
            if software!=None:
                l[1] = software
            genes.setdefault(geneID, l)
            
        elif l[2] == "mRNA":
            mRNAID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
            geneID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
            l[8] = l[8].strip()
            if software!=None:
                l[1] = software
            mRNAs.setdefault(mRNAID, l)
            if geneID in gene2mRNA:
                gene2mRNA[geneID].append(mRNAID)
            else:
                gene2mRNA.setdefault(geneID, [mRNAID])
        else:
            mRNAID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
            l[8] = l[8].strip()
            if software!=None:
                l[1] = software
            if mRNAID in features:
                features[mRNAID].append(l)
            else:
                features.setdefault(mRNAID, [l])

if outputfile == None:
    out = sys.stdout
else:
    out = open(outputfile, 'w')
    
for geneID in get_geneID_sorted(genes):

    print(*tuple(genes[geneID]), file=out, sep='\t')
    for mRNAID in gene2mRNA[geneID]:
        print(*tuple(mRNAs[mRNAID]), file=out, sep='\t')
        for feature in features[mRNAID]:
            print(*tuple(feature), file=out, sep='\t')
    print(file=out)
out.close()
