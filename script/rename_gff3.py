#!/usr/bin/env python
# coding: utf-8

import re
import sys
import gzip
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rename for gene ID from gff3 format file.',  add_help=False,
                                     epilog='date:2024/11/07 author:guisen chen email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('input', metavar='gff3', 
                          help='A input file of gff3 format.')
    optional.add_argument('-o', '--output', metavar='str', default=None,  
                          help=f'A output file of gff3 format. defualt=None')
    optional.add_argument('-p', '--prefix', metavar='str', default='',  
                          help=f'A prefix of geneID. defualt=None')
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00',  
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    prefix = args.prefix
    inputfile = args.input
    outputfile = args.output
    
genes = {}
mRNAs = {}
features = {}
gene2mRNA = {}

if inputfile == '-':
    f = sys.stdin
elif inputfile.endswith('.gz'):
    f = gzip.open(inputfile, 'rt')
else:
    f = open(inputfile, 'r')
    
for l in f:
    if not l.startswith('#') and l.strip() != '':
        l = l.split('\t')
        if l[2] == "gene":
            geneID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
            genes.setdefault(geneID, l[0:8])
            
        elif l[2] == "mRNA":
            mRNAID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
            geneID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
            mRNAs.setdefault(mRNAID, l[0:8])
            if geneID in gene2mRNA:
                gene2mRNA[geneID].append(mRNAID)
            else:
                gene2mRNA.setdefault(geneID, [mRNAID])
        else:
            mRNAID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
            if mRNAID in features:
                features[mRNAID].append(l[0:8])
            else:
                features.setdefault(mRNAID, [l[0:8]])

prefixs = {"three_prime_UTR":"utr3",
          "five_prime_UTR":"utr5"}

if outputfile == None:
    out = sys.stdout
elif outputfile.endswith('.gz'):
    out = open(outputfile, 'wt')
else:
    out = open(outputfile, 'w')
    
for index, geneID in enumerate(gene2mRNA):
    i = "{:0>8}".format(index+1)
    print(*tuple(genes[geneID]), f"ID={prefix}GENE{i};", file=out, sep='\t')
    for t, mRNAID in enumerate(gene2mRNA[geneID]):
        print(*tuple(mRNAs[mRNAID]), f"ID={prefix}MRNA{i}.{t+1};Parent={prefix}GENE{i};", file=out, sep='\t')
        u3 = 1
        u5 = 1
        e = 1
        for feature in features[mRNAID]:
            if feature[2] == 'exon':
                print(*tuple(feature), f"ID={prefix}MRNA{i}.{t+1}.exon{e};Parent={prefix}MRNA{i}.{t+1};", file=out, sep='\t')
                e += 1
            elif feature[2] == 'three_prime_UTR':
                print(*tuple(feature), f"ID={prefix}MRNA{i}.{t+1}.utr3.p{u3};Parent={prefix}MRNA{i}.{t+1};", file=out, sep='\t')
                u3 += 1
            elif feature[2] == 'five_prime_UTR':
                print(*tuple(feature), f"ID={prefix}MRNA{i}.{t+1}.utr5.p{u5};Parent={prefix}MRNA{i}.{t+1};", file=out, sep='\t')
                u5 += 1
            else:
                print(*tuple(feature), f"ID={prefix}MRNA{i}.{t+1}.{prefixs.get(feature[2], feature[2].lower())};Parent={prefix}MRNA{i}.{t+1};", file=out, sep='\t')
    print(file=out)
out.close()
