#!/usr/bin/env python
# coding: utf-8

import re
import sys
import gzip
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get an intron feature from gff3 file.', 
                                     add_help=False, epilog='Date:2025/09/21 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('gff3', metavar='gff3', help='A input file of gff3 format.')
    optional.add_argument('-o', '--output', metavar='str', help='A output file of gff3 format.')
    optional.add_argument('-r', '--retain_original_features', action='store_true', help='Retain original features.')
    optional.add_argument('-h','--help', action='help', help='Show this help message and exit.')
    optional.add_argument('-v','--version', action='version', version='v1.00', help="Show program's version number and exit.")
    args = parser.parse_args()
    inputfile = args.gff3
    outputfile = args.output
    retain_original_features = args.retain_original_features

#inputfile = "/home/chenguisen/test/Cromileptes_altivelis.longest.gff3.gz"

def get_introns(mRNA2exons):
    mRNA2introns = {}
    for mRNAID in mRNA2exons:
        mRNA2introns.setdefault(mRNAID, [])
        pos = sorted(mRNA2exons[mRNAID], key=lambda x: x[0])
        for i in range(1, len(pos)):
            start = pos[i-1][1] + 1
            end = pos[i][0] - 1
            if end - start > 0:
                mRNA2introns[mRNAID].append((start, end, pos[0][2]))
    return mRNA2introns

genes = {}
mRNAs = {}
features = {}
gene2mRNA = {}
mRNA2exons = {}

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
            genes.setdefault(geneID, l[0:8] + [l[8].strip()])
            
        elif l[2] == "mRNA":
            mRNAID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
            geneID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
            mRNAs.setdefault(mRNAID, l[0:8] + [l[8].strip()])
            if geneID in gene2mRNA:
                gene2mRNA[geneID].append(mRNAID)
            else:
                gene2mRNA.setdefault(geneID, [mRNAID])
        else:
            mRNAID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
            if mRNAID in features:
                features[mRNAID].append(l[0:8] + [l[8].strip()])
            else:
                features.setdefault(mRNAID, [l[0:8] + [l[8].strip()]])
                
            if l[2] == 'exon':
                pos = (int(l[3]), int(l[4]), l[6])
                if mRNAID in mRNA2exons:
                    mRNA2exons[mRNAID].append(pos)
                else:
                    mRNA2exons.setdefault(mRNAID, [pos])

mRNA2introns = get_introns(mRNA2exons)

if outputfile == None:
    out = sys.stdout
elif outputfile.endswith('.gz'):
    out = gzip.open(outputfile, 'wt')
else:
    out = open(outputfile, 'w')

for g in gene2mRNA:
    print('\t'.join(genes[g]), file=out)
    for m in gene2mRNA[g]:
        print('\t'.join(mRNAs[m]), file=out)
        Chr = mRNAs[m][0] 
        Source = mRNAs[m][1] 
        Strand = mRNAs[m][6]
        if retain_original_features:
            for feature in features[m]:
                print('\t'.join(feature), file=out)
        if Strand == '-':
            for n, i in enumerate(sorted(mRNA2introns[m], key=lambda x: x[1], reverse=True)):
                print(Chr, Source, 'intron', i[0], i[1], '.', i[2], '.', f'ID={m}.intron.{n+1};Parent={m};', sep='\t', file=out)
        else:
            for n, i in enumerate(mRNA2introns[m]):
                print(Chr, Source, 'intron', i[0], i[1], '.', i[2], '.', f'ID={m}.intron.{n+1};Parent={m};', sep='\t', file=out)
    print('', file=out)
