#!/usr/bin/env python
# coding: utf-8

import re
import sys
import gzip
import argparse
from collections import defaultdict

def get_longest_transcript(inputfile, outputfile=None, header=None):
    
    genes = defaultdict(list)
    mRNAs = defaultdict(list)
    gene2mRNA = defaultdict(list)
    CDS_lens = defaultdict(int)
    features = defaultdict(list)
    
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
                genes[geneID] = l
                
            #elif bool(re.search('RNA', l[2], flags=re.I)):   #include mRNA, tRNA, miRNA ...
            elif l[2] == "mRNA":
                mRNAID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
                geneID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                gene2mRNA[geneID].append(mRNAID)
                mRNAs[mRNAID] = l
            
            elif l[2] == "CDS":
                mRNAID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                if mRNAID in CDS_lens:
                    CDS_lens[mRNAID] += int(l[4]) - int(l[3]) + 1
                else:
                    CDS_lens[mRNAID] = int(l[4]) - int(l[3]) + 1
                features[mRNAID].append(l)
            
            else:
                try:
                    mRNAID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                    features[mRNAID].append(l)
                except:
                    pass
    f.close()

    no_CDS_mRNA = []  ##remove no CDS mRNA from NCBI download gff3
    for mRNA in mRNAs:
        if mRNA not in CDS_lens:
            no_CDS_mRNA.append(mRNA)
    
    if outputfile == None:
        out = sys.stdout
    elif outputfile.endswith('.gz'):
        out = gzip.open(outputfile, 'wt')
    else:
        out = open(outputfile, 'w')
    
    if header != None:
        info = ''.join(["#" + i for i in header.split('\n')])
        print(info, file=out)
        print(file=out)
    
    for gene in gene2mRNA:
        longest_CDS_len = 0
        longest_mRNA = ''
        
        no_CDS_staus = False
        for mRNA in gene2mRNA[gene]:
            if mRNA in no_CDS_mRNA:
                no_CDS_staus = True
                continue

            if longest_CDS_len < CDS_lens[mRNA]:
                longest_CDS_len = CDS_lens[mRNA]
                longest_mRNA = mRNA
        if no_CDS_staus:
            continue
        
        print(*tuple(genes[gene]), sep='\t', end='', file=out)
        print(*tuple(mRNAs[longest_mRNA]), sep='\t', end='', file=out)
        for feature in features[longest_mRNA]:
            print(*tuple(feature), sep='\t', end='', file=out)
        print('', file=out)
    out.close()
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get longest trascript from gff3 base on CDS length.',  add_help=False,
                                     epilog='Date:2024/08/04 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('input', metavar='gff3',
                          help='A input file of gff3 format.')
    optional.add_argument('-o', '--output', metavar='str', default=None,  
                          help=f'A output file of gff3 format. defualt=None')
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('--header', metavar='str', type=str, default=None,
                          help='Add header information to a gff3 file. default=None')
    optional.add_argument('-v', '--version', action='version', version='v1.00',  
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    
    get_longest_transcript(inputfile=args.input, outputfile=args.output, header=args.header)
    
