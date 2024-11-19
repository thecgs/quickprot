#!/usr/bin/env python
# coding: utf-8

import re
import argparse
from collections import defaultdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="split and filter gene model.", add_help=False, 
                                     epilog='date:2024/11/19 author:guisen chen email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', metavar='str', help='A input file of gff3 format.', required=True)
    required.add_argument('-o', '--output', metavar='str', help='A output file of gff3 format.', required=True)
    optional.add_argument('--overlap', metavar='float', type=float, default=0.8, help='If the overlap of predicted ORFs in a transcript is less than default (default: 0.8, range: 0-1), they will be dissected.')
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00', help="Show program's version number and exit.")
    args = parser.parse_args()
    infile = args.input
    outfile = args.output
    overlap = args.overlap

def get_filter_regions(regions, overlap=0.8):
    filter_regions = []
    regions = sorted(regions, key=lambda x: (x[0],x[1]))
    for region in regions:
        region = list(region)
        if len(filter_regions) == 0 or region[0] != filter_regions[-1][0] or region[1] > filter_regions[-1][2]:
            filter_regions.append(region)
        else:
            if region[3] == filter_regions[-1][3]:
                value = (region[1] - filter_regions[-1][1] +1)/(max(region[2], filter_regions[-1][2]) - min(region[1], filter_regions[-1][1]) +1)
                if overlap <= value:
                    filter_regions[-1][2] = max(region[2], filter_regions[-1][2])
                else:
                    if region[6] > filter_regions[-1][6]:
                        filter_regions[-1] = region
            else:
                if region[5] > filter_regions[-1][5]:
                    filter_regions[-1] = region
    return filter_regions

mRNA2CDS = defaultdict(list)
with open(infile, 'r') as f:
    for l in f:
        if not l.startswith('#') and l.strip() != '':
            l = l.strip().split('\t')
            if l[2] == 'mRNA':
                ID = re.search('ID=(.*?);',l[8]).group(1)
                TYPE = re.search('Name="ORF type:(.*?) \(',l[8]).group(1)
                SCORE = float(re.search(',score=(.*?)"',l[8]).group(1))
                Key = (l[0],l[6], ID, TYPE, SCORE)
            if l[2] == 'CDS':
                Value = (l[0], int(l[3]), int(l[4]), l[6], int(l[7]))
                mRNA2CDS[Key].append(Value)

new_mRNA2CDS = defaultdict(list)
for k in mRNA2CDS:
    Key = (k[0], min([v[1] for v in mRNA2CDS[k]]), max([v[2] for v in mRNA2CDS[k]]), 
           k[1], k[3], k[4], sum([v[2] - v[1] +1 for v in mRNA2CDS[k]]))
    new_mRNA2CDS[Key] = mRNA2CDS[k]
filter_regions = get_filter_regions(new_mRNA2CDS, overlap)

cluster_genes = defaultdict(list)
for g in filter_regions:
    for m in new_mRNA2CDS:
        if g[0] == m[0] and g[3] == m[3] and g[1]<=m[1] and m[2]<=g[2]:
            cluster_genes[tuple(g)].append(m)
            
out = open(outfile, 'w')
ti = 0
for gi, g in enumerate(cluster_genes):
    print(g[0], 'quickprot', 'gene', g[1], g[2], '.', g[3], '.', f'ID=QUKPGENE{gi+1};', file=out, sep='\t')
    for mi, m in enumerate(cluster_genes[g]):
        ti += 1
        print(m[0], 'quickprot', 'mRNA', m[1], m[2], m[5], m[3], '.', f'ID=QUKPMRNA{gi+1}.p{mi+1};Parent=QUKPGENE{gi+1};Type={m[4]};', file=out, sep='\t')
        for ei, e in enumerate(new_mRNA2CDS[m]):
            print(e[0], 'quickprot', 'exon', e[1], e[2], '.', e[3], '.', f'ID=QUKPMRNA{gi+1}.p{mi+1}.exon{ei+1};Parent=QUKPMRNA{gi+1}.p{mi+1};', file=out, sep='\t')
            print(e[0], 'quickprot', 'CDS', e[1], e[2], '.', e[3], e[4], f'ID=QUKPMRNA{gi+1}.p{mi+1}.cds{ei+1};Parent=QUKPMRNA{gi+1}.p{mi+1};', file=out, sep='\t')
    print(file=out)

print('Predicted number of genes:', gi+1)
print('Predicted number of mRNAs:', ti)
out.close()
