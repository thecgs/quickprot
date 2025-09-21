#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
import gzip
import argparse
import subprocess
from Bio import SeqIO
from collections import defaultdict
from argparse import RawTextHelpFormatter

def merge_region(regions):
    sorted_regions = sorted(regions, key=lambda x: x[3], reverse=False)
    merge_region = []
    for region in sorted_regions:
        if (len(merge_region) == 0) or ((region[3] - merge_region[-1][4]) > 1):
            merge_region.append(list(region))
        else:
            merge_region[-1][4] = max(merge_region[-1][4], region[4])
    return merge_region

def get_miniprot_gff3_to_gff3(gff_file, output):
    if gff_file.endswith('.gz'):
        f = gzip.open(gff_file, mode='rt')
    elif gff_file == '-':
        f = sys.stdin
    else:
        f = open(gff_file, 'r')
    
    if output == None:
        out = sys.stdout
    else:
        out = open(output, 'w')

    mRNA2CDSs = None
    for l in f:
        if not l.startswith('#') and l.strip() !='':
            l =  l.strip().split('\t')
            if l[2] == 'mRNA':
                if mRNA2CDSs != None:
                    k = list(mRNA2CDSs.keys())[0]
                    ID = k[8].split(';')[0].split('=')[1].replace('MP', '')
                    Identity = k[8].split(';')[2].split('=')[1]
                    Positive = k[8].split(';')[3].split('=')[1]
                    Target = k[8].split(';')[4].split('=')[1]
                    print('\t'.join((k[0], k[1], 'gene', k[3], k[4], '.', k[6], k[7], f"ID=MG{ID};Identity={Identity};")), file=out)
                    print('\t'.join((k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], 
                                     f"ID=MP{ID};Parent=MG{ID};Identity={Identity};Positive={Positive};Target={Target};")), 
                         file=out)
                    for v in merge_region(regions=mRNA2CDSs[k]):
                        print('\t'.join((v[0], v[1], 'exon', str(v[3]), str(v[4]), '.', v[6], '.',
                                         f"ID=MP{ID}.exon;Parent=MP{ID};")), file=out)
                        print('\t'.join((v[0], v[1], v[2], str(v[3]), str(v[4]), v[5], v[6], v[7],
                                         f"ID=MP{ID}.cds;Parent=MP{ID};")), file=out)
                    print(file=out)
                mRNA = tuple(l)
                mRNA2CDSs = {mRNA:[]}
            elif l[2] == 'CDS' or l[2] == 'stop_codon':
                l[2] = 'CDS'
                l[3] = int(l[3])
                l[4] = int(l[4])
                CDS = l
                mRNA2CDSs[mRNA].append(CDS)
    out.close()
    f.close()
    return None

def get_best_gene_model(genome, gff_file, genetic_code, identity=0.6):
    get_miniprot_gff3_to_gff3(gff_file, output='miniprot_output_tmp.gff3')
    cmd = f"{os.path.join(sys.path[0], 'extract_sequence_from_gff3.py')} miniprot_output_tmp.gff3 {genome} -G {genetic_code} --seqtype prot"
    result=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, encoding='utf-8')
    
    MPID = []
    for record in SeqIO.parse(result.stdout, 'fasta'):
        if record.seq[-1] == "*":
            record.seq = record.seq[:-1]
        if not bool(re.search('\*', str(record.seq))):
            MPID.append(record.id)
            
    mRNA2feature = defaultdict(list)
    mRNA2CDSlength = defaultdict(int)
    
    f = open("miniprot_output_tmp.gff3", 'r')
    for l in f:
        if not l.startswith('#') and l.strip() !='':
            l =  l.strip().split('\t')
            if l[2] == 'gene':
                ID = re.search("ID=(.*?);", l[8]).group(1).replace("MG", "MP")
                Identity = float(re.search('Identity=(.*?);', l[8]).group(1))
                pos = (l[0], int(l[3]), int(l[4]), l[6], l[5], Identity, ID)
                if ID in MPID:
                    mRNA2feature[pos].append('\t'.join(l))
            elif l[2] == 'mRNA' or l[2] == 'exon':
                if ID in MPID:
                    mRNA2feature[pos].append('\t'.join(l))
            elif l[2] == 'CDS':
                if ID in MPID:
                    mRNA2feature[pos].append('\t'.join(l))
                    mRNA2CDSlength[pos] += int(l[4]) - int(l[3]) + 1
    f.close()    
    tmp = sorted([i for i in mRNA2CDSlength if i[5] >= identity], key=lambda x: (x[0], x[1], x[3]))
    filter_res = []
    for i in tmp:
        if (len(filter_res) == 0) or (i[0] != filter_res[-1][0]) or (i[3] != filter_res[-1][3]) or (i[1] > filter_res[-1][2]):
            filter_res.append(i)
        else:
            if i[5] > filter_res[-1][5]:
                filter_res[-1] = i
            elif i[5] == filter_res[-1][5]:
                if mRNA2CDSlength[i] > mRNA2CDSlength[filter_res[-1]]:
                    filter_res[-1] = i      
    res = defaultdict(list)
    for i in filter_res:
        for j in mRNA2feature[i]:
            res[(i[0], i[1], i[2], i[3], mRNA2CDSlength[i])].append(j)
    return res

    
def get_pos(gff_file):
    if gff_file.endswith('.gz'):
        f = gzip.open(gff_file, mode='rt')
    else:
        f = open(gff_file, 'r')
    pos = []
    for l in f:
        if not l.startswith('#') and l.strip() !='':
            l = l.strip().split('\t')
            if l[2] == 'gene':
                pos.append((l[0], int(l[3]), int(l[4]), l[6]))
    f.close()
    return pos

def overlap(pos1, pos2):
    if (pos1[0] == pos2[0]) and (pos1[3] == pos2[3]):
        if (int(pos1[2]) - int(pos2[1]) >=0) and (int(pos2[2]) - int(pos1[1]) >=0):
            return True
    else:
        return False

def get_update(miniprot_pos1, raw_pos2):
    overlap_genes = set()
    for pos1 in miniprot_pos1:
        for pos2 in raw_pos2:
            if overlap(pos1, pos2):
                overlap_genes.add(pos1)
                break
                
    update = []
    for pos in miniprot_pos1:
        if pos not in overlap_genes:
            update.append(pos)
    return update


def main(raw_gff3, miniprot_gff3, genome, output, genetic_code, min_CDS_len=300, identity=0.6):
    raw_gff3 = os.path.realpath(raw_gff3)
    miniprot_gff3 = os.path.realpath(miniprot_gff3)
    genome = os.path.realpath(genome)
    
    raw_pos2 = get_pos(gff_file=raw_gff3)
    res = get_best_gene_model(genome=genome, gff_file=miniprot_gff3, genetic_code=genetic_code, identity=identity)
    update = get_update(miniprot_pos1=res.keys(), raw_pos2=raw_pos2)
    #print(res)
    num_gene = 0
    out = open(output, 'w')
    for k in update:
        if k[4] >= min_CDS_len:
            num_gene += 1
            #print(num_gene)
            for feature in res[k]:
                print(feature, file=out)
            print(file=out)
    out.close()
    print('update gene model: ', str(num_gene))
    os.remove('miniprot_output_tmp.gff3')
    return None
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
Update gff3 from minibusco.
    
Translate Tables/Genetic Codes:
1: The Standard
2: The Vertebrate Mitochondrial
3: The Yeast Mitochondrial
4: The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma
5: The Invertebrate Mitochondrial
6: The Ciliate, Dasycladacean and Hexamita Nuclear
9: The Echinoderm and Flatworm Mitochondrial
10: The Euplotid Nuclear
11: The Bacterial, Archaeal and Plant Plastid
12: The Alternative Yeast Nuclear
13: The Ascidian Mitochondrial
14: The Alternative Flatworm Mitochondrial
15: Blepharisma Macronuclear
16: Chlorophycean Mitochondrial
21: Trematode Mitochondrial
22: Scenedesmus obliquus Mitochondrial
23: Thraustochytrium Mitochondrial
24: Pterobranchia Mitochondrial
25: Candidate Division SR1 and Gracilibacteria
26: Pachysolen tannophilus Nuclear
27: Karyorelict Nuclear
28: Condylostoma Nuclear
29: Mesodinium Nuclear
30: Peritrich Nuclear
31: Blastocrithidia Nuclear
32: Balanophoraceae Plastid
33: Cephalodiscidae Mitochondrial
Reference website: https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes

""", add_help=False, 
    epilog='Date:2024/11/27 Author:Guisen Chen Email:thecgs001@foxmail.com', formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-r', '--raw_gff3', metavar='str', help='A file of gff3 format from raw gff3.', required=True)
    required.add_argument('-m', '--miniprot_gff3', metavar='str', help='A file gff3 format from compleasm, call named "miniprot_output.gff".', required=True)
    required.add_argument('-g', '--genome', metavar='str', help='A genome file fasta format.', required=True)
    optional.add_argument('-o', '--output', metavar='str', help='A output update file of gff3 format. default=None', default=None)
    optional.add_argument('-i', '--identity', metavar='float', type=float, default=0.6, help='Alignment identity (0-1). default=0.6')
    optional.add_argument('-l', '--min_cds_len', metavar='int', type=int, default=300, help='Minimum gene length. default=300')
    optional.add_argument('-G', '--genetic_code', metavar='int', type=int, default=1, help="Genetic code. default=1")
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.11', help="Show program's version number and exit.")
    
    args = parser.parse_args()
    main(raw_gff3=args.raw_gff3, miniprot_gff3=args.miniprot_gff3, genome=args.genome, 
         output=args.output, min_CDS_len=args.min_cds_len, identity=args.identity, genetic_code=args.genetic_code)
