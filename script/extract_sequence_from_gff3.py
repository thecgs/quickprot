#!/usr/bin/env python
# coding: utf-8

import re
import sys
import gzip
import argparse
from Bio import SeqIO, Seq
from collections import defaultdict
from argparse import RawTextHelpFormatter

def parser_genome(genome):
    genome_dict = {}
    if genome.endswith('.gz'):
        for record in SeqIO.parse(gzip.open(genome, mode='rt'), 'fasta'):
            genome_dict.setdefault(record.id, record.seq)
    else:
        for record in SeqIO.parse(genome, 'fasta'):
            genome_dict.setdefault(record.id, record.seq)
    return genome_dict

def parser_gff3(gff3):
    #genes = defaultdict(list)
    mRNAs = defaultdict(list)
    #gene2mRNA = defaultdict(list)
    mRNA2gene = defaultdict(str)
    CDSs = defaultdict(list)
    exons = defaultdict(list)
    
    if gff3.endswith('.gz'):
        f = gzip.open(gff3, mode='rt')
    else:
        f = open(gff3, 'r')
    
    for l in f:
        if not l.startswith('#') and l.strip() != '':
            l = l.split('\t')
            #if l[2] == "gene":
            #    geneID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
            #    genes[geneID] = (l[0], min(int(l[3]), int(l[4])), max(int(l[3]), int(l[4])), l[6])
                
            #elif bool(re.search('RNA', l[2], flags=re.I)):   #include mRNA, tRNA, miRNA ...
            if l[2] == "mRNA":
                mRNAID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
                try:
                    geneID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                    mRNA2gene[mRNAID] = geneID
                except:
                    pass
                mRNAs[mRNAID] = (l[0], min(int(l[3]), int(l[4])), max(int(l[3]), int(l[4])), l[6])
            elif l[2] in ["CDS", "stop_codon"]:
                mRNAID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                CDSs[mRNAID].append((l[0], min(int(l[3]), int(l[4])), max(int(l[3]), int(l[4])), l[6]))
            elif l[2] == "exon":
                mRNAID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                exons[mRNAID].append((l[0], min(int(l[3]), int(l[4])), max(int(l[3]), int(l[4])), l[6]))
    f.close()
    #return genes, mRNAs, gene2mRNA, CDSs, exons
    return mRNAs, mRNA2gene, CDSs, exons

def extract_sequence(pos, genome_dict, cut=False):
    sequence = ""
    for p in sorted(pos, key=lambda x: x[1]):
        #print(p)
        sequence += genome_dict[p[0]][p[1]-1:p[2]]
    if p[3] == "-":
        sequence = sequence.reverse_complement()
    
    if cut:
        if len(sequence) % 3 == 1:
            sequence = sequence[:-1]
        elif len(sequence) % 3 == 2:
            sequence = sequence[:-2]
    return sequence

def main(genome, gff3, output=None, seqtype=['prot', 'CDS'][0], genetic_code=1, remove_stop_codon=False, stop_codon="*"):
    genome_dict = parser_genome(genome)
    mRNAs, mRNA2gene, CDSs, exons = parser_gff3(gff3)
    if output == None:
        out = sys.stdout
    else:
        out = open(output, 'w')
        
    if seqtype == 'CDS':
        for m in mRNAs:
            sequence = extract_sequence(CDSs[m], genome_dict)
            pos_str = f"{mRNAs[m][0]}:{mRNAs[m][1]}-{mRNAs[m][2]}({mRNAs[m][3]})"
            if m in mRNA2gene:
                print(f'>{m} Gene={mRNA2gene[m]} Length={len(sequence)} Position={pos_str}', file=out)
            else:
                print(f'>{m} Length={len(sequence)} Position={pos_str}', file=out)
            print(sequence, file=out)
            
    if seqtype == 'prot':
        for m in mRNAs:
            #print(m)
            sequence = extract_sequence(CDSs[m], genome_dict, cut=True).translate(table=genetic_code, stop_symbol=stop_codon, to_stop=remove_stop_codon)
            #print(sequence)
            if sequence[-1] == stop_codon:
                prot_len = len(sequence) - 1
            else:
                prot_len = len(sequence)
            pos_str = f"{mRNAs[m][0]}:{mRNAs[m][1]}-{mRNAs[m][2]}({mRNAs[m][3]})"
            if m in mRNA2gene:
                print(f'>{m} Gene={mRNA2gene[m]} Length={prot_len} Position={pos_str}', file=out)
            else:
                print(f'>{m} Length={prot_len} Position={pos_str}', file=out)
            print(sequence, file=out)
    out.close()
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
Extract sequence from gff3 file.

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

""", add_help=False, epilog='Date:2025/09/08 Author:Guisen Chen Email:thecgs001@foxmail.com', formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('gff3', metavar='gff3', 
                          help='A file of gff3 format.')
    required.add_argument('genome', metavar='genome', 
                          help='A file of genome fasta format.')
    optional.add_argument('-o', '--output', metavar='str', default=None, 
                          help="A file of output. default=None")
    optional.add_argument('-G', '--genetic_code', metavar='int', type=int, default=1, 
                          help="Genetic code. default=1")
    optional.add_argument('-t', '--seqtype', metavar='str', default='prot', choices=['prot', 'CDS'], 
                          help="Sequence type [prot|CDS]. default=prot")
    optional.add_argument('--remove_stop_codon', action='store_true', 
                          help="Remove stop codon. default=False")
    optional.add_argument('--stop_codon', metavar='str', default='*', 
                          help="Stop codon symbol. default=*")
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.11', 
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    main(genome=args.genome, gff3=args.gff3, output=args.output, seqtype=args.seqtype, genetic_code=args.genetic_code,
         remove_stop_codon=args.remove_stop_codon, stop_codon=args.stop_codon)
