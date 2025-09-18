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

def parser_gtf(gtf):
    transcripts = defaultdict(list)
    transcripts2genes = defaultdict(str)
    if gtf.endswith('.gz'):
        f = gzip.open(gtf, mode='rt')
    else:
        f = open(gtf, 'r')
    
    for l in f:
        if not l.startswith('#') and l.strip() != '':
            l = l.split('\t')
            if l[2] == "exon":
                gene_id = re.search('gene_id "(.*?)";', l[8]).group(1)
                transcript_id = re.search('transcript_id "(.*?)";', l[8]).group(1)
                transcripts2genes.setdefault(transcript_id, gene_id)
                transcripts[transcript_id].append((l[0], min(int(l[3]), int(l[4])), max(int(l[3]), int(l[4])), l[6]))
    f.close()
    return transcripts, transcripts2genes

def extract_sequence(pos, genome_dict):
    sequence = ""
    for p in sorted(pos, key=lambda x: x[1]):
        sequence += genome_dict[p[0]][p[1]-1:p[2]]
    if p[3] == "-":
        sequence = sequence.reverse_complement()
    return sequence

def main(genome, gtf, output=None):
    if output == None:
        out = sys.stdout
    elif output.endswith('.gz'):
        out = gzip(output, 'wt')
    else:
        out = open(output, 'w')
        
    genome_dict = parser_genome(genome)
    transcripts, transcripts2genes = parser_gtf(gtf)
    for t in transcripts:
        print(f">{t} {transcripts2genes[t]}", file=out)
        print(extract_sequence(transcripts[t], genome_dict), file=out)
    out.close()
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get cDNA from genome, similar to gtf_genome_to_cdna_fasta.pl script.",
                                     add_help=False, epilog='date:2025/09/08 author:guisen chen email:thecgs001@foxmail.com',
                                     formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('gtf', metavar='gtf', 
                          help='A file of gtf format.')
    required.add_argument('genome', metavar='genome', 
                          help='A file of genome fasta format.')
    optional.add_argument('-o', '--output', metavar='str', default=None, 
                          help="A file of output. default=None")
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.11', 
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    main(genome=args.genome, gtf=args.gtf, output=args.output)
