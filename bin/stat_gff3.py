#!/usr/bin/env python3
# coding: utf-8

import os
import re
import sys
import argparse
from collections import defaultdict

def get_Data(gff3_file):
    
    gene_lens = defaultdict(int)
    mRNA_lens = defaultdict(int)
    CDS_lens = defaultdict(int)
    exon_lens = defaultdict(int)
    intron_lens = defaultdict(int)
    gene2mRNA = defaultdict(list)
    
    CDS_nums = defaultdict(int)
    exon_nums = defaultdict(int)
    intron_nums = defaultdict(int)
    
    with open(gff3_file, 'r')  as f:
        for l in f:
            if not l.startswith('#') and l.strip() != '':
                l = l.split('\t')
                if l[2] == "gene":
                    geneID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
                    gene_lens[geneID] = int(l[4]) - int(l[3]) + 1

                elif l[2] == "mRNA":
                #elif bool(re.search('RNA', l[2], flags=re.I)):   #include mRNA, tRNA, miRNA ...
                    mRNAID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
                    geneID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                    gene2mRNA[geneID].append(mRNAID)
                    mRNA_lens[mRNAID] = int(l[4]) - int(l[3]) + 1

                elif l[2] == "CDS":
                    mRNAID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                    if mRNAID in CDS_lens:
                        CDS_lens[mRNAID] += int(l[4]) - int(l[3]) + 1
                    else:
                        CDS_lens[mRNAID] = int(l[4]) - int(l[3]) + 1
                    CDS_nums[mRNAID] += 1
                    
                elif l[2] == "exon":
                    mRNAID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                    if mRNAID in exon_lens:
                        exon_lens[mRNAID] += int(l[4]) - int(l[3]) + 1
                    else:
                        exon_lens[mRNAID] = int(l[4]) - int(l[3]) + 1
                    exon_nums[mRNAID] += 1
                    
    for mRNAID in exon_nums:
        intron_nums[mRNAID] = exon_nums[mRNAID] + 1
        
    for mRNAID in mRNA_lens:
        intron_lens[mRNAID] = mRNA_lens[mRNAID] - exon_lens[mRNAID]
    return gene_lens, mRNA_lens, CDS_lens, exon_lens, intron_lens, CDS_nums, exon_nums, intron_nums, gene2mRNA

def stat_gff3(gff3_file):
    gene_lens, mRNA_lens, CDS_lens, exon_lens, intron_lens, CDS_nums, exon_nums, intron_nums, gene2mRNA = get_Data(gff3_file)

    gene_number = len(gene_lens)
    mRNA_number = len(mRNA_lens)

    Average_CDS_length = sum(CDS_lens.values())/sum(CDS_nums.values())
    Average_exon_length = sum(exon_lens.values())/sum(exon_nums.values())
    Average_intron_length = sum(intron_lens.values())/sum(intron_nums.values())

    Average_gene_length = sum(gene_lens.values())/gene_number
    Average_mRNA_length = sum(mRNA_lens.values())/mRNA_number

    per_gene_Average_CDS_length = sum(CDS_lens.values())/mRNA_number
    per_gene_Average_exon_length = sum(exon_lens.values())/mRNA_number
    per_gene_Average_intron_length = sum(intron_lens.values())/mRNA_number

    Average_CDS_length_per_gene = sum(CDS_nums.values())/mRNA_number
    Average_exon_length_per_gene = sum(exon_nums.values())/mRNA_number
    Average_intron_length_per_gene = sum(intron_nums.values())/mRNA_number

    result = {"Gene number": str(gene_number),
     "Transcript number": str(mRNA_number),
     "CDS Average number per transcript": "{:.2f}".format(Average_CDS_length_per_gene),
     "Exon Average number per transcript": "{:.2f}".format(Average_exon_length_per_gene),
     "Intron Average number per transcript": "{:.2f}".format(Average_intron_length_per_gene),
     "CDS average length(bp)": "{:.2f}".format(Average_CDS_length),
     "Exon average length(bp)": "{:.2f}".format(Average_exon_length),
     "Intron average length(bp)": "{:.2f}".format(Average_intron_length),
     "Gene average length(bp)": "{:.2f}".format(Average_gene_length),
     "Transcript average length(bp)": "{:.2f}".format(Average_mRNA_length),
     "CDS average length for per transcript(bp)": "{:.2f}".format(per_gene_Average_CDS_length),
     "Exon average length for per transcript(bp)": "{:.2f}".format(per_gene_Average_exon_length),
     "Intron average length for per transcript(bp)": "{:.2f}".format(per_gene_Average_intron_length),

    }
    return result

def main(gff3_files, output=None):
    if output==None:
        out = sys.stdout
    else:
        out = open(output, 'w')
    
    outerr = sys.stderr
    print("""# Note:
#
# Gene number: Total number of genes
# Transcript number: Total number of transcripts
# Total length of introns = Total number of transcripts - Total length of exons
# CDS Average number per transcript = Total number of CDSs / Total number of transcripts
# Exon Average number per transcript = Total number of exons / Total number of transcripts
# Intron Average number per transcript = Total number of introns / Total number of transcripts
# CDS average length(bp) = Total length of CDSs / Total number of CDSs
# Exon average length(bp) = Total length of exons / Total number of exons
# Intron average length(bp) = Total length of introns / Total number of introns
# Gene average length(bp)  = Total length of genes / Total number of genes
# Transcript average length(bp) = Total length of transcripts / Total number of transcripts
# CDS average length for per transcript(bp) = Total length of CDSs / Total number of transcripts
# Exon average length for per transcript(bp) = Total length of exons / Total number of transcripts
# Intron average length for per transcript(bp) = Total length of introns / Total number of transcripts
""", file=outerr)

    print("Name", "Gene number", "Transcript number", 
          "CDS Average number per transcript", "Exon Average number per transcript", "Intron Average number per transcript",
          "CDS average length(bp)", "Exon average length(bp)", "Intron average length(bp)",
          "Gene average length(bp)", "Transcript average length(bp)",
          "CDS average length for per transcript(bp)",
          "Exon average length for per transcript(bp)",
          "Intron average length for per transcript(bp)",
          sep='\t', file=out)
    
    for gff3_file in gff3_files:
        prefix = os.path.splitext(os.path.basename(gff3_file))[0]
        stat = stat_gff3(gff3_file)
        print(prefix, 
              stat["Gene number"],
              stat["Transcript number"],
              stat["CDS Average number per transcript"],
              stat["Exon Average number per transcript"],
              stat["Intron Average number per transcript"],
              stat["CDS average length(bp)"],
              stat["Exon average length(bp)"],
              stat["Intron average length(bp)"],
              stat["Gene average length(bp)"],
              stat["Transcript average length(bp)"],
              stat["CDS average length for per transcript(bp)"],
              stat["Exon average length for per transcript(bp)"],
              stat["Intron average length for per transcript(bp)"],
              sep='\t', file=out)

    out.close()
    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='stat gff3 file.',  add_help=False,
                                     epilog='date:2024/11/08 author:guisen chen email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('input', metavar='gff3', nargs="+", 
                          help='A or multiple files of gff3 format.')
    optional.add_argument('-o', '--output', metavar='str', default=None,  
                          help=f'A file of tsv format. defualt=None')
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00',  
                          help="Show program's version number and exit.")
    args = parser.parse_args()
    main(gff3_files=args.input, output=args.output)


