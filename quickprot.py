#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
import glob
import shutil
import argparse
import subprocess
from collections import defaultdict
from argparse import RawTextHelpFormatter

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Quick genome annotation base on protein.", add_help=False, 
                                     epilog='date:2024/11/19 author:guisen chen email:thecgs001@foxmail.com', 
                                     formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-q', '--query', metavar='str', help='A file of query protein fasta format.', required=True)
    required.add_argument('-g', '--genome', metavar='str', help='A file of genome fasta format.', required=True)
    optional.add_argument('-p', '--prefix', metavar='str', default='quickprot', help='Prefix of a output file. default=quickprot')
    optional.add_argument('-i', '--identity', metavar='float', type=float, default=0.95, help='Alignment identity (0-1). default=0.95')
    optional.add_argument('--outs', metavar='float', type=float, default=0.99, help='Output score at least bestScore (0-1). default=0.99')
    optional.add_argument('--overlap', metavar='float', type=float, default=0.8, help="""If the overlap of predicted ORFs in a transcript is less than default value (0-1). default=0.8, 
they will be dissected.""")
    optional.add_argument('-t', '--thread', metavar='int', type=int, default=os.cpu_count(), help=f'Thread number of run miniprot sortware. defualt={os.cpu_count()}')
    optional.add_argument('-G', '--genetic_code', metavar='str', type=str, default='Universal', help=f"""Genetic Codes (derived from: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). defualt=Universal
The supported genetic codon tables are Acetabularia, Candida, Ciliate, Dasycladacean, Euplotid, Hexamita,
Mesodinium, Mitochondrial-Ascidian, Mitochondrial-Chlorophycean, Mitochondrial-Echinoderm, Mitochondrial-Flatworm,
Mitochondrial-Invertebrates, Mitochondrial-Protozoan, Mitochondrial-Pterobranchia, Mitochondrial-Scenedesmus_obliquus,
Mitochondrial-Thraustochytrium, Mitochondrial-Trematode, Mitochondrial-Vertebrates, Mitochondrial-Yeast, 
Pachysolen_tannophilus,Peritrich, SR1_Gracilibacteria, Tetrahymena, and Universal.
""")
    optional.add_argument('-s', '--skip_align', action='store_true', help="Skip run miniprot step. default=False")
    optional.add_argument('-m', '--mask', action='store_true', help="Soft-masked (dna_sm) genome convert to masked(dna_rm) genome. default=False")
    optional.add_argument('-n', '--noclean', action='store_true', help="Do not delete intermediate files. default=False")
    optional.add_argument('-b', '--single_best_only', action='store_true', help="""Retain only the single best orf per transcript. default=False
It is not recommended to use it because when two reference proteins overlap during alignment, 
it can lead to fusion during transcript assembly. If a transcript is not set with only one ORF,
the fused ORF will be split in subsequent analysis.""")
    optional.add_argument('-h', '--help', action='help', help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.11', help="Show program's version number and exit.")
    args = parser.parse_args()
    genetic_code = args.genetic_code
    query_file = os.path.realpath(args.query)
    genome_file = os.path.realpath(args.genome)
    thread = args.thread
    mask = args.mask
    identity = args.identity
    prefix = args.prefix
    single_best_only = args.single_best_only
    skip_align = args.skip_align
    noclean = args.noclean
    outs = args.outs
    overlap = args.overlap

def rm(file):
    if os.path.isdir(file):
        shutil.rmtree(file, ignore_errors=True)
    else:
        try:
            os.remove(file)
        except FileNotFoundError:
            pass
    return None

def run_miniprot(query_file, genome_file, output, thread, mask, skip_align, outs):
    if mask:
        cmd = f"{os.path.join(sys.path[0], 'bin', 'sm2rmForFasta.py')} -i {genome_file} -o {genome_file}.tmp"
        subprocess.run(cmd, shell=True)
        genome_file = os.path.realpath(genome_file+'.tmp')
    
    miniprot  = os.path.join(sys.path[0], 'bin', 'miniprot')
    cmd = f'{miniprot} -I --outs={outs} -t {thread} --aln --gff {genome_file} {query_file} > {output}'
    if skip_align:
        pass
    else:
        subprocess.run(cmd, shell=True)
    
    if mask:
        try:
            os.remove(genome_file+'.tmp')
        except:
            pass
    return None

def merge_region(regions):
    sorted_regions = sorted(list(regions), key=lambda x: (x[3], x[0], x[1]), reverse=False)
    merge_region = []
    for region in sorted_regions:
        if (len(merge_region) == 0) or (region[0] != merge_region[-1][0]) or (region[3] != merge_region[-1][3]) \
        or (region[1] > merge_region[-1][2]):
            merge_region.append(list(region))
        else:
            merge_region[-1][2] = max(merge_region[-1][2], region[2])
            #merge_region[-1][1] = min(merge_region[-1][1], region[1])
    return merge_region


miniprot_output = 'miniprot_output.gff3'
run_miniprot(query_file=query_file, genome_file=genome_file, output=miniprot_output, 
             thread=thread, mask=mask, skip_align=skip_align, outs=args.outs)

transcript_regions = set()
exon_regions = set()
with open(miniprot_output, 'r') as f:
    for l in f:
        if not l.startswith('#') and l.strip() !='':
            l = l.strip().split('\t')
            if l[2] == 'mRNA':
                _status = False
                if float(re.search('Identity=(.*?);', l[8]).group(1)) >= identity:
                    transcript_region = (l[0], int(l[3]), int(l[4]), l[6])
                    transcript_regions.add(transcript_region)
                else:
                    _status = True
            if _status == False:
                if l[2] == 'CDS':
                    #if float(re.search('Identity=(.*?);', l[8]).group(1)) >= identity:
                    exon_region = (l[0], int(l[3]), int(l[4]), l[6])
                    exon_regions.add(exon_region)

transcript_regions = merge_region(transcript_regions)
exon_regions = merge_region(exon_regions)
transcript_regions = sorted(list(transcript_regions), key=lambda x: (x[0], x[1]), reverse=False)
exon_regions = sorted(list(exon_regions), key=lambda x: (x[0], x[1]), reverse=False)

clusters = defaultdict(list)
for exon_region in exon_regions:
    for transcript_region in transcript_regions:
        if transcript_region[0] == exon_region[0] and transcript_region[3] == exon_region[3] and transcript_region[1] <= exon_region[1] and exon_region[2] <= transcript_region[2]:
            clusters[tuple(transcript_region)].append(tuple(exon_region))
            break
print('Assemble {} transcripts'.format(len(clusters)))

out = open('quickprot.transcript.gtf', 'w')
for index, transcript in enumerate(clusters):
    print(transcript[0], 'quickprot', 'transcript', transcript[1], transcript[2], '.', transcript[3], '.',
          'gene_id "QUKPGENE{}"; transcript_id "QUKPMRNA{}";'.format(index+1, index+1), sep='\t', file=out)
    for exon in clusters[transcript]:
        print(exon[0], 'quickprot', 'exon', exon[1], exon[2], '.', exon[3], '.',
              'gene_id "QUKPGENE{}"; transcript_id "QUKPMRNA{};'.format(index+1, index+1), sep='\t', file=out)
out.close()

cmd = f"{os.path.join(sys.path[0], 'bin/TransDecoder-5.7.1/util/gtf_to_alignment_gff3.pl')} quickprot.transcript.gtf > quickprot.transcript.gff3"
subprocess.run(cmd, shell=True, capture_output=True)
cmd = f"{os.path.join(sys.path[0], 'bin/TransDecoder-5.7.1/util/gtf_genome_to_cdna_fasta.pl')} quickprot.transcript.gtf {genome_file} > quickprot.transcript.fasta"
subprocess.run(cmd, shell=True, capture_output=True)
cmd = f"{os.path.join(sys.path[0], 'bin/TransDecoder-5.7.1/TransDecoder.LongOrfs')} -t quickprot.transcript.fasta --genetic_code {genetic_code}"
subprocess.run(cmd, shell=True, capture_output=True)
cmd = f"{os.path.join(sys.path[0], 'bin/TransDecoder-5.7.1/TransDecoder.Predict')} -t quickprot.transcript.fasta --genetic_code {genetic_code}"
if single_best_only == True:
    cmd += ' --single_best_only'
subprocess.run(cmd, shell=True, capture_output=True)

cmd = f"{os.path.join(sys.path[0], 'bin/TransDecoder-5.7.1/util/cdna_alignment_orf_to_genome_orf.pl')} \
quickprot.transcript.fasta.transdecoder.gff3 quickprot.transcript.gff3 quickprot.transcript.fasta > quickprot.transcript.genome.gff3"
subprocess.run(cmd, shell=True, capture_output=True)

cmd = f"{os.path.join(sys.path[0], 'bin/split_and_filter_gene_model.py')} \
-i quickprot.transcript.genome.gff3 -o {prefix}.gff3 --overlap {overlap}"
subprocess.run(cmd, shell=True)

cmd = f"{os.path.join(sys.path[0], 'bin/TransDecoder-5.7.1/util/gff3_file_to_proteins.pl')} \
--gff3 {prefix}.gff3 --fasta {genome_file} --genetic_code {genetic_code} --seqType prot > {prefix}.pep.fasta"
subprocess.run(cmd, shell=True, capture_output=True)

cmd = f"{os.path.join(sys.path[0], 'bin/TransDecoder-5.7.1/util/gff3_file_to_proteins.pl')} \
--gff3 {prefix}.gff3 --fasta {genome_file} --genetic_code {genetic_code} --seqType CDS > {prefix}.cds.fasta"
subprocess.run(cmd, shell=True, capture_output=True)

cmd = f"{os.path.join(sys.path[0], 'bin/get_longest_transcript_gff3.py')} {prefix}.gff3 -o {prefix}.longest.gff3"
subprocess.run(cmd, shell=True, capture_output=True)

cmd = f"{os.path.join(sys.path[0], 'bin/TransDecoder-5.7.1/util/gff3_file_to_proteins.pl')} \
--gff3 {prefix}.longest.gff3 --fasta {genome_file} --genetic_code {genetic_code} --seqType prot > {prefix}.longest.pep.fasta"
subprocess.run(cmd, shell=True, capture_output=True)

cmd = f"{os.path.join(sys.path[0], 'bin/TransDecoder-5.7.1/util/gff3_file_to_proteins.pl')} \
--gff3 {prefix}.longest.gff3 --fasta {genome_file} --genetic_code {genetic_code} --seqType CDS > {prefix}.longest.cds.fasta"
subprocess.run(cmd, shell=True, capture_output=True)

intermediate_files = [#os.path.realpath('miniprot_output.gff3'),
                      os.path.realpath('quickprot.transcript.gtf'),
                      os.path.realpath('quickprot.transcript.gff3'),
                      os.path.realpath('quickprot.transcript.genome.gff3'),
                      *tuple(glob.glob(os.path.realpath('quickprot.transcript.fasta')+'*'))]
if noclean:
    pass
else:
    for file in intermediate_files:
        print(file)
        rm(file)

