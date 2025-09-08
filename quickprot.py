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
    parser = argparse.ArgumentParser(description="""
Quick genome annotation base on protein.

Translate Tables/Genetic Codes:
1: The Standard
2: The Vertebrate Mitochondrial
3: The Yeast Mitochondrial
4: The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma
5: The Invertebrate Mitochondrial
6: The Ciliate, Dasycladacean and Hexamita Nuclear
9: The Echinoderm and Flatworm Mitochondrial
10: The Euplotid Nuclear
11: The Bacterial, Archaeal and Plant Plastid        ## TransDecoder not supported
12: The Alternative Yeast Nuclear
13: The Ascidian Mitochondrial
14: The Alternative Flatworm Mitochondrial           ## TransDecoder not supported
15: Blepharisma Macronuclear
16: Chlorophycean Mitochondrial
21: Trematode Mitochondrial                          
22: Scenedesmus obliquus Mitochondrial                 
23: Thraustochytrium Mitochondrial                   
24: Pterobranchia Mitochondrial                      
25: Candidate Division SR1 and Gracilibacteria       
26: Pachysolen tannophilus Nuclear                   
27: Karyorelict Nuclear                              ## TransDecoder not supported
28: Condylostoma Nuclear                             ## TransDecoder not supported
29: Mesodinium Nuclear
30: Peritrich Nuclear
31: Blastocrithidia Nuclear                          ## TransDecoder not supported
32: Balanophoraceae Plastid                          ## TransDecoder not supported
33: Cephalodiscidae Mitochondrial                    ## TransDecoder not supported
Reference website: https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes
""", add_help=False, epilog='date:2024/11/19 author:guisen chen email:thecgs001@foxmail.com', formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-q', '--query', metavar='str', help='A file of query protein fasta format, supports .gz compressed files.', required=True)
    required.add_argument('-g', '--genome', metavar='str', help='A file of genome fasta format, supports .gz compressed files.', required=True)
    optional.add_argument('-p', '--prefix', metavar='str', default='quickprot', help='Prefix of a output file. default=quickprot')
    optional.add_argument('-i', '--identity', metavar='float', type=float, default=0.8, help='Alignment identity (0-1). default=0.8')
    optional.add_argument('--outs', metavar='float', type=float, default=0.95, help='Output score at least bestScore (0-1). default=0.95')
    optional.add_argument('--overlap', metavar='float', type=float, default=0.8, help="""If the overlap of predicted ORFs in a transcript is less than default value (0-1). default=0.8, 
they will be dissected.""")
    optional.add_argument('-t', '--thread', metavar='int', type=int, default=os.cpu_count(), help=f'Thread number of run miniprot sortware. defualt={os.cpu_count()}')
    optional.add_argument('-G', '--genetic_code', metavar='int', type=int, default=1, help="Genetic code. default=1")
    optional.add_argument('-s', '--skip_align', action='store_true', help="Skip run miniprot step. default=False")
    optional.add_argument('-m', '--mask', action='store_true', help="Soft-masked (dna_sm) genome convert to masked(dna_rm) genome. default=False")
    optional.add_argument('-n', '--noclean', action='store_true', help="Do not delete intermediate files. default=False")
    optional.add_argument('-b', '--single_best_only', action='store_true', help="""Retain only the single best orf per transcript. default=False
It is not recommended to use it because when two reference proteins overlap during alignment, 
it can lead to fusion during transcript assembly. If a transcript is not set with only one ORF,
the fused ORF will be split in subsequent analysis.""")
    optional.add_argument('-miniprot_PATH', '--miniprot_PATH', metavar='str', help="miniprot PATH default=auto.", default=None)
    optional.add_argument('-TransDecoder_PATH', '--TransDecoder_PATH', metavar='str', help="TransDecoder PATH default=auto.", default=None)
    optional.add_argument('-ORFSoftware', '--ORFSoftware', metavar='str', 
                          help="Tool for selecting predicted ORFs, TransDecoder or TD2. default=TransDecoder", default="TransDecoder")
    optional.add_argument('--debug_info', action='store_true', help="Display all software output details. default=False")
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
    TransDecoder_PATH = args.TransDecoder_PATH
    miniprot_PATH = args.miniprot_PATH
    ORFSoftware = args.ORFSoftware
    debug_info = args.debug_info

NCBI2TransDecoder_genetic_code = {1: "Universal",
                                  2: "Mitochondrial-Vertebrates",
                                  3: "Mitochondrial-Yeast",
                                  4: "Mitochondrial-Protozoan",
                                  5: "Mitochondrial-Invertebrates",
                                  6: "Ciliate",                    # Acetabularia, Dasycladacean, Hexamita, Tetrahymena
                                  9: "Mitochondrial-Echinoderm",   # Mitochondrial-Flatworm
                                  10: "Euplotid",
                                  12: "Candida",
                                  13: "Mitochondrial-Ascidian",
                                  16: "Mitochondrial-Chlorophycean",
                                  21: "Mitochondrial-Trematode",
                                  22: "Mitochondrial-Scenedesmus_obliquus",
                                  23: "Mitochondrial-Thraustochytrium",
                                  24: "Mitochondrial-Pterobranchia",
                                  25: "SR1_Gracilibacteria",
                                  26: "Pachysolen_tannophilus",
                                  29: "Mesodinium",
                                  30: "Peritrich",
}

def check_dependencies(ORFSoftware, miniprot_PATH=None, TransDecoder_PATH=None):
    print("Check dependencies...")
    if miniprot_PATH==None:
        if os.path.exists(os.path.join(sys.path[0], 'bin/miniprot-0.12/miniprot')):
            miniprot_PATH = os.path.join(sys.path[0], 'bin/miniprot-0.12/miniprot')
        else:
            for PATH in os.environ['PATH'].split(':'):
                if os.path.exists(os.path.join(PATH, 'miniprot')):
                    miniprot_PATH = os.path.join(PATH, 'miniprot')
                    break

        if miniprot_PATH == None:
            print('Error: miniprot does not exist!')
        else:
            print(f'miniprot PATH: {miniprot_PATH}')
            
    if TransDecoder_PATH == None:
        if os.path.exists(os.path.join(sys.path[0], 'bin/TransDecoder-5.7.1/TransDecoder.Predict')):
            TransDecoder_PATH = os.path.join(sys.path[0], 'bin/TransDecoder-5.7.1/')
        else:
            for PATH in os.environ['PATH'].split(':'):
                if os.path.exists(os.path.join(PATH, 'TransDecoder.Predict')):
                    TransDecoder_PATH = PATH
                    break
        if TransDecoder_PATH == None:
            print('Error: TransDecoder does not exist!')
        else:
            print(f'TransDecoder PATH: {TransDecoder_PATH}')
        
    if miniprot_PATH == None or TransDecoder_PATH == None:
        sys.exit()
        
    if ORFSoftware=="TransDecoder":
        pass
    elif ORFSoftware=="TD2":
        try:
            import TD2
            print("TD2 installed.")
        except ModuleNotFoundError:
            print('Plase install TD2. Install command: `pip3 install TD2`.')

    try:
        from Bio import SeqIO
        print("biopython installed.")
    except ModuleNotFoundError:
        print('Plase install biopython. Install command: `pip3 install biopython`.')
    return miniprot_PATH, TransDecoder_PATH

def run_cmd(cmd, jobname=None, capture_output=False):
    
    import sys
    import time
    import subprocess
    
    def htime(t):
        if 0 <= round(t,2) < 60:
            res = '{:.2f}{}'.format(t, 's')
        elif 60 <= round(t,2) < 60*60:
            res = '{:.2f}{}'.format(t/60, 'm')
        elif 60*60 <= round(t,2) < 60*60*60:
            res = '{:.2f}{}'.format(t/60/60, 'h')
        else:
            res = '{:.2f}{}'.format(t/60/60/24, 'd')
        return res
    
    log = sys.stderr
    
    #start_time = time.gmtime()
    start_time = time.localtime()
    ftime = time.strftime('%Y-%m-%d %H:%M:%S', start_time)
    print(f"\033[033m[Running {ftime}]: {cmd}\033[0m", file=log)
    if jobname==None:
        res = subprocess.run(cmd, shell=True, capture_output=capture_output)
        #end_time = time.gmtime()
        end_time = time.localtime()
        ftime = time.strftime('%Y-%m-%d %H:%M:%S', end_time)
        if res.returncode == 0:
            print(f"\033[032m[Finished {ftime} {htime(time.mktime(end_time)-time.mktime(start_time))}]: {cmd}\033[0m", file=log)
        else:
            print(f"\033[031m[Error {ftime} {htime(time.mktime(end_time)-time.mktime(start_time))}]: {cmd}\033[0m", file=log)
    else:
        if not os.path.exists(f'{jobname}.done'):
            res = subprocess.run(cmd, shell=True, capture_output=capture_output)
            #end_time = time.gmtime()
            end_time = time.localtime()
            if res.returncode == 0:
                out = open(f'{jobname}.done', 'w')
                out.close()
                print(f"\033[032m[Finished {ftime} {htime(time.mktime(end_time)-time.mktime(start_time))}]: {cmd}\033[0m", file=log)
            else:
                print(f"\033[031m[Error {ftime} {htime(time.mktime(end_time)-time.mktime(start_time))}]: {cmd}\033[0m", file=log)
        else:
            print(f"\033[032m[Skip {ftime}]: {cmd}\033[0m", file=log)
    return None
    
def rm(file):
    if os.path.isdir(file):
        shutil.rmtree(file, ignore_errors=True)
    else:
        try:
            os.remove(file)
        except FileNotFoundError:
            pass
    return None

def mkdir(path):
    try:
        os.makedirs(path)
    except FileExistsError:
        pass
    return None
    
def run_miniprot(query_file, genome_file, output, thread, mask, skip_align, outs):
    if mask:
        cmd = f"{os.path.join(sys.path[0], 'script', 'sm2rmForFasta.py')} -i {genome_file} -o {genome_file}.tmp"
        #subprocess.run(cmd, shell=True)
        run_cmd(cmd, jobname=None, capture_output=debug_info)
        genome_file = os.path.realpath(genome_file+'.tmp')
    
    #miniprot  = os.path.join(sys.path[0], 'bin', 'miniprot')
    miniprot = miniprot_PATH
    cmd = f'{miniprot} -I --outs={outs} -t {thread} --aln --gff {genome_file} {query_file} > {output}'
    if skip_align:
        pass
    else:
        #subprocess.run(cmd, shell=True)
        run_cmd(cmd, jobname=None, capture_output=debug_info)
    
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
    
miniprot_PATH, TransDecoder_PATH = check_dependencies(ORFSoftware, miniprot_PATH, TransDecoder_PATH)

if os.path.dirname(prefix) != '':
    mkdir(os.path.realpath(os.path.dirname(prefix)))
    
miniprot_output = prefix + '.miniprot_output.gff3'
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

out = open(f'{prefix}.transcript.gtf', 'w')
for index, transcript in enumerate(clusters):
    print(transcript[0], 'quickprot', 'transcript', transcript[1], transcript[2], '.', transcript[3], '.',
          'gene_id "QUKPGENE{}"; transcript_id "QUKPMRNA{}";'.format(index+1, index+1), sep='\t', file=out)
    for exon in clusters[transcript]:
        print(exon[0], 'quickprot', 'exon', exon[1], exon[2], '.', exon[3], '.',
              'gene_id "QUKPGENE{}"; transcript_id "QUKPMRNA{}";'.format(index+1, index+1), sep='\t', file=out)
out.close()

cmd = f"{os.path.join(TransDecoder_PATH, 'util/gtf_to_alignment_gff3.pl')} {prefix}.transcript.gtf > {prefix}.transcript.gff3"
#subprocess.run(cmd, shell=True, capture_output=False)
run_cmd(cmd, jobname=None, capture_output=debug_info)
#cmd = f"{os.path.join(TransDecoder_PATH, 'util/gtf_genome_to_cdna_fasta.pl')} {prefix}.transcript.gtf {genome_file} > {prefix}.transcript.fasta"
cmd = f"{os.path.join(sys.path[0], 'script/gtf_genome_to_cdna_fasta.py')} {prefix}.transcript.gtf {genome_file} > {prefix}.transcript.fasta"
#subprocess.run(cmd, shell=True, capture_output=True)
run_cmd(cmd, jobname=None, capture_output=debug_info)

if os.path.dirname(prefix) == '':
    output_dir = os.getcwd()
else:
    output_dir = os.path.realpath(os.path.dirname(prefix))

if ORFSoftware=="TransDecoder":
    cmd = f"{os.path.join(TransDecoder_PATH, 'TransDecoder.LongOrfs')} -t {prefix}.transcript.fasta --genetic_code {NCBI2TransDecoder_genetic_code[genetic_code]} --output_dir {output_dir}"
    #subprocess.run(cmd, shell=True, capture_output=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    cmd = f"{os.path.join(TransDecoder_PATH, 'TransDecoder.Predict')} -t {prefix}.transcript.fasta --genetic_code {NCBI2TransDecoder_genetic_code[genetic_code]} --output_dir {output_dir}"
    if single_best_only == True:
        cmd += ' --single_best_only'
    #subprocess.run(cmd, shell=True, capture_output=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    cmd = f"{os.path.join(TransDecoder_PATH, 'util/cdna_alignment_orf_to_genome_orf.pl')} \
    {prefix}.transcript.fasta.transdecoder.gff3 {prefix}.transcript.gff3 {prefix}.transcript.fasta > {prefix}.transcript.genome.gff3"
    #subprocess.run(cmd, shell=True, capture_output=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    cmd = f"{os.path.join(sys.path[0], 'script/split_and_filter_gene_model.py')} \
    -i {prefix}.transcript.genome.gff3 -o {prefix}.gff3 --overlap {overlap}"
    #subprocess.run(cmd, shell=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)

    #cmd = f"{os.path.join(TransDecoder_PATH, 'util/gff3_file_to_proteins.pl')} \
    #--gff3 {prefix}.gff3 --fasta {genome_file} --genetic_code {NCBI2TransDecoder_genetic_code[genetic_code]} --seqType prot > {prefix}.pep.fasta"
    cmd = f"{os.path.join(sys.path[0], 'script/extract_sequence_from_gff3.py')} {prefix}.gff3 {genome_file} -G {genetic_code} --seqtype prot -o {prefix}.pep.fasta"
    #subprocess.run(cmd, shell=True, capture_output=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    #cmd = f"{os.path.join(TransDecoder_PATH, 'util/gff3_file_to_proteins.pl')} \
    #--gff3 {prefix}.gff3 --fasta {genome_file} --genetic_code {NCBI2TransDecoder_genetic_code[genetic_code]} --seqType CDS > {prefix}.cds.fasta"
    cmd = f"{os.path.join(sys.path[0], 'script/extract_sequence_from_gff3.py')} {prefix}.gff3 {genome_file} -G {genetic_code} --seqtype CDS -o {prefix}.cds.fasta"
    
    #subprocess.run(cmd, shell=True, capture_output=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    cmd = f"{os.path.join(sys.path[0], 'script/get_longest_transcript_gff3.py')} {prefix}.gff3 -o {prefix}.longest.gff3"
    #subprocess.run(cmd, shell=True, capture_output=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    #cmd = f"{os.path.join(TransDecoder_PATH, 'util/gff3_file_to_proteins.pl')} \
    #--gff3 {prefix}.longest.gff3 --fasta {genome_file} --genetic_code {NCBI2TransDecoder_genetic_code[genetic_code]} --seqType prot > {prefix}.longest.pep.fasta"
    cmd = f"{os.path.join(sys.path[0], 'script/extract_sequence_from_gff3.py')} {prefix}.gff3 {genome_file} -G {genetic_code} --seqtype prot -o {prefix}.longest.pep.fasta"
    #subprocess.run(cmd, shell=True, capture_output=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    #cmd = f"{os.path.join(TransDecoder_PATH, 'util/gff3_file_to_proteins.pl')} \
    #--gff3 {prefix}.longest.gff3 --fasta {genome_file} --genetic_code {NCBI2TransDecoder_genetic_code[genetic_code]} --seqType CDS > {prefix}.longest.cds.fasta"
    cmd = f"{os.path.join(sys.path[0], 'script/extract_sequence_from_gff3.py')} {prefix}.gff3 {genome_file} -G {genetic_code} --seqtype CDS -o {prefix}.longest.cds.fasta"
    #subprocess.run(cmd, shell=True, capture_output=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    intermediate_files = [#os.path.realpath(f'{prefix}.miniprot_output.gff3'),
                          os.path.realpath(f'{prefix}.transcript.gtf'),
                          os.path.realpath(f'{prefix}.transcript.gff3'),
                          os.path.realpath(f'{prefix}.transcript.genome.gff3'),
                          *tuple(glob.glob(os.path.realpath(f'{prefix}.transcript.fasta')+'*'))]
    if noclean:
        pass
    else:
        for file in intermediate_files:
            #print(file)
            rm(file)
            
elif ORFSoftware=="TD2":
    cmd = f"TD2.LongOrfs --precise -t {prefix}.transcript.fasta -G {genetic_code} -O {output_dir} -@ {thread}"
    #subprocess.run(cmd, shell=True, capture_output=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    cmd = f"TD2.Predict --precise -t {prefix}.transcript.fasta -G {genetic_code} -O {output_dir}"
    if single_best_only == True:
        pass
    else:
        #pass
        cmd += ' --all-good'
    #subprocess.run(cmd, shell=True, capture_output=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    cmd = f"{os.path.join(TransDecoder_PATH, 'util/cdna_alignment_orf_to_genome_orf.pl')} \
    {prefix}.transcript.fasta.TD2.gff3 {prefix}.transcript.gff3 {prefix}.transcript.fasta > {prefix}.transcript.genome.gff3"
    #subprocess.run(cmd, shell=True, capture_output=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    cmd = f"{os.path.join(sys.path[0], 'script/split_and_filter_gene_model.py')} \
    -i {prefix}.transcript.genome.gff3 -o {prefix}.gff3 --overlap {overlap}"
    #subprocess.run(cmd, shell=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    #cmd = f"{os.path.join(TransDecoder_PATH, 'util/gff3_file_to_proteins.pl')} \
    #--gff3 {prefix}.gff3 --fasta {genome_file} --genetic_code {NCBI2TransDecoder_genetic_code[genetic_code]} --seqType prot > {prefix}.pep.fasta"
    #subprocess.run(cmd, shell=True, capture_output=True)
    cmd = f"{os.path.join(sys.path[0], 'script/extract_sequence_from_gff3.py')} {prefix}.gff3 {genome_file} -G {genetic_code} --seqtype prot -o {prefix}.pep.fasta"
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    #cmd = f"{os.path.join(TransDecoder_PATH, 'util/gff3_file_to_proteins.pl')} \
    #--gff3 {prefix}.gff3 --fasta {genome_file} --genetic_code {NCBI2TransDecoder_genetic_code[genetic_code]} --seqType CDS > {prefix}.cds.fasta"
    #subprocess.run(cmd, shell=True, capture_output=True)
    cmd = f"{os.path.join(sys.path[0], 'script/extract_sequence_from_gff3.py')} {prefix}.gff3 {genome_file} -G {genetic_code} --seqtype CDS -o {prefix}.cds.fasta"
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    cmd = f"{os.path.join(sys.path[0], 'script/get_longest_transcript_gff3.py')} {prefix}.gff3 -o {prefix}.longest.gff3"
    #subprocess.run(cmd, shell=True, capture_output=True)
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    #cmd = f"{os.path.join(TransDecoder_PATH, 'util/gff3_file_to_proteins.pl')} \
    #--gff3 {prefix}.longest.gff3 --fasta {genome_file} --genetic_code {NCBI2TransDecoder_genetic_code[genetic_code]} --seqType prot > {prefix}.longest.pep.fasta"
    #subprocess.run(cmd, shell=True, capture_output=True)
    cmd = f"{os.path.join(sys.path[0], 'script/extract_sequence_from_gff3.py')} {prefix}.gff3 {genome_file} -G {genetic_code} --seqtype prot -o {prefix}.longest.pep.fasta"
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    #cmd = f"{os.path.join(TransDecoder_PATH, 'util/gff3_file_to_proteins.pl')} \
    #--gff3 {prefix}.longest.gff3 --fasta {genome_file} --genetic_code {NCBI2TransDecoder_genetic_code[genetic_code]} --seqType CDS > {prefix}.longest.cds.fasta"
    #subprocess.run(cmd, shell=True, capture_output=True)
    cmd = f"{os.path.join(sys.path[0], 'script/extract_sequence_from_gff3.py')} {prefix}.gff3 {genome_file} -G {genetic_code} --seqtype CDS -o {prefix}.longest.cds.fasta"
    run_cmd(cmd, jobname=None, capture_output=debug_info)
    
    intermediate_files = [#os.path.realpath(f'{prefix}.miniprot_output.gff3'),
                          os.path.realpath(f'{prefix}.transcript.gtf'),
                          os.path.realpath(f'{prefix}.transcript.gff3'),
                          os.path.realpath(f'{prefix}.transcript.genome.gff3'),
                          os.path.realpath(os.path.join(os.path.dirname(prefix),'longest_orfs.cds')),
                          os.path.realpath(os.path.join(os.path.dirname(prefix),'longest_orfs.pep')),
                          os.path.realpath(os.path.join(os.path.dirname(prefix),'longest_orfs.gff3')),
                          os.path.realpath(os.path.join(os.path.dirname(prefix),'psauron_score.csv')),
                          *tuple(glob.glob(os.path.realpath(f'{prefix}.transcript.fasta')+'*'))]
    #print(intermediate_files)
    
    if noclean:
        pass
    else:
        for file in intermediate_files:
            rm(file)
