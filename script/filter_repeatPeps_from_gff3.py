#!/usr/bin/env python
# coding: utf-8

import re
import os
import sys
import argparse
from collections import defaultdict

def run_cmd(cmd, jobname=None):
    
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
    
    start_time = time.gmtime()
    ftime = time.strftime('%Y-%m-%d %H:%M:%S', start_time)
    print(f"\033[033m[Running {ftime}]: {cmd}\033[0m", file=log)
    if jobname==None:
        res = subprocess.run(cmd, shell=True)
        end_time = time.gmtime()
        ftime = time.strftime('%Y-%m-%d %H:%M:%S', end_time)
        if res.returncode == 0:
            print(f"\033[032m[Finished {ftime} {htime(time.mktime(end_time)-time.mktime(start_time))}]: {cmd}\033[0m", file=log)
        else:
            print(f"\033[031m[Error {ftime} {htime(time.mktime(end_time)-time.mktime(start_time))}]: {cmd}\033[0m", file=log)
    else:
        if not os.path.exists(f'{jobname}.done'):
            res = subprocess.run(cmd, shell=True)
            end_time = time.gmtime()
            if res.returncode == 0:
                out = open(f'{jobname}.done', 'w')
                out.close()
                print(f"\033[032m[Finished {ftime} {htime(time.mktime(end_time)-time.mktime(start_time))}]: {cmd}\033[0m", file=log)
            else:
                print(f"\033[031m[Error {ftime} {htime(time.mktime(end_time)-time.mktime(start_time))}]: {cmd}\033[0m", file=log)
        else:
            print(f"\033[032m[Skip {ftime}]: {cmd}\033[0m", file=log)
    return None

def run_diamond(prefix, db, query, thread=os.cpu_count(), evalue='1e-5', model=2, subcmd='blastp', diamond_path='diamond', identity=60, query_cover=60):
    model_dict = {1: '--faster', 
                  2: '--fast',
                  3: '--mid-sensitive',
                  4: '--sensitive',
                  5: '--more-sensitive',
                  6: '--very-sensitive',
                  7: '--ultra-sensitive'
                 }
    
    cmd = f"diamond makedb --threads {thread} --db RepeatPeps --in {db}"
    run_cmd(cmd, jobname=None)
    
    cmd = f"{diamond_path} {subcmd} {model_dict.get(model, '--sensitive')} --query-cover {query_cover} --id {identity} --threads {thread} --outfmt 6 --quiet --no-parse-seqids --db RepeatPeps.dmnd --query {query} --evalue {evalue} --out {prefix}.blastout.tsv"
    run_cmd(cmd, jobname=None)
    return None

def filter_gff3(prefix, gff3_file, outfile1=None, outfile2=None):
    filter_mRNA_IDs = set()
    with open(f"{prefix}.blastout.tsv", "r") as f:
        for l in f:
            filter_mRNA_IDs.add(l.strip().split('\t')[0])
            
    genes = defaultdict(list)
    mRNAs = defaultdict(list)
    gene2mRNA = defaultdict(list)
    features = defaultdict(list)
    
    with open(gff3_file, 'r')  as f:
        for l in f:
            if not l.startswith('#') and l.strip() != '':
                l = l.split('\t')
                if l[2] == "gene":
                    geneID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
                    genes[geneID] = l
                    
                elif bool(re.search('RNA', l[2], flags=re.I)):   #include mRNA, tRNA, miRNA ...
                    mRNAID = re.search('ID=(.*?)[;,\n]', l[8]).group(1)
                    geneID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                    gene2mRNA[geneID].append(mRNAID)
                    mRNAs[mRNAID] = l
                    
                elif l[2] == "CDS":
                    mRNAID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                    features[mRNAID].append(l)
                    
                else:
                    mRNAID = re.search('Parent=(.*?)[;,\n]', l[8]).group(1)
                    features[mRNAID].append(l)
                    
    if outfile1 == None:
        out1 = sys.stdout
    else:
        out1 = open(outfile1, 'w')
        
    if outfile2 == None:
        out2 = sys.stdout
    else:
        out2 = open(outfile2, 'w')
    
    for gene in gene2mRNA:
        status = True 
        for mRNA in gene2mRNA[gene]:
            if mRNA in filter_mRNA_IDs:
                status = False
        if status:
            print(*tuple(genes[gene]), sep='\t', end='', file=out1)
            for mRNA in gene2mRNA[gene]:
                print(*tuple(mRNAs[mRNA]), sep='\t', end='', file=out1)
                for feature in features[mRNA]:
                    print(*tuple(feature), sep='\t', end='', file=out1)
            print('', file=out1)
        else:
            print(*tuple(genes[gene]), sep='\t', end='', file=out2)
            for mRNA in gene2mRNA[gene]:
                print(*tuple(mRNAs[mRNA]), sep='\t', end='', file=out2)
                for feature in features[mRNA]:
                    print(*tuple(feature), sep='\t', end='', file=out2)
            
            print('', file=out2)
    out1.close()
    out2.close()
    return None

prefix = "Query2RepeatPeps"
db = os.path.join(os.path.dirname(sys.path[0]), "lib", "RepeatPeps.lib")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter repeat protein from gff3 file.', add_help=False,
                                     epilog='Date:2025/09/06 Author:Guisen Chen Email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument('-q', '--query',  metavar='str', required=True,
                          help='A file of fasta format.')
    required.add_argument('-g', '--gff3', metavar='str', required=True, 
                          help='A input file of gff3 format.')
    optional.add_argument('-o1', '--output1', metavar='str', default="retain.gff3",
                          help=f'A retain output file of gff3 format. defualt=retain.gff3')
    optional.add_argument('-o2', '--output2', metavar='str', default="discard.gff3",  
                          help=f'A discard output file of gff3 format. defualt=discard.gff3')
    optional.add_argument('-d', '--db', metavar='str', default=db,
                          help=f'A path of repeat protein database. default={db}')
    optional.add_argument('-c', '--subcmd', metavar='str', default='blastp', choices=['blastp', 'blastx'],
                          help='A subcommand of diamond software. defualt=blastp')
    optional.add_argument('-e', '--evalue', metavar='str', default='1e-10',
                          help='A maximum evalue of report alignments. defualt=1e-10')
    optional.add_argument('-m', '--model',  metavar='int', default=2, choices=range(1,8), type=int,
                          help='A running model of diamond software, 1: faster; 2: fast; 3: mid-sensitive; 4: sensitive; 5: more-sensitive; 6: very-sensitive; 7: ultra-sensitive. default=2')
    optional.add_argument('-t', '--thread', metavar='int', default=os.cpu_count(), type=int,
                          help='A thread number of diamond software. defualt='+str(os.cpu_count()))
    optional.add_argument('--diamond_path', metavar='str', default='diamond',
                          help='A path of diamond software. defualt=diamond')                   
    optional.add_argument('--identity', metavar='int', default=60,
                          help='minimum identity (percentage) to report an alignment. defualt=60')
    optional.add_argument('--query_cover', metavar='int', default=60,
                          help='minimum query cover (percentage) to report an alignment. defualt=60')
    optional.add_argument('-h', '--help', action='help', 
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.00',  
                          help="Show program's version number and exit.")
    
    args = parser.parse_args()
    run_diamond(prefix=prefix,
                db=args.db,
                query=args.query,
                thread=args.thread,
                evalue=args.evalue,
                model=args.model,
                subcmd=args.subcmd,
                diamond_path=args.diamond_path,
                identity=args.identity,
                query_cover=args.query_cover)
    filter_gff3(prefix=prefix, gff3_file=args.gff3, outfile1=args.output1, outfile2=args.output2)
