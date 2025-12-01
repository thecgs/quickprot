

# QuickProt User Guide

## Update

- 2025/09/16
  1. Added add_type_gff3.py script.
  2. Optimized prediction of stop codon.
  3. Four new options have been added, namely "-c", "-ps", "-ms", and "-an", for quality control of protein mapping, rational use to reduce pseudogenes.

- 2025/09/08
  1.  Added gtf_genome_to_cdna_fasta.py and gtf_genome_to_cdna_fasta.py script.
  2. Input files now support .gz compressed files.
  3. Now QuickProt can output more running details.

- 2025/09/06
  1. Provides information about the running process.
  2. Added filter_repeatPeps_from_gff3.py script for removing repeat proteins from gff3 file.
- 2025/05/28
  1. Provide -ORFSoftware TD2 option, you can use TD2 as a tool for ORF prediction.
  2. Optimization of genetic code options.

##  What is QuickProt？

The QuickProt algorithm is a homology-based method for predicting gene models across entire genomes, designed to rapidly construct a non-redundant set of gene models. As illustrated in Figure 1, its core principle is analogous to the blotting method. It primarily employs [miniprot](https://github.com/lh3/miniprot/) (v0.18), to align homologous protein sequences to the genome, delineates high-alignment regions to assemble pseudo-transcripts (lacking UTR regions), and predicts coding regions within these pseudo-transcripts using [TransDecoder](https://github.com/TransDecoder/TransDecoder) (v5.7.1). Subsequently, low-quality gene models are filtered out and chimeric gene models are dissected, ultimately generating a high-accuracy, non-redundant gene set.

![Schema of quickprot algorithm](./docs/Schema_of_quickprot_algorithm.png#pic_center)

<center>Fig1. Schema of QuickProt algorithm</center>

## Installation:

Before use, you need to install Perl, Python, and biopython.

Python3 >= 3.8, perl >= 5 

For ease of use, [miniprot](https://github.com/lh3/miniprot/) (v0.18) and  [TransDecoder](https://github.com/TransDecoder/TransDecoder) (v5.7.1)  software are integrated into QuickProt.

```
wget https://github.com/thecgs/quickprot/archive/refs/tags/quickprot-v1.8.0.tar.gz
tar -zxvf quickprot-v1.8.0.tar.gz
cd quickprot-v1.8.0
./quickprot -h
```

Note:

```
# if you need to use --mask optional of qucikprot.py script, and you need to install biopython
pip install biopython

# if you need to use sort_gff3.py script, and you need to install natsort.
pip install natsort

# if you need to use -ORFSoftware TD2, and you need to install TD2
pip install TD2
```

## Usage:

To  quickly  run QuickProt software. like this, 

```
./quickprot -q protein.fasta -g genome.fasta
```

This pipeline can improve busco missing result, but you need to download [compleasm](https://github.com/huangnengCSU/compleasm) software.

```
## step1. running quickprot software
./quickprot.py -q protein.fasta -g genome.fasta -p quickprot.raw

## step2. running compleasm software
compleasm.py run -a genome.fasta -o ./ -l your_lineage

## step3. to update raw gff3 of step1 from compleasm result
./script/update_gff3_from_minibusco.py -r quickprot.raw.longest.gff3 -m ./your_lineage/miniprot_output.gff -g genome.fasta -o improve_busco.gff3

## step4. merge step1 and step3 gff3 result
cat quickprot.raw.longest.gff3 improve_busco.gff3 > genome.longest.gff.tmp

## step5. to sort by chromosomes or scaffold and gene start position and to rename gff3
./script/sort_gff3.py genome.longest.gff.tmp | ./script/rename_gff3.py - -o genome.longest.gff3 -p QUICKPROT; rm genome.longest.gff.tmp

## step6. extract protein sequence from genome and gff file
./bin/TransDecoder-5.7.1/util/gff3_file_to_proteins.pl --gff3 genome.longest.gff3 --fasta genome.fasta --seqType prot  > genome.longest.pep.fasta

## step7. extract CDS sequence from genome and gff file
./bin/TransDecoder-5.7.1/util/gff3_file_to_proteins.pl --gff3 genome.longest.gff3 --fasta genome.fasta --seqType CDS  > genome.longest.cds.fasta
```

This step can help you remove repeat proteins (e.g. ENV, Gag, Pol, RT, RH, INT, etc.), but you need to download [diamond](https://github.com/bbuchfink/diamond).

```
./script/filter_repeatPeps_from_gff3.py -q genome.longest.pep.fasta -g genome.longest.gff3

## results
## retain.gff3 —— Gene model without repeat proteins
## discard.gff3 —— Gene model of repeat proteins
```

## Run with Singularity

Download the Singularity image [here](https://zenodo.org/records/17718142/files/quickprot.v1.8.sif?download=1)

```
singularity exec -B PATH -e quickprot.v1.8.sif quickprot.py -q protein.fasta -g genome.fasta
```

## Cite QuickProt:

If you use QuickProt, please cite:

> Guisen Chen, Hehe Du, Zhenjie Cao, Ying Wu, Chen Zhang, Yongcan Zhou, Jingqun Ao, Yun Sun, Zihao Yuan. QuickProt: A Fast and Accurate Homology-Based Protein Annotation Tool for Non-Model Organism Genomes to Advance Comparative Genomics.

