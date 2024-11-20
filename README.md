# Quick genome annotation base on protein (quickprot)


## Installation

Before use, you need to install Python and Perl.

Python3 >= 3.8, perl >= 5

In running quickprot, protein alignment is done using [miniprot](https://github.com/lh3/miniprot/) (v0.12), and ORFs prediction is done using [TransDecoder](https://github.com/TransDecoder/TransDecoder) (v5.7.1). For ease of use, these two software are integrated into quickprot.

```
wget https://github.com/thecgs/quickprot/archive/refs/tags/quickprot-v1.11.tar.gz
tar -zxvf quickprot-v1.11.tar.gz
cd quickprot-v1.11
./quickprot -h
```

Note

```
# if you used --mask optional of qucikprot.py script, you has install biopython
pip install biopython

# if you used sort_gff3.py script has install natsort.
pip install natsort
```

## Algorithm

![Schema of quickprot algorithm](./docs/Schema_of_quickprot_algorithm.png#pic_center)

Schema of quickprot algorithm


## Usage

To run quickprot, use

```
./quickprot -q protein.fasta -g genome.fasta
```

quickprot optional

```
./quickprot.py -h
usage: quickprot.py -q str -g str [-p str] [-i float] [--outs float] [--overlap float] [-t int] [-G str] [-s] [-m] [-n] [-b] [-h] [-v]

Quick genome annotation base on protein.

required arguments:
  -q str, --query str   A file of query protein fasta format.
  -g str, --genome str  A file of genome fasta format.

optional arguments:
  -p str, --prefix str  Prefix of a output file. default=quickprot
  -i float, --identity float
                        Alignment identity (0-1). default=0.95
  --outs float          Output score at least bestScore (0-1). default=0.99
  --overlap float       If the overlap of predicted ORFs in a transcript is less than default value (0-1). default=0.8, 
                        they will be dissected.
  -t int, --thread int  Thread number of run miniprot sortware. defualt=24
  -G str, --genetic_code str
                        Genetic Codes (derived from: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). defualt=Universal
                        The supported genetic codon tables are Acetabularia, Candida, Ciliate, Dasycladacean, Euplotid, Hexamita,
                        Mesodinium, Mitochondrial-Ascidian, Mitochondrial-Chlorophycean, Mitochondrial-Echinoderm, Mitochondrial-Flatworm,
                        Mitochondrial-Invertebrates, Mitochondrial-Protozoan, Mitochondrial-Pterobranchia, Mitochondrial-Scenedesmus_obliquus,
                        Mitochondrial-Thraustochytrium, Mitochondrial-Trematode, Mitochondrial-Vertebrates, Mitochondrial-Yeast, 
                        Pachysolen_tannophilus,Peritrich, SR1_Gracilibacteria, Tetrahymena, and Universal.
  -s, --skip_align      Skip run miniprot step. default=False
  -m, --mask            Soft-masked (dna_sm) genome convert to masked(dna_rm) genome. default=False
  -n, --noclean         Do not delete intermediate files. default=False
  -b, --single_best_only
                        Retain only the single best orf per transcript. default=False
                        It is not recommended to use it because when two reference proteins overlap during alignment, 
                        it can lead to fusion during transcript assembly. If a transcript is not set with only one ORF,
                        the fused ORF will be split in subsequent analysis.
  -h, --help            Show program's help message and exit.
  -v, --version         Show program's version number and exit.

date:2024/11/19 author:guisen chen email:thecgs001@foxmail.com
```
