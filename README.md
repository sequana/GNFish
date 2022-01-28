# GNFish
Suite of small Python programs following the pipeline detailed in (DOI of article yet to be published) for:
(1) “fishing” specific gene sequences of interest from genomic data available in NCBI databases
(2) processing and depuration of retrieved sequences
(3) production of a multiple sequence alignment
(4) selection of best-fit model of evolution
(5) solid reconstruction of a phylogenetic tree.

There are 12 scripts:

align_sequences.py -> Align sequences using MAFFT software
blast.py -> Perfoms BLAST searches against the download genomes
class_list_files.py -> Class for accessing directories. Not need to be run, but the other scripts need it.
decompress_genomes.py -> Decompress genome files.
get_combined_seqs.py -> Merge sequences for creating a whole alingment.
get_genomes.py -> Download genomes from NCBI databases
get_protein_query_.py -> Download a dataset of protein sequences that act as query for BLAST searches
get_RAW_sequences.py -> Extracts the sequences mapped after BLAST searches
get_unique_hits.py -> Gets unique hits from Blast output files based on genomes IDs
iqtree.py -> Runs IQ-Tree program for phylogentic inference
translate_seq.py -> Translate nucleotide sequences to protein.
trimal.py -> Runs trimAl program for trimming alignments.

**This software is released under the license GNU GPLv3.**

**This software is provided as is without warranty of any kind.** 

People wishing to contribute to the software, report issues or seek support can contact Hector Lorente Martinez at hlorente@ucm.es


# Download

Download scripts from https://github.com/hectorloma/GNFish

In linux terminal

```
git clone https://github.com/hectorloma/GNFish
```

If git is not installed visit https://github.com/git-guides/install-git  
which is available for any OS

Make scripts executable if you want to run it from the very same terminal

```
cd GNFish
```

```
chmod +x *.py
```

# Dependencies

This software uses Biopython module in python 3 (tested version 1.78-2). It can be installed in your system using pip:

```
pip install biopython 
```

# get_genomes.py

**DESCRIPTION:**

Program for dowloading genomes from NCBI Databases through Entrez. Requires internet conection.
It needs your e-mail and a file with your queries. By default it will download genomic data, but you can also add protein or RNA, using --protein or --rna respectively.
It creates a directory named Data and three subdirectories named Genomic, Rna and Protein where it downloads the genomes.
Genomes file are compressed you must descompress for working with them (see decompress genome file).

Type on terminal get_genomes.py -h for further information.

**USAGE:**


**get_genomes.py 'e-mail' './Example/query.txt' -c convertfile_Y/N > assembly.fas**


**PARAMETERS:**

**--email** -> mandatory e-mail for NCBI searches
	parser.add_argument('query', help='file with the queries', type=argparse.FileType('r'))
	parser.add_argument('--genome', help='downloads whole genomic data', action='store_true')
	parser.add_argument('--rna', help='download protein annotation data', action='store_true')
	parser.add_argument('--protein', help='download protein annotation data', action='store_true')
	parser.add_argument('--exclusive', help='download just protein or rna annotation data if available.', action='store_true')
	parser.add_argument('--retmax', help='number of NCBI records reported for every query. Default value equal to 200',nargs='?', const=200, type=int, default=200)
	parser.add_argument('--filters', help='add

**-m MITOFILE --mitofile** -> Input file with mitogenome sequence in fasta format, as submitted to MITOS2

**-g GENESFILE, --genesfile** -> Input file with MITOS2 output with individual genes in fasta format

**-c CONVERTFILE, --convertfile** -> Gene names (fasta headers) from MITOS2 will be simplified and made compliant with aln2tbl. Yes=Y No=N

**EXAMPLE:**



./Code/get_genomes.py -hlorente@ucm.es -g ./example/input/Hyalella_solida_genes_mitos2.fas -c Y > ./example/input/Hyalella_solida_assembly.fas


