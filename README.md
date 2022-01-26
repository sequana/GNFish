# GNFish
Suite of small Python programs following the pipeline detailed in (DOI of article yet to be published) for:
(1) “fishing” specific gene sequences of interest from genomic data available in NCBI databases
(2) processing and depuration of retrieved sequences
(3) production of a multiple sequence alignment
(4) selection of best-fit model of evolution
(5) solid reconstruction of a phylogenetic tree.

There are 10 scripts:

align_sequences.py -> Align sequences using MAFFT software
blast.py -> Perfoms BLAST searches against the download genomes
class_list_files.py -> Class for accessing directories. Not need to be run, but the other scripts need it. 
get_genomes.py -> Download genomes from NCBI databases
get_protein_query_.py -> Download a dataset of protein sequences that act as query for BLAST searches
get_RAW_sequences.py -> Extracts the sequences mapped after BLAST searches
get_unique_hits.py -> Gets unique hits from Blast output files based on genomes IDs
iqtree.py -> 
translate_seq.py -> 
trimal.py -> 

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

This software uses Biopython module in python 3 (tested version 1.78-2). They can be installed in your system using pip:

```
pip install biopython 
```

**DESCRIPTION:**


**This software is released under the license GNU GPLv3.**

**This software is provided as is without warranty of any kind.** 

People wishing to contribute to the software, report issues or seek support can contact Hector Lorente at hlorente@ucm.es
