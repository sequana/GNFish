# GNFish
Suite of small Python programs following the pipeline detailed in (DOI of article yet to be published) for:<br />
(1) “fishing” specific gene sequences of interest from genomic data available in NCBI databases<br />
(2) processing and depuration of retrieved sequences<br />
(3) production of a multiple sequence alignment<br />
(4) selection of best-fit model of evolution<br />
(5) solid reconstruction of a phylogenetic tree.<br />

There are 12 scripts:<br />

align_sequences.py -> Align sequences using MAFFT software<br />
blast.py -> Perfoms BLAST searches against the download genomes<br />
class_list_files.py -> Class for accessing directories. Not need to be run, but the other scripts need it.<br />
decompress_genomes.py -> Decompress genome files.<br />
get_combined_seqs.py -> Merge sequences for creating a whole alingment.<br />
get_genomes.py -> Download genomes from NCBI databases<br />
get_protein_query_.py -> Download a dataset of protein sequences that act as query for BLAST searches<br />
get_RAW_sequences.py -> Extracts the sequences mapped after BLAST searches<br />
get_unique_hits.py -> Gets unique hits from Blast output files based on genomes IDs<br />
iqtree.py -> Runs IQ-Tree program for phylogentic inference<br />
translate_seq.py -> Translate nucleotide sequences to protein.<br />
trimal.py -> Runs trimAl program for trimming alignments.<br />

**This software is released under the license GNU GPLv3.**<br />

**This software is provided as is without warranty of any kind.**<br />

People wishing to contribute to the software, report issues or seek support can contact Hector Lorente Martinez at hlorente@ucm.es<br />


## Download<br />

Download scripts from https://github.com/hectorloma/GNFish<br />
In linux terminal

```
git clone https://github.com/hectorloma/GNFish
```
If git is not installed visit https://github.com/git-guides/install-git which is available for any OS<br />
<br />
Make scripts executable if you want to run it from the very same terminal
```
cd GNFish
```
```
chmod +x *.py
```


## Dependencies<br />

This software uses Biopython module in python 3 (tested version 1.78-2). It can be installed in your system using pip:
```
pip install biopython 
```


## get_genomes.py<br />

### DESCRIPTION:<br />

Program for dowloading genomes from NCBI Databases through Entrez. Requires internet conection.<br />
It needs your e-mail and a file with your queries. By default it will download genomic data, but you can also add protein or RNA, using --protein or --rna respectively.<br />
It creates a directory named Data and three subdirectories named Genomic, Rna and Protein where it downloads the genomes.<br />
Genomes file are compressed you must descompress for working with them (see decompress genome file).<br />

Type on terminal get_genomes.py -h for further information.<br />

### USAGE:<br />

**get_genomes.py 'e-mail' 'query.txt' -c convertfile_Y/N > assembly.fas**<br />

### PARAMETERS:<br />

**--email** -> mandatory e-mail for NCBI searches<br />
**--query**-> file with the queries. Usually simple taxa names (species, group). Field tags or filters can be added to each query. See examples below or look at Examples directory for examples of query files.<br />
Optional parameters
**--genomic** -> downloads whole genomic data<br />
**--rna** -> downloads protein annotation data<br />
**--protein** ->download protein annotation data', action='store_true')<br />
**--exclusive', help='download just protein or rna annotation data if available.', action='store_true')
**--retmax', help='number of NCBI records reported for every query. Default value equal to 200',nargs='?', const=200, type=int, default=200)
**--refine', help='add

**-m MITOFILE --mitofile** -> Input file with mitogenome sequence in fasta format, as submitted to MITOS2

**-g GENESFILE, --genesfile** -> Input file with MITOS2 output with individual genes in fasta format

**-c CONVERTFILE, --convertfile** -> Gene names (fasta headers) from MITOS2 will be simplified and made compliant with aln2tbl. Yes=Y No=N

### EXAMPLES:

Examples files are availabe at Example directory.

#### 1.- Plain search
Query txt with simple searches, just taxa nor filters or field tags. No refine argument.
Creates ./Data directory and ./Data/Genomic, ./Data/Rna and ./Data/Protein subdirectories. 
Nor filters or field tags applied. Not curated and redundant genomes (one genom for more than one species is used)
Use this when you do not care very much about filtering.

Genomic donwload
```
./Code/get_genomes.py -hlorente@ucm.es Example/query.txt
```
Protein donwload, look for genomic data as backup (For transcrits use --rna)
```
./Code/get_genomes.py -hlorente@ucm.es Example/query.txt --protein
```
Exclusive protein download. No backup. (For transcrits use --rna)
```
./Code/get_genomes.py -hlorente@ucm.es Example/query.txt --protein --exlusive
```
Protein and Rna downloading. Genomic search as backup
```
./Code/get_genomes.py -hlorente@ucm.es Example/query.txt --protein --rna
```
Rna donwload changing number of records downloaded. Will donwload until 1000 genomes if available. By default 200 are donwloaded.
```
./Code/get_genomes.py -hlorente@ucm.es Example/query.txt --rna --retmax 1000
```
#### 2.- Refine search using --refine argument. **Recomended**
Default. Applies Representative (just one genome for species), Latest, Not Anomalous.
Just protein example. For rna, genomic, exclusive or retmax argument see 1 above and them to this command line.
```
./Code/get_genomes.py -hlorente@ucm.es Example/query.txt --protein --refine
```
Own filters or field tags. 
```
./Code/get_genomes.py -hlorente@ucm.es Example/query.txt --protein --refine ''
```
More info about filters and field tags at x and reading y.

#### 3.- Refine applied to each query
As each contains is own filters or field tags this will be applied exclusively. See ./Example/query_filters.txt for a scheme.
Just protein example. For rna, genomic, exclusive or retmax argument see 1 above and them to this command line.
```
./Code/get_genomes.py -hlorente@ucm.es Example/query_filters.txt --protein
```
