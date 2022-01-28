# GNFish
Suite of small Python programs following the pipeline detailed in (DOI of article yet to be published) for:<br />
(1) “fishing” specific gene sequences of interest from genomic data available in NCBI databases<br />
(2) processing and depuration of retrieved sequences<br />
(3) production of a multiple sequence alignment<br />
(4) selection of best-fit model of evolution<br />
(5) solid reconstruction of a phylogenetic tree.<br />

There are 12 scripts:<br />

**align_sequences.py** -> Align sequences using MAFFT software<br />
**blast.py** -> Perfoms BLAST searches against the download genomes<br />
**class_list_files.py** -> Class for accessing directories. Not need to be run, but the other scripts need it.<br />
**decompress_genomes.py** -> Decompress genome files.<br />
**get_combined_seqs.py** -> Merge sequences for creating a whole alingment.<br />
**get_genomes.py** -> Download genomes from NCBI databases<br />
**get_protein_query_.py** -> Download a dataset of protein sequences that act as query for BLAST searches<br />
**get_RAW_sequences.py** -> Extracts the sequences mapped after BLAST searches<br />
**get_unique_hits.py** -> Gets unique hits from Blast output files based on genomes IDs<br />
**iqtree.py** -> Runs IQ-Tree program for phylogentic inference<br />
**translate_seq.py** -> Translate nucleotide sequences to protein.<br />
**trimal.py** -> Runs trimAl program for trimming alignments.<br />

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
It needs your e-mail and a file with your queries.<br />
Queries are usallly taxa names (see query_genomes examples files at ./Examples)<br />
By default it will download genomic data, but you can also add protein or RNA, using --protein or --rna respectively.<br />
Retmax set a maximum number of downloaded genomes for query.<br />
It creates a directory named Data and three subdirectories named Genomic, Rna and Protein where it downloads the genomes.<br />
Genomes file are compressed you must descompress for working with them (see decompress genome file).<br />
Genome file name will be "Genus_species_assemblyID_datatype.f[n,a]a".


Type on terminal get_genomes.py -h for further information.<br />

### USAGE:<br />

**get_genomes.py 'e-mail' 'query.txt' **<br />

### PARAMETERS:<br />

**--email** -> mandatory e-mail for NCBI searches<br />
**--query** -> file with the queries. Usually simple taxa names (species, group). Field tags or filters can be added to each query. See examples below or look at Examples directory for examples of query files<br />

Optional parameters:<br />
**--genomic** -> downloads whole genomic data<br />
**--rna** -> downloads protein annotation data<br />
**--protein** -> downloads protein annotation data<br />
**--exclusive** -> downloads just protein or rna annotation data if available. No genomic backup<br />
**--retmax** -> number of NCBI records reported for every query. Default value equal to 200<br />
**--refine** -> adds filter or field information to all queries. Constant value AND (latest[filter] AND "representative genome"[filter] AND all[filter] NOT anomalous[filter])<br />

### EXAMPLES:

Examples files are availabe at Example directory.

#### 1.- Plain search
Query txt with simple searches, just taxa nor filters or field tags. No refine argument.<br />
See ./Examples/query_genome_3.txt for examples of queries.<br />
Creates ./Data directory and ./Data/Genomic, ./Data/Rna and ./Data/Protein subdirectories.<br />
Nor filters or field tags applied. Not curated and redundant genomes (one genom for more than one species is used).<br />
Use this when you do not care very much about filtering.<br />

Genomic donwload
```
./Code/get_genomes.py hlorente@ucm.es Example/query.txt
```
Protein donwload, look for genomic data as backup (For transcrits use --rna)
```
./Code/get_genomes.py hlorente@ucm.es Example/query.txt --protein
```
Exclusive protein download. No backup. (For transcrits use --rna)
```
./Code/get_genomes.py hlorente@ucm.es Example/query.txt --protein --exlusive
```
Protein and Rna downloading. Genomic search as backup
```
./Code/get_genomes.py hlorente@ucm.es Example/query.txt --protein --rna
```
Rna donwload changing number of records downloaded. Will donwload until 1000 genomes if available. By default 200 are donwloaded.
```
./Code/get_genomes.py hlorente@ucm.es Example/query.txt --rna --retmax 1000
```
#### 2.- Refine search using --refine argument. **Recomended**
Default. Applies Representative (just one genome for species), Latest, Not Anomalous.<br />
Just protein example. For rna, genomic, exclusive or retmax argument see 1 above and them to this command line.<br />
```
./Code/get_genomes.py hlorente@ucm.es Example/query.txt --protein --refine
```
Own filters or field tags. 
```
./Code/get_genomes.py hlorente@ucm.es Example/query.txt --protein --refine ''
```
More info about filters and field tags at x and reading y.

#### 3.- Refine applied to each query
You must add a the field tags or filters after your query.<br />
See ./Examples/query_genome_3.txt for examples of queries.<br />
As each contains is own filters or field tags this will be applied exclusively.<br />
Just protein example. For rna, genomic, exclusive or retmax argument see 1 above and them to this command line.<br />
```
./Code/get_genomes.py hlorente@ucm.es Example/query_filters.txt --protein
```

## get_query_seqs.py<br />

### DESCRIPTION:<br />

Program for dowloading protein or nucleotide from NCBI Databases through Entrez. Requires internet conection.<br />
It needs your e-mail and a file with your queries.<br />
Queries must follow "Gene name (filters). If no filters type () after the gene name.<br />
By default it will download protein data, but you can change using --nucleotide parameter.<br />
Retmax set a maximum number of downloaded genomes for query.<br />
Curated refine the search or protein, just those entries which include the name of the protein in Protein Feature.
Stores sequences at ./Data/Query_seqs.<br />
Output file will be named "genename_query_seqs.fas.<br />

Type on terminal get_genomes.py -h for further information.<br />

### USAGE:<br />

**get_query_seqs_.py 'e-mail' 'query.txt' **<br />

### PARAMETERS:<br />

**--email** -> mandatory e-mail for NCBI searches<br />
**--query** -> file with the queries. Usually simple taxa names (species, group). Field tags or filters can be added to each query. See examples below or look at Examples directory for examples of query files<br />

Optional parameters:<br />
**--nucleotide** -> for nucleotide downloading<br />
**--curated** -> downloads just protein or rna annotation data if available. No genomic backup<br />
**--retmax** -> number of NCBI records reported for every query. Default value equal to 200<br />
**--refine** -> adds filter or field tags to the query. Constant value refseq[filter]. Follow constant value structure for your custom filters, begin with "AND"<br />

### EXAMPLES:

Examples files are availabe at Example directory.

#### 1.- Simple search
Query txt with simple searches, just gene name and filters equal to ().<br />
Creates ./Data/Query_seqs directory.<br />
Peforms well for protein for nucleotide it is better to see 2.

Protein download
```
./Code/get_query_seqs_.py hlorente@ucm.es Example/query_seqs.txt
```
Protein download, restricting maximum number to 50
```
./Code/get_query_seqs_.py hlorente@ucm.es Example/query_seqs.txt --retmax 50
```
Protein download, restricting maximum number to 50 and curating
```
./Code/get_query_seqs_.py hlorente@ucm.es Example/query_seqs.txt --retmax 50 --curated
```
Nucleotide downloading, restricting maximum number to 50. Better see 2.
```
./Code/get_query_seqs_.py hlorente@ucm.es Example/query_seqs.txt --retmax 50 --nucleotide
```

#### 2.- Refine search using --refine argument.
Default. Applies Representative (just one genome for species), Latest, Not Anomalous.<br />
Just protein example. For rna, genomic, exclusive or retmax argument see 1 above and them to this command line.<br />
```
./Code/get_genomes.py hlorente@ucm.es Example/query.txt --protein --refine
```
Own filters or field tags. 
```
./Code/get_genomes.py hlorente@ucm.es Example/query.txt --protein --refine ''
```
More info about filters and field tags at x and reading y.

#### 3.- Refine applied to each query
As each contains is own filters or field tags this will be applied exclusively. See ./Example/query_filters.txt for a scheme.<br />
Just protein example. For rna, genomic, exclusive or retmax argument see 1 above and them to this command line.<br />
```
./Code/get_genomes.py hlorente@ucm.es Example/query_filters.txt --protein
```
