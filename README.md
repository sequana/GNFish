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
Queries are usallly taxa names (see query_genomes examples files at ./Example)<br />
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
**--query** -> file with the queries. Usually simple taxa names (species, group). Field tags or filters can be added to each query. See examples below or look at Example directory for examples of query files<br />

Optional parameters:<br />
**--genomic** -> downloads whole genomic data<br />
**--rna** -> downloads protein annotation data<br />
**--protein** -> downloads protein annotation data<br />
**--exclusive** -> downloads just protein or rna annotation data if available. No genomic backup<br />
**--retmax** -> number of NCBI records reported for every query. Default value equal to 200<br />
**--refine** -> adds filter or field information to all queries. Constant value AND (latest[filter] AND "representative genome"[filter] AND all[filter] NOT anomalous[filter])<br />

### EXAMPLES:

Examples files are availabe at Example directory.

#### 1.- Simple search
Query txt with simple searches, just taxa nor filters or field tags. No refine argument.<br />
See ./Example/query_genome_2.txt for examples of queries.<br />
Creates ./Data directory and ./Data/Genomic, ./Data/Rna and ./Data/Protein subdirectories.<br />
Nor filters or field tags applied. Not curated and redundant genomes (one genom for more than one species is used).<br />
Use this when you do not care very much about filtering.<br />

Genomic
```
./Code/get_genomes.py hlorente@ucm.es Example/query.txt
```
Protein. Look for genomic data as backup (For transcrits use --rna)
```
./Code/get_genomes.py hlorente@ucm.es Example/query.txt --protein
```
Exclusive protein. No backup. (For transcrits use --rna)
```
./Code/get_genomes.py hlorente@ucm.es Example/query.txt --protein --exlusive
```
Protein and Rna. Genomic search as backup
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
See ./Example/query_genome_2.txt for examples of queries.<br />
As each contains is own filters or field tags this will be applied exclusively.<br />
Just protein example. For rna, genomic, exclusive or retmax argument see 1 above and them to this command line.<br />
```
./Code/get_genomes.py hlorente@ucm.es Example/query_filters.txt --protein
```

#### 4.- Refine to each query and general refine
You must add a the field tags or filters after your query.<br />
See ./Example/query_genome_2.txt for examples of queries.<br />
As each contains is own filters or field tags this will be applied exclusively.<br />
Just protein example. For rna, genomic, exclusive or retmax argument see 1 above and them to this command line.<br />
```
./Code/get_genomes.py hlorente@ucm.es Example/query_filters.txt --protein --refine ''
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

Type on terminal get_query_seqs.py -h for further information.<br />

### USAGE:<br />

**get_query_seqs.py 'e-mail' 'query.txt' **<br />

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

Protein
```
./Code/get_query_seqs.py hlorente@ucm.es Example/query_seqs.txt
```
Protei. Downloading restricting maximum number to 50
```
./Code/get_query_seqs.py hlorente@ucm.es Example/query_seqs.txt --retmax 50
```
Protei. Downloading restricting maximum number to 50 and curating
```
./Code/get_query_seqs.py hlorente@ucm.es Example/query_seqs.txt --retmax 50 --curated
```
Nucleotide. Downloading restricting maximum number to 50. Better see 2.
```
./Code/get_query_seqs.py hlorente@ucm.es Example/query_seqs.txt --retmax 50 --nucleotide
```

#### 2.- Refine search using --refine argument. **Recommended**
Default. Applies refseq[filter].<br />
Protein download, restricting maximum number to 50, curating, refine default just sequences from RefSeq.<br />
```
./Code/get_query_seqs.py hlorente@ucm.es Example/query_seqs.txt --retmax 50 --curated --refine
```
Protein. Protein. Downloading restricting maximum number to 50, curating and just sequences from RefSeq (refine argument).
```
./Code/get_query_seqs.py hlorente@ucm.es Example/query_seqs.txt --retmax 50 --curated --refine 
```
Nucleotide **recommended** use. Dowloading restricting maximum number to 50 and just transcripts (refine argument).
```
./Code/get_query_seqs.py hlorente@ucm.es Example/query_seqs.txt --retmax 50 --curated --refine 'biomol_mrna[PROP]'
```

#### 3.- Refine applied to each query
You must add a the field tags or filters after your query.<br />
As each query contains is own filters or field tags these will be applied exclusively to every query. See ./Example/query_seqs_2.txt and quer_seqs_3.txt for a scheme.<br />
Protein
```
./Code/get_query_seqs.py hlorente@ucm.es Example/query_seqs_3.txt
```
Nucleotide
```
./Code/get_query_seqs.py hlorente@ucm.es Example/query_seqs_4.txt
```

#### 4.- Refine to each query and general refine
Protein
```
./Code/get_query_seqs.py hlorente@ucm.es Example/query_seqs_3.txt --refine 'AND vertebrate'
```
Nucleotide
```
./Code/get_query_seqs.py hlorente@ucm.es Example/query_seqs_4.txt --refine 'AND vertebrate'
```
More info about filters and field tags at x and reading y.

## decompress_genomes.py<br />

### DESCRIPTION:<br />

Program for decompress downloaded genomes.<br />
It searchs for '.gz' files. Can be applied to any other '.gz' file.<br />
With data type (genomic, rna, protein) argument the program will look at ./Data/Datat_type.<br />
With directory argument  the program will iterativily within the detail path.<br />

Type on terminal decompress_genomes.py -h for further information.<br />

### USAGE:<br />

**decompress_genomes.py --data_type or --directory 'path' **<br />

### PARAMETERS:<br />

Optional parameters:<br />
**--directory** -> path to folders enclosing files<br />
**--genomic** -> looks at ./Data/Genomic<br />
**--rna** -> looks at ./Data/Rna<br />
**--protein** -> looks at ./Data/Protein<br />

## blast.py<br />

### DESCRIPTION:<br />

BLAST against genome .fna and .faa files.<br />
It requires BLAST download at https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/.

Type on terminal blast.py -h for further information.<br />

### USAGE:<br />

**blast.py --data_type or --directory 'path' **<br />

### PARAMETERS:<br />

**blast_path** -> path to bin ncbi-blast directory<br />
**query_file** -> path to your query file with your reference sequences (Fasta file with prot or nucl seqs) <br />
**query_type** -> data type of query sequences. prot or nucl <br />

Optional parameters:<br />
**--directory** -> path to folders enclosing genomes<br />
**--genomic** -> looks at ./Data/Genomic<br />
**--rna** -> looks at ./Data/Rna<br />
**--protein** -> looks at ./Data/Protein<br />
**--evalue** -> E-value threshold for search. Default value equal to 1e-10<br />
**--outfmt** -> Output alignment options. Default value 6<br />
**--out_exten** -> Extension of the BLAST output file. Default "_out.tsv"<br />


### EXAMPLES:

Examples files are availabe at Example directory.

#### 1.- Query and database (genome) are proteins

Protein directory
```
./Code/blast.py 'path_to_bin_directory' ./Example/protein_query.fas 'prot' --protein
```
Custom directory
```
./Code/blast.py 'path_to_bin_directory' ./Example/protein_query.fas 'prot' --directory 'path_to_custom_directory' --protein
```

#### 2.- Query and database (genome) are nucleotides
Nucleotide directory
```
./Code/blast.py 'path_to_bin_directory' ./Example/nucleotide_query.fas 'nucl' --rna
```
Custom directory
```
./Code/blast.py 'path_to_bin_directory' ./Example/nucleotide_query.fas 'nucl' --directory 'path_to_custom_directory' --rna
```

#### 3.- Query is protein and database (genome) is nucleotide
Nucleotide directory
```
./Code/blast.py 'path_to_bin_directory' ./Example/protein_query.fas 'prot' --rna
```
Custom directory
```
./Code/blast.py 'path_to_bin_directory' ./Example/protein_query.fas 'prot' --directory 'path_to_custom_directory' --rna
```

#### 4.- Query is nucleotide and database (genome) is protein
Protein directory
```
./Code/blast.py 'path_to_bin_directory' ./Example/protein_query.fas 'nucl' --protein
```
Custom directory
```
./Code/blast.py 'path_to_bin_directory' ./Example/protein_query.fas 'nucl' --directory 'path_to_custom_directory' --protein
```

## get_unique_hits.py<br />

### DESCRIPTION:<br />

Gets unique hits from BLAST output files based on genomes IDs.<br />
The program search for '.tsv' files. Change this with --pattern argument
Generates output with 'unique.tsv' extension.
Outputs will be stored at Genomic, Rna, Protein or custom directory.

Type on terminal get_unique_hits.py -h for further information.<br />

### USAGE:<br />

**get_unique_hits.py  --data_type or --directory 'path' **<br />

### PARAMETERS:<br />

Optional parameters:<br />
**--directory** -> path to folders enclosing genomes<br />
**--genomic** -> looks at ./Data/Genomic<br />
**--rna** -> looks at ./Data/Rna<br />
**--protein** -> looks at ./Data/Protein<br />
**--pattern** -> custom pattern to find Blast output files. Default ".tsv"<br />


### EXAMPLES:

Examples files are availabe at Example directory.

#### 1.- Default

Protein directory. Output will be stored at ./Data/Protein/Species_directory
```
./Code/get_unique_hits.py --protein
```
Custom directory. Output will be stored at custom directory (where the file is located)
```
./Code/get_unique_hits.py --directory 'path_to_custom_directory'
```

#### 2.- New pattern
Using '.txt'. Program will serach for '.txt' files.
```
./Code/get_unique_hits.py --protein --pattern '.txt'
```
