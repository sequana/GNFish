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
**--directory** -> sets path to custom folder<br />
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
**--directory** -> sets path to custom folder<br />
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
Rna directory
```
./Code/blast.py 'path_to_bin_directory' ./Example/nucleotide_query.fas 'nucl' --rna
```
Custom directory
```
./Code/blast.py 'path_to_bin_directory' ./Example/nucleotide_query.fas 'nucl' --directory 'path_to_custom_directory' --rna
```

#### 3.- Query is protein and database (genome) is nucleotide
Rna directory
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
The program search for '.tsv' files. Change this with --pattern argument.<br />
Generates output with 'unique.tsv' extension.<br />
Outputs will be stored at Genomic, Rna, Protein or custom directory.<br />

Type on terminal get_unique_hits.py -h for further information.<br />

### USAGE:<br />

**get_unique_hits.py  --data_type or --directory 'path' **<br />

### PARAMETERS:<br />

Optional parameters:<br />
**--directory** -> sets path to custom folder<br />
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

## get_RAW_sequences.py <br />

### DESCRIPTION:<br />

Gets RAW sequences from genomes based on 'unique.tsv' files.<br />
The program search for 'unique.tsv' and '.f[a,n]a files. Change this with --blast_pattern an--genome_pattern arguments.<br />
Genome and Blast unique output file must have a "Genus_species_assemblyID" structure or the name provided by NCBI when downloaded.<br />

Type on terminal get_unique_hits.py -h for further information.<br />

### USAGE:<br />

**get_RAW_sequences.py  --data_type or --directory 'path' **<br />

### PARAMETERS:<br />

Optional parameters:<br />
**--directory** -> sets path to custom folder<br />
**--genomic** -> looks at ./Data/Genomic<br />
**--rna** -> looks at ./Data/Rna<br />
**--protein** -> looks at ./Data/Protein<br />
**--blast_pattern** -> custom pattern to find BLAST output files. Default ".tsv"<br />
**--genome_pattern** -> pattern to find genome files. Default ".f[a,n]a"<br />
**--in_len** -> Number of sites extracted upstream and downstream from the blast hit. Default 10000. A whole sequence of at least 20000 sites if exists<br />
**--query_seqs** -> Use it when you want to attach some query sequences according to BLAST results for future alignments<br />
**--query_seqs_num** -> Maximum number of query sequences extracted. Use it when you want to attach some sequences to genomic sequences for future alignments. 5 by default<br />
**--email** -> Mandatory when you want to download some sequences to complete the RAW files for future alignments. As in get_query_seqs<br />


### EXAMPLES:

Examples files are availabe at Example directory.

#### 1.- Default
Protein directory. Output will be stored at ./Data/Protein/Species_directory
```
./Code/get_RAW_sequences.py --protein
```
Genomic directory, sequences of 20,000 positions if possible. Output will be stored at ./Data/Genomic/Species_directory
```
./Code/get_RAW_sequences.py --genomic
```
Custom directory. Output will be stored at custom directory (where the file is located)
```
./Code/get_RAW_sequences.py --directory 'path_to_custom_directory'
```

#### 2.- Adding seqs for future alingment. Recommended for Genomic sequences
Genomic directory, sequences of 20,000 positions if possible. Output will be stored at ./Data/Genomic/Species_directory<br />
20 sequences from the dataset of query_seqs.txt.
```
./Code/get_RAW_sequences.py --genomic --in_len --query_seqs ./Example/query_seqs.txt --query_seqs_num 20 --email hlorente@ucm.es
```

#### 3.- Changing genome and BLAST pattern
Genomic directory, sequences of 20,000 positions if possible. Output will be stored at ./Data/Genomic/Species_directory<br />
20 sequences from the dataset of query_seqs.txt.
Program will look for 'unique.txt' and for '.fas'.
```
./Code/get_RAW_sequences.py --genomic --in_len --blast_pattern 'unique.txt' --genome_pattern '.fas' --query_seqs ./Example/query_seqs.txt --query_seqs_num 20 --email hlorente@ucm.es
```
#### 4.- Changing length of genomic extracted sequence
Genomic directory, sequences of 30,000 positions if possible. Output will be stored at ./Data/Genomic/Species_directory
```
./Code/get_RAW_sequences.py --genomic --in_len 15000
```
## align_sequences.py <br />

### DESCRIPTION:<br />

Aligns files enclosed at Data folder or your custom file.<br />
You must have used --query_seqs argument in get_RAW_sequences.py or have manually added some sequences to raw files.<br />
See https://mafft.cbrc.jp/ for information about algorithms.<br />
Type on terminal get_unique_hits.py -h for further information.<br />

### USAGE:<br />

**align_sequences.py  --data_type or --directory 'path' **<br />

### PARAMETERS:<br />

Optional parameters:<br />
**--directory** -> sets path to custom folder<br />
**--genomic** -> looks at ./Data/Genomic<br />
**--rna** -> looks at ./Data/Rna<br />
**--protein** -> looks at ./Data/Protein<br />
**--pattern** -> sets custom pattern to find files for alignment. Default 'RAW.fas'. Backup search RAW.fas<br />
**--algorithm** -> sets MAFFT algorithm. Default mafft (automatic)<br />


### EXAMPLES:

Examples files are availabe at Example directory.

#### 1.- Automatic
Genomic directory. Output will be stored at ./Data/Genomic/Species_directory
```
./Code/align_sequences.py --genomic
```
Custom directory. Output will be stored at custom directory (where the file is located)
```
./Code/align_sequences.py --directory 'path_to_custom_directory'
```

#### 2.- Custom directory, e-insi algorithma and custom pattern
Custom directory. Output will be stored at custom directory (where the file is located)
```
./Code/align_sequences.py --directory 'path_to_custom_directory' --algorithm einsi --pattern 'RAW_mod.fas'
```

## translate_sequences.py  <br />

### DESCRIPTION:<br />

Translates files.<br />
Open code file of the program for checking for different genetic codes.<br />

### USAGE:<br />

**translate_sequences.py  --data_type or --directory 'path' **<br />

### PARAMETERS:<br />

Optional parameters:<br />
**--directory** -> sets path to custom folder<br />
**--genomic** -> looks at ./Data/Genomic<br />
**--rna** -> looks at ./Data/Rna<br />
**--pattern** -> custom pattern to find files for translation. Defaul '_final.fas' and 'RAW.fas' as backup<br />
**--genetic_code** -> Default 1 -Standard. Look at code file for more genetics codes<br />
**--out_exten** -> extension of the output file. Default '_transtaled.fas'.<br />

### EXAMPLES:

#### 1.- Default
Rna directory. Output will be stored at ./Data/Rna/Species_directory
```
./Code/translate_sequences.py --rna
```
Custom directory. Output will be stored at custom directory (where the file is located)
```
./Code/translate_sequences.py --directory 'path_to_custom_directory'
```

#### 2.- Changing genetic code
Custom directory. Output will be stored at custom directory (where the file is located)
```
./Code/get_RAW_sequences.py --directory 'path_to_custom_directory' --genetic_code 4
```

## get_combined_seqs.py  <br />

### DESCRIPTION:<br />

Combines FASTA files in just one matrix.<br />
You can combine genomic, RNA and protein data. But be careful with duplications. Besides genomic and RNA should have been translated (see above) in the correct reading framework before translation (see above).<br />
Note that the program concatenate sequences to same file if you do not change file_name argument.<br />
After the pipeline (see **article**) you may have:<br />
1.-Untouched sequences. Typically protein or RNA sequence (valid for future nucleotide matrix if you do not care about UTRs; not valid for future protein matrix, because of lost of reading framework). Extension 'RAW.fas'.<br />
2.- Aligned and trimmed (mismatched positions) sequences. Mandatory for genomic sequences and in some cases for RNA sequences (if you are going to uses for protein matrix after translation, but in this case they will follow point 3). Extension '_final.fas'. Add the extension to your file after visualizing the alingment even when you do not trimmed any position of the alingment.<br />
3.- Translated sequences. Mandatory for genomic sequences and RNA sequences<br />
4.- Raw translated sequences. If the RNA dataset you are working with is just the coding sequence and you do not need to trimm it. Extension '_RAW_translated.fas<br />

### USAGE:<br />

**get_combined_seqs.py  --data_type or --directory 'path' **<br />

### PARAMETERS:<br />
**output_name** -> name for the file of the combined matrix<br />
**output_path** -> path were to store the combined matrix<br />

Optional parameters:<br />
**--directory** -> sets path to custom folder<br />
**--genomic** -> looks at ./Data/Genomic<br />
**--rna** -> looks at ./Data/Rna<br />
**--protein** -> looks at ./Data/Protein<br />
**--pattern** -> custom pattern to find files for translation. Defaul '_translated.fas' and 'RAW.fas' as backup'<br />
**--out_exten** -> extension of the output file. Default "_all_combined.fas"<br />

### EXAMPLES:

#### 1.- Default
Protein and Genomic directory. Output will be stored at ./Data/Rna/Species_directory
```
./Code/get_combined_seqs.py --protein --genomic
```
Custom directory. Output will be stored at custom directory (where the file is located)
```
./Code/get_combined_seqs.py --directory 'path_to_custom_directory'
```

#### 2.- Changing pattern
Protein and Genomic directory. Output will be stored at ./Data/Rna/Species_directory
```
./Code/get_combined_seqs.py --protein --genomic --pattern 'modified.fas'
```

## alignment_trimming.py  <br />

### DESCRIPTION:<br />

Trims alignment.<br />
See http://trimal.cgenomics.org/use_of_the_command_line_trimal_v1.2 for more information about trimAL and running algorithms.<br />

### USAGE:<br />

**alignment_trimming.py trimal_path 'path'  input_file 'path' output_file 'path' **<br />

### PARAMETERS:<br />
**trimal_path** -> path to trimAl directory<br />
**input_file** -> alignment (FASTA file with prot or nucl aligned seqs)<br />
**output_file** -> output file name<br />

Optional parameters:<br />
**--trimal_command** -> command to run trimal options. Default -gt 0.1<br />

### EXAMPLES:

Examples files are availabe at Example directory.

#### 1.- Default
Remove all positions in the alignment with gaps in 90% 
```
./Code/alignment_trimming.py trimal_path 'path_to_trimAl_source_directory' --input_file 'path_to_alignment' --output_file 'path_to_ouput_file'
```

#### 2.- Changing algorithm
Remove all positions in the alignment with gaps in 50% 
```
./Code/alignment_trimming.py trimal_path 'path_to_source_directory' --input_file 'path_to_alignment' --output_file 'path_to_ouput_file' --trimal_command -gt 0.5
```

## phylogenetic_inference.py  <br />

### DESCRIPTION:<br />

Phylogenetic analysis using IQ-TREE.<br />
See http://www.iqtree.org/doc/ for more information about IQ-TREE and running algorithms.<br />

### USAGE:<br />

**phylogenetic_inference.py  --data_type or --directory 'path' **<br />

### PARAMETERS:<br />
**iqtree_path** -> path to IQ-TREE bin folder<br />
**input_file** -> alignment (FASTA file with prot or nucl aligned seqs)<br />

Optional parameters:<br />
**--trimal_command** -> command to run IQ-Tree, without output declaration, see next parameter. Default TEST, ULTRAFAST BOOSTRAP 1000, ALRT 1000<br />
**--output_suffix** -> suffixes added to input file. Recommended, model and type of support. Default _TEST_UFBS_alrt<br />

### EXAMPLES:

#### 1.- Default
Phylogenetic inference running TEST for calculating model of evolutoon, 1000 replicates of Ultrafast Boostrap and 1000 replicates of alrt as support.
```
./Code/alignment_trimming.py iqtree_path 'path_to_IQTREE_bin_directory' --input_file 'path_to_alignment'
```

#### 2.- Changing algorithm
Phylogenetic inference running LG as model of evolution, 1000 replicates of Boostrap and 1000 replicates of alrt as support.
```
./Code/alignment_trimming.py iqtree_path 'path_to_IQTREE_bin_directory' --input_file 'path_to_alignment' --trimal_command '-m LG -b 1000 -alrt 1000
```
