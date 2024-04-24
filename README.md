# Pronufa
Pronufa is the module for processing and working with biological sequences and data. 

There are files `pronufa.py`, `bio_files_processor.py`, `Showcases.ipynb`, `custom_random_forest.py` and `test_pronufa.py`

## Table of contents
  * [The structure of repository](#structure)
  * [Brief description of files in repo](#description_1)
    

## The structure of repository <a name="structure"></a>
``` python
-/
 |- pronufa.py 
 |- bio_files_processor.py
 |- Showcases.ipynb
 |- custom_random_forest.py
 |- test_pronufa.py
 |- requirements.txt
 |- README.md
 |- data/
       |- example_fasta.fasta
       |- example_gbk.gbk
       |- sequence (2).fasta
       |- test_file.fastq
```

## Brief description of files in repo <a name="description_1"></a>
**Bio_files_processor** - tool for working with files
Bio_files_processor contains functions: `convert_multiline_fasta_to_oneline`, `get_gene_and_translation`, `select_genes_from_gbk_to_fasta`, `change_fasta_start_pos`, `parse_blast_output`, dataclass `FastaRecord` and class `OpenFasta`

**Pronufa** - toolkit provides a set of classes and functions for working with biological sequences, including DNA, RNA, and amino acid sequences. It offers functionalities such as sequence type checking, complement generation, GC content calculation, sequence transcription, reverse transcription, and isoelectric point calculation for amino acid sequences. 
Class `GenscanOutput` and function `run_genscan` allows to get GenScan prediction for the given sequence or sequence file. 

**Custom_random_forest** contains a custom Random Forest implementation with the ability to train a random forest in parallel and use parallelism for predictions.

**Test_pronufa** file with some tests for `pronufa` and `bio_files_processor`

**Showcases** notebook which provides examples of the module's functionality


