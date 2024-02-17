# Pronufa
There are two files `pronufa.py` and `bio_files_processor.py`

## Table of contents
  * [Brief description of the bio_files_processor.py](#description_1)
  * [Brief description of the pronufa.py](#description_2)
  * [The structure of repository](#structure)
  * [Tool for working with fastq-sequences](#fastq)
  * [Classes for working with biological acids sequences](#na)
  * [Usage examples](#example)

## Brief description of the bio_files_processor.py <a name="description_1"></a>
Bio_files_processor - tool for working with files
Bio_files_processor contains following functions:
- `convert_multiline_fasta_to_oneline`: Converts a multiline format FASTA file to a single-line format and saves it to a new FASTA file
- `get_gene_and_translation`: Extracts gene names and their corresponding translations from a list of lines containing gene information
- `select_genes_from_gbk_to_fasta`: Extracts neighbours to specific genes and their translations from a GenBank (GBK) file and saves them in a FASTA file
- `change_fasta_start_pos`: Shifts the starting position of sequences in a FASTA file and saves the modified sequences to a new FASTA file
- `parse_blast_output`: Parses a BLAST output file to extract and sort information about proteins

## Brief description of the pronufa.py <a name="description_2"></a>
This toolkit provides a set of classes and functions for working with biological sequences, including DNA, RNA, and amino acid sequences. It offers functionalities such as sequence type checking, complement generation, GC content calculation, sequence transcription, reverse transcription, and isoelectric point calculation for amino acid sequences.
**✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ:**


## The structure of repository <a name="structure"></a>
Directiry `my_tools` contains files `dna_rna_tools.py`,  `fastq_tool.py`,  `protein_tool.py` with some functions which are required for correct operation of functions described above.
``` python
-/
 |- pronufa.py # (impots and 3 functions)
 |- bio_files_processor.py
 |- README.md
 |- my_tools/
       |- dna_rna_tools.py
       |- protein_tool.py
       |- fastq_tool.py
```
**✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ:**

## Tool for working with fastq-sequences <a name="fastq"></a>
The function `filter_fastq`is designed to select sequences that satisfy the given conditions. The function takes the following arguments:
-  `input_path`: the path to the file with fastq sequences
- `gc_bounds`: composition GC interval (in percent) for filtering. It can take a tuple of two values, or one number - the upper limit of the interval (by default it is equal to (0, 100), i.e. all reads are saved)
- `length_bounds`: length interval for filtering. Can take a tuple of two values, or one number - the upper limit of the interval. Default is (0, 2* *32)
- `quality_threshold`: threshold value of average read quality for filtering, default is 0 (phred33 scale). Reads with average quality across all nucleotides below the threshold are discarded.
- `output_filename`: name of output file

The function write sequences that satisfy all conditions to the file

**✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ: ✧･ﾟ:**

## Classes for working with biological acids sequences <a name="na"></a>
`DNASequence, RNASequence, AminoAcidSequence` are designed to represent biological sequences
  
## Usage examples <a name="example"></a>
1. fastq-sequences selection
``` python
EXAMPLE_FASTQ = {'@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079804:1:SRR292678:1:1101:24563:24563': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    '@SRX079804:1:SRR292678:1:1101:30161:30161': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD')}
filter_fastq(EXAMPLE_FASTQ, 55, (10,100), 32)
# {'@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'), '@SRX079804:1:SRR292678:1:1101:30161:30161': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD')}
```

2. Processing nucleic acids sequences
``` python
# Create an instance of DNASequence
dna_seq = DNASequence("ATCG")

# Transcribe DNA sequence into RNA
rna_seq = dna_seq.transcribe()

# Calculate GC content of DNA sequence
gc_content = dna_seq.gc_content()

print(f"Transcribed RNA sequence: {rna_seq}")      # UAGC
print(f"GC content of DNA sequence: {gc_content}") # 0.5
```
3. Processing amino acid sequences
``` python
prot = AminoAcidSequence('LPTQQ')
prot.calculate_pi() # 5.885
```
