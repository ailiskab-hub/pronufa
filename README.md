# Pronufa
## Brief description
Pronufa (*Protein, nucleic, fastq*) - tool for dealing with fastq sequences and proteins or nucleic acids sequenses
Pronufa contains following functions:
- `run_fastq_tool`: tool for selection sequences
- `run_dna_rna_tools`: tool for processing DNA and RNA sequences
- `run_protein_tool`: tool for processing amini acid sequences

Directiry `my_tools` contains files `dna_rna_tools.py`,  `fastq_tool.py`,  `protein_tool.py` with some functions which are required for correct operation of functions described above.
## The structure of repository
``` python
-/
 |- pronufa.py # (impots and 3 functions)
 |- README.md
 |- my_tools/
       |- dna_rna_tools.py
       |- protein_tool.py
       |- fastq_tool.py
```

## Tool for working with fastq-sequences
The function `run_fastq_tool`is designed to select sequences that satisfy the given conditions. The function takes the following arguments:
- `seqs`: a dictionary consisting of fastq sequences with the following structure: key - string, sequence name, value - tuple of two strings: sequence and quality
- `gc_bounds`: composition GC interval (in percent) for filtering. It can take a tuple of two values, or one number - the upper limit of the interval (by default it is equal to (0, 100), i.e. all reads are saved)
- `length_bounds`: length interval for filtering. Can take a tuple of two values, or one number - the upper limit of the interval. Default is (0, 2* *32)
- `quality_threshold`: threshold value of average read quality for filtering, default is 0 (phred33 scale). Reads with average quality across all nucleotides below the threshold are discarded.

The function returns a dictionary with a structure similar to the one passed, but containing only sequences that satisfy all conditions

## Tool for working with nucleic acids sequences
The `run_dna_rna_tools` is designed to perform various operations with RNA and DNA sequences.

The `run_dna_rna_tools` function takes as input an arbitrary number of arguments with DNA or RNA sequences, as well as the name of the procedure to be performed (the last argument). After this, the command performs the specified action on all transmitted sequences

**List of procedures:**
- `transcribe` — print the transcribed sequence, accepts only DNA as input
- `reverse` — print the reversed sequence
- `complement` — print the complementary sequence
- `reverse_complement` — print the reverse complementary sequence
- `reverse_transcribe` — print the transcribed sequence, accepts only RNA as input
- `mutate` — print a sequence with one random replacement. Any nucleotide is replaced with A, T, C, G if it is DNA, and with A, U, C, G if it is RNA. This preserves the register of the original sequence
- `deletion` - prints a sequence with an arbitrary deletion if the sequence contains more than three nucleotides

## Tool for working with amino acids sequences
The function `protein_tool' is designed to perform various operations with amino acid sequences. It takes as input the name of the procedure and the sequence of amino acids, or two sequences, in the case of some procedures.


**List of procedures:**
- `calculate_amino_acid_percentages` - calculation of the relative amino acid composition
- `classify_amino_acid` - counting the relative number of amino acids by class
- `find_amino_acid_indices' - get indexes of all occurrences of AK in protein
- `counting_point_mutations` - counting point mutations and counting the percentage of similarity (sequence length/number of point mutations)
- `counting_molecular_weight' - counting molecular weight
- `get_occurrences` - find the number of occurrences in a sequence of another sequence, indexes and the number of occurrences
- `count_variant_rna` - counting RNA variants that could encode a given sequence
- `determine_total_protein_charge` - determination of the total charge of the protein
- `calculate_pi` - calculation of the approximate isoelectric point

## Usage examples
1. fastq-sequences selection
``` python
EXAMPLE_FASTQ = {'@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079804:1:SRR292678:1:1101:24563:24563': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    '@SRX079804:1:SRR292678:1:1101:30161:30161': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD')}
run_fastq_tool(EXAMPLE_FASTQ, 55, (10,100), 32)
# {'@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'), '@SRX079804:1:SRR292678:1:1101:30161:30161': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD')}
```

2. Processing nucleic acids sequences
``` python
run_dna_rna_tools('ATGccaT', 'reverse_complement') #AtggCAT
run_dna_rna_tools('ATaaTCCGattG', 'transcribe') #AUaaUCCGauuG
run_dna_rna_tools('AtCccTTG', 'mutate') #AtCtcTTG
run_dna_rna_tools('AtCGgTTTTTTTTTTcTTAAg', 'deletion') #AtCGgTTTAg
```
3. Processing amino acid sequences
``` python
run_protein_tool('ASQGAMQR', 'counting_molecular_weight') # '847'
run_protein_tool('TKKKKTDDDA', 'calculate_pI') # '7.225555555555555'
run_protein_tool('TATAQQQWRVVTDDDA', 'count_variant_rna') # '25165824'
run_protein_tool('ASQRGARWQRMQR', 'QR', 'get_occurrences') # 'Number of occurrences: 3; indexes: 3, 9, 12'
```
