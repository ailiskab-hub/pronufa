import re
from typing import List, Tuple, NoReturn


GENE_NAME_PATTERN = r'(?<=\").*(?=\")'
PATTERN_TO_NAME_GBK = r'(?<=[\\/])\w+(?=\.gbk)'
PATTERN_TO_NAME_BLAST = r'(?<=[\\/])\w+\.txt'
PATTERN_TO_SPLIT = r'\s{2,}'
BEG = '>'
END = '\n'


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> NoReturn:
    """
    Converts a multiline format FASTA file to a single-line format and saves it to a new FASTA file.

    Args:
    - input_fasta (str): The filename of the input FASTA file in multiline format.
    - output_fasta (str, optional): The filename for the output FASTA file.
        If not provided, a default filename is generated based on the input filename.
    """
    with open(input_fasta) as fa_file:
        lines = fa_file.readlines()
        seqs = {}
        curr_seq = []
        name = lines[0]
        for line in lines[1:]:
            if line.startswith('>'):
                seqs[name] = ''.join(curr_seq)
                name = line
            else:
                curr_seq.append(line.strip())

    if not output_fasta:
        output_fasta = input_fasta.rstrip('.fasta')

    with open(f'{output_fasta}.fasta', mode='w') as new_file:
        for name, sequence in seqs.items():
            new_file.write(name)
            new_file.write(sequence + '\n')


def get_gene_and_translation(lines: List[str]) -> Tuple[Dict[str, str], Dict[str, int]]:
    """
    Extracts gene names and their corresponding translations from a list of lines containing gene information.

    Args:
    - lines (List[str]): A list of lines containing gene information in a specific format.

    Returns:
    - gene_translation (Dict[str, str]): A dictionary where gene names are keys, and their translations are values.
    - gene_numeration (Dict[str, int]): A dictionary where gene names are keys, and their order of appearance is values.
    """
    index = 0
    gene_translation = {}
    gene_numeration = {}
    gene_number = 0
    while index < len(lines):
        if 'gene' in lines[index]:  # find the 'gene' in line
            gene_name = re.search(GENE_NAME_PATTERN, lines[index])[0]  # parse gene name
            index += 1
            while 'translation' not in lines[index]:  # look for beginning of the sequence
                index += 1
            if not lines[index].endswith('"\n'):  # if the sequence is more than one line
                curr_seq = [lines[index].strip().lstrip('/translation="')]
                index += 1
                while not lines[index].endswith('"\n'):
                    curr_seq.append(lines[index].strip())
                    index += 1
                curr_seq.append(lines[index].strip().rstrip('"'))
            else:  # if the sequence one line only or less
                curr_seq = [lines[index].strip().lstrip('/translation="').rstrip('"\n')]
            gene_translation[gene_name] = ''.join(curr_seq)
            gene_numeration[gene_name] = gene_number
            gene_number += 1
        index += 1
    return gene_translation, gene_numeration


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: List[str], n_before: int = 1, n_after: int = 1,
                                   output_fasta: str = None) -> NoReturn:
    """
    Extracts neighbours of specific genes and their translations from a GenBank (GBK) file and saves them in a FASTA file.

    Args:
    - input_gbk (str): The path to the GenBank file to extract genes from.
    - genes (List[str]): A list of gene names to extract.
    - n_before (int, optional): The number of genes to include before the selected gene (default is 1).
    - n_after (int, optional): The number of genes to include after the selected gene (default is 1).
    - output_fasta (str, optional): The filename for the output FASTA file.
        If not provided, a default filename is generated.
    """
    with open(input_gbk) as gbk_file:
        lines = gbk_file.readlines()
        gene_translation, gene_numeration = get_gene_and_translation(lines)
        genes_only = list(gene_translation.keys())  # for making slices
        max_num = len(genes_only) - 1
        genes_interests = []

        for gene in genes:
            gene_number = gene_numeration.get(gene)
            genes_interests_start = gene_number - n_before if gene_number > n_before else 0
            genes_interests_finish = gene_number + n_after
            if genes_interests_finish > max_num:
                genes_interests_finish = max_num

            genes_interests.extend(genes_only[genes_interests_start:gene_number])
            genes_interests.extend(genes_only[gene_number + 1:genes_interests_finish + 1])

        genes_interests = list(dict.fromkeys(genes_interests))  # collect all genes of interest and left unique only

    if not output_fasta:  # make the filename
        output_fasta = re.search(PATTERN_TO_NAME_GBK, input_gbk)[0]

    with open(f'{output_fasta}.fasta', mode='w') as new_file:
        for interesting_gene in genes_interests:
            new_file.write(BEG + interesting_gene + END)
            new_file.write(gene_translation.get(interesting_gene) + END)


def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta: str = None) -> NoReturn:
    """
    Shifts the starting position of sequences in a FASTA file and saves the modified sequences to a new FASTA file.

    Args:
    - input_fasta (str): The filename of the input FASTA file.
    - shift (int): The number of positions to shift the start of sequence.
    - output_fasta (str, optional): The filename for the output FASTA file.
        If not provided, a default filename is generated.
    """
    with open(input_fasta) as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                name = line
            else:
                sequence = line.strip()

        beginning = sequence[:shift]
        ending = sequence[shift:]

        changed_seq = ending + beginning

    if not output_fasta:
        output_fasta = input_fasta.rstrip('.fasta')

    with open(f'{output_fasta}.fasta', mode='w') as new_file:
        new_file.write(name)
        new_file.write(changed_seq + '\n')


def parse_blast_output(input_file: str, output_file: str = None) -> NoReturn:
    """
    Parses a BLAST output file to extract and sort information about some proteins.

    Args:
    - input_file (str): The filename of the input BLAST output file in TXT format.
    - output_file (str, optional): The filename for the output file. If not provided, a default filename is generated.
    """
    with open(input_file) as blast_res:
        query_num = blast_res.read().count('Query #')

    with open(input_file) as blast_res:
        num = 0
        line = blast_res.readline()
        proteins_info = []
        while num < query_num:
            num += 1
            while not line.startswith('Sequences producing significant alignments:'):
                line = blast_res.readline()

            while not line.startswith('Description'):
                line = blast_res.readline()
            line = blast_res.readline().strip()
            proteins_info.append(re.split(PATTERN_TO_SPLIT, line)[0])

        proteins_info.sort(key=lambda string: string[0].lower())

    if not output_file:
        output_file = re.search(PATTERN_TO_NAME_BLAST, input_file)[0]

    with open(output_file, mode='w') as new_file:
        for protein in proteins_info:
            new_file.write(protein + END)
