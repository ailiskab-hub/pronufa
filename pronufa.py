import re
import sys
from typing import Tuple, Dict, NoReturn

from my_tools import dna_rna_tools as na
from my_tools import fastq_tool as fq
from my_tools import protein_tool as pt

PATTERN_TO_NAME = r'(?<=/)\w+(?=\.fastq)'


def run_fastq_tool(input_path: str, gc_bounds: Tuple[int, int] = (0, 100),
                   length_bounds: Tuple[int, int] = (0, 2 ** 32), quality_threshold: int = 0,
                   output_filename: str = None) -> NoReturn:
    """
    Main function that is used to select sequences that satisfy the conditions

    Arguments:
    - input_path (str): the path to the file with fastq sequences
    - gc_bounds (Tuple[int, str] or int): composition GC interval or upper limit (in percent)
    - length_bounds (Tuple[int, str] or int): length interval for filtering
    - quality_threshold (int): threshold value of average read quality for filtering
    - output_filename (str): name of output file with
    """
    names, seqs, comments, qualities = fq.read_file(input_path)

    number_of_sequence = len(names)

    selected_seqs = {}
    try:
        min_gc, max_gc = gc_bounds
    except TypeError:
        min_gc = 0
        max_gc = gc_bounds
    try:
        min_length, max_length = length_bounds
    except TypeError:
        min_length = 0
        max_length = length_bounds

    for index in range(number_of_sequence):
        name = names[index]
        seq = seqs[index]
        comment = comments[index]
        q_seq = qualities[index]
        length = len(seq)
        mean_quality = fq.count_quality(q_seq) / length
        if fq.check_length(min_length, length, max_length) and fq.count_check_gc(seq, min_gc, max_gc, length) and fq.check_quality(q_seq, quality_threshold):
                selected_seqs[name] = (seq, comment, q_seq)

    if not output_filename:
        output_filename = re.search(PATTERN_TO_NAME, input_path)[0]
    fq.write_file(selected_seqs, output_filename)


def run_dna_rna_tools(*args: str) -> str or list[str]:
    """
    Main function that is used to get sequence(s) and command.
    It performs a given action with the entered sequence

    Arguments:
    - args (str): nucleotide sequence(s) and command.
    The last element of the string must be the command

    Returns:
    - result (str or list[str]): the result of a given sequence processing
    """
    *sequences, action = args
    command = na.commands[action]
    for seq in sequences:
        if not (na.is_dna(seq) or na.is_rna(seq)):
            print(f'Sequence {seq} is not RNA or DNA', file=sys.stderr)
            sys.exit(1)
    return [command(seq) for seq in sequences] if len(sequences) > 1 else command(sequences[0])


def run_protein_tool(*args: str) -> str:
    """
    Main function that is used to get sequence(s) and command. It performs a given action with the entered sequence

    Arguments:
    - args (str): amino acid sequence(s) and command.
    The input must use the single letter amino acid code
    The last element of the string must be the command

    Returns:
    - result (str): the result of a given sequence processing
    """
    *sequences, action = args
    sequences = [seq.upper() for seq in sequences]
    command = pt.commands[action]
    for seq in sequences:
        if not pt.is_protein(seq):
            print('Sequence is not protein', file=sys.stderr)
            sys.exit(1)
    return str(command(sequences[0], sequences[1])) if len(sequences) == 2 else str(command(sequences[0]))
