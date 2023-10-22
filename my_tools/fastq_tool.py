from collections import Counter
from typing import Tuple, Dict, List, NoReturn
import os


END = '\n'


def check_length(min_length: int, length: int, max_length: int) -> bool:
    """
    Check whether the length bigger than minimal length but less than maximum length

    Arguments:
    - min_length (int): minimal length
    - length (int): length of the sequence that we want to check
    - max_length (int): maximum length

    Return:
    - bool: True if the sequence is matches the condition, False if it's not
    """
    return min_length <= length <= max_length


def count_check_gc(seq: str, min_gc: int, max_gc: int, length: int) -> bool:
    """
    Check whether the sequence is bigger than minimal gc bound but less than maximum gc bound

    Arguments:
    - seq (str): sequence
    - min_gc (int): minimal gc bound
    - max_gc (int): maximum gc bound
    - length (int): length of the sequence that we want to check

    Return:
    - bool: True if the sequence is matches the condition, False if it's not
    """
    cnt = Counter(seq)
    gc = 100 * (cnt['G'] + cnt['C']) / length
    return min_gc <= gc <= max_gc


def count_quality(q_seq: str) -> int:
    """
    Count quality of the sequence

    Arguments:
    - q_seq (str): the quality of each nucleotide in the sequence

    Return:
    - int: total value of sequence quality
    """
    res = 0
    for score in q_seq:
        res += ord(score) - 33
    return res


def read_file(input_path: str) -> Tuple[List[str], List[str], List[str], List[str]]:
    """
    Read file with sequences

    Arguments:
    - input_path (str): the path to the file with fastq sequences

    Return:
    - Tuple[List[str], List[str], List[str], List[str]]: four lists with names, sequences, comments and quality og the sequences
    """
    names = []
    seqs = []
    comments = []
    qualities = []
    with open(input_path, 'r') as fastq_file:
        lines = fastq_file.readlines()
        for line_ind in range(0, len(lines), 4):
            names.append(lines[line_ind].strip())
            seqs.append(lines[line_ind + 1].strip())
            comments.append(lines[line_ind + 2].strip())
            qualities.append(lines[line_ind + 3].strip())
    return names, seqs, comments, qualities


def write_file(selected_seqs: Dict[str, Tuple[str, str, str]], output_filename: str) -> NoReturn:
    """
    Function that is used to write selected sequences to the file

    Arguments:
    - seqs (Dict[str, Tuple[str, str]]): a dictionary consisting of selected sequences
    - output_filename (str): name for the file
    """

    if not os.path.isdir("fastq_filtrator_resuls"):
        os.mkdir("fastq_filtrator_resuls")

    with open(f'fastq_filtrator_resuls/{output_filename}.fastq', mode='w') as new_file:
        for name, content in selected_seqs.items():
            new_file.write(name+END)
            new_file.write(content[0] + END)
            new_file.write(content[1] + END)
            new_file.write(content[2] + END)
