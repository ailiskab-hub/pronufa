from collections import Counter


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
