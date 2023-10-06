import random
import sys


DNA_NUCLEOTIDES = set('ATGC')
RNA_NUCLEOTIDES = set('AUGC')
PAIRS_DNA = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
PAIRS_RNA = {'a': 'u', 'u': 'a', 'c': 'g', 'g': 'c'}


def is_dna(seq: str, nucleotides: set = None) -> bool:
    """
    Checks whether a sequence is DNA or not

    Arguments:
    - seq (str): nucleotide sequence that will be checked

    Returns:
    - bool: True if the sequence is DNA, False if it's not
    """
    if nucleotides is None:
        nucleotides = DNA_NUCLEOTIDES
    values = set(seq.upper())
    return values <= nucleotides


def is_rna(seq):
    """
    Checks whether a sequence is RNA or not

    Arguments:
    - seq (str): nucleotide sequence that will be checked

    Returns:
    - bool: True if the sequence is RNA, False if it's not
    """
    return is_dna(seq, nucleotides=RNA_NUCLEOTIDES)


def transcribe(seq: str) -> str:
    """
    Imitate the transcription process. All "T" are replaced by "U"

    Arguments:
    - seq (str): DNA sequence that will be transcribed

    Returns:
    - res (str): transcribed sequence
    """
    if not is_dna(seq):
        print(f'Sequence {seq} is not DNA', file=sys.stderr)
        sys.exit(1)
    res = seq.replace('T', 'U').replace('t', 'u')
    return res


def reverse(seq: str) -> str:
    """
    Convert the sequence to reverse sequence

    Arguments:
    - seq (str): nucleotide sequence

    Returns:
    - res (str): reverse sequence
    """
    return seq[::-1]


def complement(seq: str) -> str:
    """
    Return complement sequence

    Arguments:
    - seq (str): nucleotide sequence

    Returns:
    - res (str): complement sequence
    """
    if is_dna(seq):
        res = []
        for base in seq:
            res.append(PAIRS_DNA[base] if base.islower() else PAIRS_DNA[base.lower()].upper())
    else:
        res = []
        for base in seq:
            res.append(PAIRS_RNA[base] if base.islower() else PAIRS_RNA[base.lower()].upper())
    return ''.join(res)


def reverse_complement(seq: str) -> str:
    """
    Return complement and reversed sequence

    Arguments:
    - seq (str): nucleotide sequence

    Returns:
    - res (str): complement and reversed sequence
    """
    return complement(reverse(seq))


def reverse_transcribe(seq: str) -> str:
    """
    imitate the process of reverse transcription. All "U" will be replaced with "T"

    Arguments:
    - seq (str): RNA sequence

    Returns:
    - res (str): reverse transcribed sequence
    """
    if not is_rna(seq):
        print(f'Sequence {seq} is not RNA', file=sys.stderr)
        sys.exit(1)
    res = seq.replace('U', 'T').replace('u', 't')
    return res


def mutate(seq: str) -> str:
    """
    Sequence with one random replacement.
    Any nucleotide is replaced with A, T, C, G if it is DNA, and with A, U, C, G if it is RNA.
    This preserves the register of the original sequence

    Arguments:
    - seq (str): nucleotide sequence

    Returns:
    - res (str): sequence with one replacement
    """
    mut_site = random.randint(1, len(seq)-1)
    mut_nucl = seq[mut_site]
    if is_dna(seq):
        change = random.choice(list(DNA_NUCLEOTIDES))
    else:
        change = random.choice(list(RNA_NUCLEOTIDES))
    res = seq.replace(mut_nucl, change if mut_nucl.isupper() else change.lower(), 1)
    return res


def deletion(seq: str) -> str:
    """
    Produces a random deletion if the sequence contains more than three nucleotides

    Arguments:
    - seq (str): nucleotide sequence

    Returns:
    - res (str): sequence with random deletion
    """
    del_start = random.randint(0, len(seq) - 2)
    del_end = random.randint(del_start, len(seq) - 1)
    return seq if len(seq) <= 3 else seq[:del_start] + seq[del_end:]


commands = {'transcribe': transcribe, 'reverse': reverse, 'complement': complement,
            'reverse_complement': reverse_complement, 'reverse_transcribe': reverse_transcribe,
            'mutate': mutate, 'deletion': deletion}
