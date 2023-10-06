from collections import Counter


NUMBER_CODONS = {'F': 2, 'L': 6, 'I': 3, 'M': 1, 'V': 4, 'S': 6, 'P': 4, 'T': 4,
                 'A': 4, 'Y': 2, 'H': 2, 'Q': 2, 'N': 2, 'K': 2, 'D': 2, 'E': 2,
                 'C': 2, 'W': 1, 'R': 6, 'G': 4
                 }
NEG_CHARGED = ['D', 'E']
POS_CHARGED = ['H', 'K', 'R']
PK1 = {'F': 2.2, 'L': 2.36, 'I': 2.36, 'M': 2.28,
       'V': 2.32, 'S': 2.21, 'P': 1.99, 'T': 2.71,
       'A': 2.34, 'Y': 2.2, 'H': 1.82, 'Q': 2.17,
       'N': 2.02, 'K': 2.18, 'D': 1.88, 'E': 2.19,
       'C': 1.71, 'W': 2.38, 'R': 2.17, 'G': 2.34
       }
PK2 = {'F': 9.09, 'L': 9.6, 'I': 9.68, 'M': 9.21,
       'V': 9.62, 'S': 9.15, 'P': 10.96, 'T': 9.62,
       'A': 9.69, 'Y': 9.11, 'H': 9.17, 'Q': 9.13,
       'N': 9.8, 'K': 8.95, 'D': 9.6, 'E': 9.67,
       'C': 8.33, 'W': 9.39, 'R': 9.04, 'G': 9.6
       }
PK3 = {'Y': 10.07, 'H': 6.0, 'K': 10.53,
       'C': 10.78, 'D': 3.65, 'E': 4.25, 'R': 12.48
       }
AA_SET = set('FLIMVSPTAYHQNKDECWRG')
DICT_MOLECULAR_MASS = {
    'G': 75, 'A': 89, 'V': 117, 'L': 131, 'I': 131, 'P': 115,
    'F': 165, 'Y': 181, 'W': 204, 'S': 105, 'T': 119, 'C': 121,
    'M': 149, 'N': 132, 'Q': 146, 'D': 133, 'E': 147, 'K': 146,
    'R': 174, 'H': 155
}


def calculate_amino_acid_percentages(seq: str) -> str:
    pass
#     """
#     Calculating the percentage of amino acids in protein.
#
#     Arguments:
#     - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code.
#
#     Returns:
#     - output (str): percentage of amino acids in protein in descending order.
#     """
#     aa_count = {}  # counting amino acids in sequence
#     for amino_acid in seq:
#         if amino_acid in aa_count:
#             aa_count[amino_acid] += 1
#         else:
#             aa_count[amino_acid] = 1
#     composition_rates = {}
#     for aa in aa_count:
#         composition_rates[aa] = aa_count[aa] / len(seq) * 100
#     output = ', '.join([f'{key}: {round(value, 2)}' for key, value in sorted(composition_rates.items(),
#                                                                              key=lambda item: -item[1])])
#     return output


def classify_amino_acid(seq: str) -> str:
    pass
#     """
#     Determine the percentage of acidic, basic and neutral amino acids in protein.
#
#     Arguments:
#     - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code.
#
#     Returns
#     - output (str): percentage of neutral, acidic and basic amino acids in protein.
#     """
#     amino_acid_counts = {'acidic': 0, 'neutral': 0, 'basic': 0}
#     for amino_acid in seq:
#         if amino_acid in ['D', 'E']:
#             amino_acid_counts['acidic'] += 1
#         elif amino_acid in ['A', 'N', 'C', 'Q', 'G', 'I', 'L', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']:
#             amino_acid_counts['neutral'] += 1
#         elif amino_acid in ['H', 'R', 'K']:
#             amino_acid_counts['basic'] += 1
#     acidic_percentage = round(amino_acid_counts['acidic'] / len(seq) * 100, 2)
#     neutral_percentage = round(amino_acid_counts['neutral'] / len(seq) * 100, 2)
#     basic_percentage = round(amino_acid_counts['basic'] / len(seq) * 100, 2)
#     output = f'neutral: {neutral_percentage}, acidic: {acidic_percentage}, basic: {basic_percentage}'
#     return output


def counting_point_mutations(seq1: str, seq2: str) -> int:
    """
    Counts the number of mutations - amino acid substitutions in the sequence seq2 relative to seq1.
    Input sequences must have the same length.

    Arguments:
    - seq1 (str): sequence to compare with
    - seq2 (str): sequence to compare to

    Return:
    - output (int): number of amino acid substitutions
    """
    output = 0
    for number_amino_acid in range(len(seq1)):
        if seq1[number_amino_acid] != seq2[number_amino_acid]:
            output += 1
    return output


def counting_molecular_weight(seq: str) -> int:
    """
    Counts the molecular mass of a protein sequence seq

    Arguments:
    - seq (str): sequence to count the molecular weight

    Return:
    - output (int): molecular weight value
    """
    output = 0
    for amino_acid in seq:
        output += DICT_MOLECULAR_MASS[amino_acid]
    return output - 18 * (len(seq) - 1)


def get_occurrences(seq1: str, seq2: str) -> str:
    """
    Counting the number of occurrences of string seq2 in string seq1.
    Getting indexes of occurrences of string seq2 in string seq1.

    Arguments:
    - seq1 (str): sequence in which search
    - seq2 (str): sequence to search in

    Return:
    - output (str): str, first element is the number of occurrences (int).
      All subsequent elements are indexes of occurrences of seq2 in seq1.
    """
    output = [seq1.count(seq2)]
    for i in range(len(seq1)):
        if seq1.startswith(seq2, i):
            output.append(i + 1)
    return (f'Number of occurrences: {output[0]}; '
            f'indexes: {", ".join(str(element) for element in output[1:])}')


def find_amino_acid_indices(seq: str, amino_acid: str) -> str:
    pass
#     """
#     Finds the amino acid indices specified in the input.
#
#     Arguments:
#     - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code.
#     - amino_acid (str): amino acid for which indices need to be found.
#
#     Returns:
#     - output (str): all found indices in the protein for which the entered amino acid corresponds to.
#     """
#     indices = []
#     if amino_acid not in seq:
#         raise ValueError('Amino acid not found')
#     for index, aa in enumerate(seq):
#         if aa == amino_acid:
#             indices.append(index + 1)
#     output = ', '.join(str(i) for i in indices)
#     return output


def count_variant_rna(seq: str) -> int:
    """
    Сounts the number of RNAs that can be a template for the synthesis of the entered sequence

    Arguments:
    - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code

    Returns:
    - output (int): number of RNAs that can be a template for the synthesis of the entered sequence
    """
    output = 1
    for i in seq:
        output = output * NUMBER_CODONS[i]
    return output


def determine_total_protein_charge(seq) -> str:
    """
    Determine whether the protein has positive, negative or neutral charge in neutral pH

    Arguments:
    - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code

    Returns:
    - output (str): positive, negative or neutral charge of protein in neutral pH
    """
    seq_list = list(seq.strip())
    aa_cnt = Counter(seq_list)
    number_of_pos = sum([aa_cnt[aa] for aa in POS_CHARGED])
    number_of_neg = sum([aa_cnt[aa] for aa in NEG_CHARGED])
    if number_of_pos == number_of_neg:
        return 'neutral'
    return 'positive' if number_of_pos > number_of_neg else 'negative'


def calculate_pi(seq: str) -> float:
    """
    Calculation pI of the protein in neutral pH

    Arguments:
    - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code

    Returns:
    - output (float): approximate value of the pI of the protein in neutral pH
    """
    seq_list = list(seq.strip())
    first_aa = seq_list[0]
    last_aa = seq_list[-1]
    aa_cnt = Counter(seq_list)

    summ_charge = [PK2[first_aa], PK1[last_aa]]

    for key, value in aa_cnt.items():
        try:
            summ_charge.extend([PK3[key] for _ in range(value)])
        except KeyError:
            pass

    return sum(summ_charge) / len(summ_charge)


def is_protein(seq: str) -> bool:
    """
    Check whether the transmitted sequence consists of amino acids

    Arguments:
    - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code
    - aa_set(set): all amino acid that we use. May be replaced with extra amino acids or it's modifications.
    But other functions are not intended to work with unusual amino acids

    Returns:
    - bool: True if the sequence is an amino acid sequence, False if it contains other symbols
    """
    return True if set(seq) <= AA_SET else False


commands = {'calculate_amino_acid_percentages': calculate_amino_acid_percentages,
            'classify_amino_acid': classify_amino_acid,
            'find_amino_acid_indices': find_amino_acid_indices,
            'counting_point_mutations': counting_point_mutations,
            'counting_molecular_weight': counting_molecular_weight,
            'get_occurrences': get_occurrences,
            'count_variant_rna': count_variant_rna,
            'determine_total_protein_charge': determine_total_protein_charge,
            'calculate_pi': calculate_pi}
# def protein_tool(*args: str) -> str:
#     """
#     Main function that is used to get sequence(s) and command. It performs a given action with the entered sequence
#
#     Arguments:
#     - args (str): amino acid sequence(s) and command.
#     The input must use the single letter amino acid code
#     The last element of the string must be the command
#
#     Returns:
#     - result (str): the result of a given sequence processing
#     """
#
#     *sequences, action = args
#     sequences = [seq.upper() for seq in sequences]
#     commands = {'calculate_amino_acid_percentages': calculate_amino_acid_percentages,
#                 'classify_amino_acid': classify_amino_acid,
#                 'find_amino_acid_indices': find_amino_acid_indices,
#                 'counting_point_mutations': counting_point_mutations,
#                 'counting_molecular_weight': counting_molecular_weight,
#                 'get_occurrences': get_occurrences,
#                 'count_variant_rna': count_variant_rna,
#                 'determine_total_protein_charge': determine_total_protein_charge,
#                 'calculate_pI': calculate_pI}
#     command = commands[action]
#     for seq in sequences:
#         if not is_protein(seq):
#             print('Sequence is not protein', file=sys.stderr)
#             sys.exit(1)
#
#     return str(command(sequences[0], sequences[1])) if len(sequences) == 2 else str(command(sequences[0]))

# Проверка
# c = ['calculate_amino_acid_percentages', 'classify_amino_acid', 'find_amino_acid_indices', 'counting_point_mutations',
#      'counting_molecular_weight', 'get_occurrences', 'count_variant_rna', 'determine_total_protein_charge',
#      'calculate_pI']
#
# for com in c:
#     try:
#         print(com, protein_tool('ASQRGARWQRMQR', 'ASQRGARWARMQR', com))
#     except:
#         try:
#             print(com, protein_tool('ASQRGARWQRMQR', 'A', com))
#
#         except:
#             print(com, protein_tool('ASQRGARWQRMQR', com))
