import os
import re
import requests
from Bio import SeqIO, SeqUtils
from typing import Tuple, NoReturn
from collections import Counter
from dataclasses import dataclass
from bs4 import BeautifulSoup

AA_SET = set('FLIMVSPTAYHQNKDECWRG')
DNA_NUCLEOTIDES = set('ATGC')
RNA_NUCLEOTIDES = set('AUGC')
PAIRS_DNA = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
PAIRS_RNA = {'a': 'u', 'u': 'a', 'c': 'g', 'g': 'c'}
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


def filter_fastq(input_path: str, gc_bounds: Tuple[int, int] = (0, 100),
                 length_bounds: Tuple[int, int] = (0, 2 ** 32), quality_threshold: int = 0,
                 output_filename: str = None) -> NoReturn:
    """
    Filters a FASTQ file based on specified criteria and writes the filtered sequences to a new FASTQ file.

    Args:
        input_path (str): Path to the input FASTQ file.
        gc_bounds (Tuple[int, int], optional): Tuple specifying the lower and upper bounds for GC content percentage.
                Defaults to (0, 100).
        length_bounds (Tuple[int, int], optional): Tuple specifying the lower and upper bounds for sequence length.
                Defaults to (0, 2 ** 32).
        quality_threshold (int, optional): Minimum quality threshold for the sequence. Defaults to 0.
        output_filename (str, optional): Name of the output FASTQ file. If None, default name will be used.
                Defaults to None.

    Returns:
        NoReturn: This function does not return any value.
    """
    if not os.path.isdir("fastq_filtrator_resuls"):
        os.mkdir("fastq_filtrator_resuls")
    with open(f'fastq_filtrator_resuls/{output_filename}.fastq', mode='w'):
        pass

    try:
        min_length, max_length = length_bounds
    except TypeError:
        min_length = 0
        max_length = length_bounds

    try:
        min_gc, max_gc = gc_bounds
    except TypeError:
        min_gc = 0
        max_gc = gc_bounds

    for record in SeqIO.parse(open(input_path), "fastq"):
        name, seq, description, quality = record.id, record.seq, record.description, record.letter_annotations[
            "phred_quality"]
        length = len(seq)
        q_seq = sum(quality) / length
        gc = SeqUtils.gc_fraction(seq) * 100
        if min_length <= length <= max_length and \
                min_gc <= gc <= max_gc and \
                q_seq >= quality_threshold:
            with open(f'fastq_filtrator_resuls/{output_filename}.fastq', mode='a') as new_file:
                new_file.write(f'{record.format("fastq")} \n')


class BiologicalSequence(str):
    """
    Represents a biological sequence, such as DNA, RNA, or amino acid sequence.

    Attributes:
        seq (str): The biological sequence.
        mol_type (str): The type of molecule in the sequence (e.g., DNA, RNA, AA_seq).
    """

    def __init__(self, seq):
        self.seq = seq
        self.mol_type = None

    def check_type(self, assumed_type):
        """
        Checks if the sequence conforms to a given type of molecule.

        Args:
            assumed_type (str): The assumed type of molecule (e.g., DNA, RNA, AA_seq).

        Returns:
            bool: True if the sequence matches the assumed type, False otherwise.

        Raises:
            KeyError: If the assumed type is not recognized.
        """

        mol_types = {'DNA': DNA_NUCLEOTIDES,
                     'RNA': RNA_NUCLEOTIDES,
                     'AA_seq': AA_SET}

        if assumed_type not in mol_types:
            raise KeyError('Unknown type suggested')
        values = set(self.seq.upper())
        return values <= mol_types[assumed_type]

    def __repr__(self):
        return f'Sequence: {self.seq}'

    def __str__(self):
        return self.seq


class TypeDefinedSequence(BiologicalSequence):
    """
    Represents a biological sequence with a defined type (DNA, RNA, or amino acid sequence).

    Inherits from BiologicalSequence class.

    Attributes:
        seq (str): The biological sequence.
        mol_type (str): The type of molecule in the sequence (DNA, RNA, AA_seq).
    """

    def __init__(self, seq):
        super().__init__(seq)
        if self.check_type('DNA'):
            self.mol_type = 'DNA'
        elif self.check_type('RNA'):
            self.mol_type = 'RNA'
        elif self.check_type('AA_seq'):
            self.mol_type = 'AA_seq'
        else:
            raise TypeError(f'Sequence {self.seq} cannot be interpreted')


class NucleicAcidSequence(TypeDefinedSequence):
    """
    Represents a nucleic acid sequence, which can be either DNA or RNA.

    Inherits from TypeDefinedSequence class.

    Attributes:
        seq (str): The nucleic acid sequence.
        mol_type (str): The type of molecule in the sequence (DNA or RNA).
        complement_pairs (dict): Dictionary containing complement pairs for the sequence.
    """

    def __init__(self, seq):
        super().__init__(seq)
        if self.mol_type == 'DNA':
            self.complement_pairs = PAIRS_DNA
        elif self.mol_type == 'RNA':
            self.complement_pairs = PAIRS_RNA
        else:
            raise TypeError('Incorrect type of the sequence')

    def complement(self):
        """
        Generates the complement of the sequence.

        Returns:
            NucleicAcidSequence: The complement sequence.
        """
        res = []
        for base in self.seq:
            res.append(self.complement_pairs[base] if base.islower() else self.complement_pairs[base.lower()].upper())
        return NucleicAcidSequence(''.join(res))

    def gc_content(self):
        """
        Calculates the GC content of the sequence.

        Returns:
            float: The GC content of the sequence.
        """
        cnt = Counter(self.seq.upper())
        return round((cnt['C'] + cnt['G']) / len(self.seq), 4)


class DNASequence(NucleicAcidSequence):
    """
    Represents a DNA sequence.

    Inherits from NucleicAcidSequence class.

    Attributes:
        seq (str): The DNA sequence.
        mol_type (str): The type of molecule in the sequence (DNA).
        complement_pairs (dict): Dictionary containing complement pairs for DNA.
    """

    def __init__(self, seq):
        super().__init__(seq)
        if not self.mol_type == 'DNA':
            raise TypeError(f'Sequence {self.seq} is not DNA')

    def transcribe(self):
        """
        Transcribes the DNA sequence into an RNA sequence.

        Returns:
            RNASequence: The transcribed RNA sequence.
        """
        res = self.seq.replace('T', 'U').replace('t', 'u')
        return RNASequence(res)


class RNASequence(NucleicAcidSequence):
    """
    Represents an RNA sequence.

    Inherits from NucleicAcidSequence class.

    Attributes:
        seq (str): The RNA sequence.
        mol_type (str): The type of molecule in the sequence (RNA).
        complement_pairs (dict): Dictionary containing complement pairs for RNA.

    """

    def __init__(self, seq):
        super().__init__(seq)
        if not self.mol_type == 'RNA':
            raise TypeError(f'Sequence {self.seq} is not RNA')

    def reverse_transcribe(self):
        """
        Reverse transcribes the RNA sequence into a DNA sequence.

        Returns:
            DNASequence: The reverse transcribed DNA sequence.
        """
        res = self.seq.replace('U', 'T').replace('u', 't')
        return DNASequence(res)


class AminoAcidSequence(TypeDefinedSequence):
    """
    Represents an amino acid sequence.

    Inherits from TypeDefinedSequence class.

    Attributes:
        seq (str): The amino acid sequence.
        mol_type (str): The type of molecule in the sequence (AA_seq).
    """

    def __init__(self, seq):
        super().__init__(seq)
        if not self.mol_type == 'AA_seq':
            raise TypeError(f'Sequence {self.seq} is not protein')

    def calculate_pi(self):
        """
        Calculates the isoelectric point (pI) of the amino acid sequence.

        Returns:
            float: The isoelectric point (pI) of the sequence.
        """
        seq_list = list(self.seq.strip())
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


class ConnectionError(Exception):
    pass


@dataclass
class GenscanOutput:
    status: str
    cds_list: list
    exon_list: list


def run_genscan(sequence=None, sequence_file=None, organism="Vertebrate", exon_cutoff=1.00, sequence_name=""):
    site_url = 'http://hollywood.mit.edu'

    form_data = {
        "-o": organism,
        "-e": exon_cutoff,
        "-n": sequence_name,
        "-p": "Predicted peptides only"
    }
    if not sequence:
        with open(sequence_file, mode='rb') as file:
            files = {"-u": file}
            response = requests.post(site_url + "/cgi-bin/genscanw_py.cgi", data=form_data, files=files)
    else:
        form_data["-s"] = sequence
        response = requests.post(site_url + "/cgi-bin/genscanw_py.cgi", data=form_data)

    if not response.status_code == 200:
        raise ConnectionError('Check Internet Connection')

    soup = BeautifulSoup(response.text, "html.parser")

    status = response.status_code


    predictions = soup.find('pre').get_text()  # predictions = soup.find('pre').get_text().split('\n'*6)
    ind = predictions.find('Predicted peptide sequence(s):')
    cds = predictions[ind + 30:].strip().split('>')
    cds_list = []
    for i in cds[1:]:
        cds_list.append(''.join(i.split("\n\n")[1:]))

    predictions_splitted = predictions.split('\n' * 6)
    ind_1st_tab = predictions_splitted[0].find('Predicted genes/exons:')
    exones = predictions_splitted[0][ind_1st_tab + len('Predicted genes/exons:'):].strip().split('\n' * 4)[1:]
    n = 0
    exon_list = []
    for i in exones:
        for j in i.strip().split('\n' * 2):
            curr = j.strip().split()
            exon_list.append({'Number': n, 'Type': curr[1], 'Start': int(curr[3]), 'End': int(curr[4])})
            n += 1

    return GenscanOutput(status, cds_list, exon_list)


@dataclass
class FastaRecord:
    id: str
    seq: str
    description: str = ""


class OpenFasta:
    def __init__(self, name):
        self.name = name
        self.file = None

    def __enter__(self):
        self.file = open(self.name)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file:
            self.file.close()

    def __iter__(self):
        return self

    def __next__(self):
        record = self.read_record()

        for rec in record:
            if rec is None:
                raise StopIteration
            return rec


    def read_record(self):
        pattern = r'(?<=\>)[A-Z\d\.]+?(?=\s)'
        seq_id = None
        seq_description = ''
        curr_seq = ''


        for line in self.file:
            line = line.strip()
            match = re.search(pattern, line)
            if match is not None:
                if seq_id is not None:
                    res = FastaRecord(seq_id, curr_seq, seq_description)

                    yield res

                seq_id = match[0]
                seq_description = line[match.end()+1:]
                curr_seq = ''

            else:
                curr_seq += line
        yield FastaRecord(seq_id, curr_seq, seq_description) if seq_id is not None else None

    def read_records(self):
        recs = []
        rec = self.read_record()
        for i in rec:
            recs.append(i)
        return recs
