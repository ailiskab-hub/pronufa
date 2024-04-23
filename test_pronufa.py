import unittest
import os
from pronufa import filter_fastq, BiologicalSequence, DNASequence, RNASequence, AminoAcidSequence
from bio_files_processor import convert_multiline_fasta_to_oneline, select_genes_from_gbk_to_fasta, change_fasta_start_pos


class TestModulePronufa(unittest.TestCase):
    def test_filter_fastq_file_io(self, test_file_path='data/test_file.fastq'):
        output_file_path = "filtered_fastq"
        filter_fastq(test_file_path, output_filename=output_file_path)

        self.assertTrue(os.path.isfile(f'fastq_filtrator_resuls/{output_file_path}.fastq'))

        os.remove(f'fastq_filtrator_resuls/{output_file_path}.fastq')

    def test_check_type(self):
        dna_seq = BiologicalSequence("TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG")
        self.assertTrue(dna_seq.check_type('DNA'))

        rna_seq = BiologicalSequence("UACUUUUUUGGCAGCUGUAUCGCCGCUAGC")
        self.assertTrue(rna_seq.check_type('RNA'))

        aa_seq = BiologicalSequence("MTLPVSGVILNADLSELIGTKEELSNIDLV")
        self.assertTrue(aa_seq.check_type('AA_seq'))

        # check mistake
        rna_seq = BiologicalSequence("UACC")
        self.assertTrue(rna_seq.check_type('DNA'))

    def test_transcribe_reverse_transcribe(self):
        dna_seq = DNASequence("AGCT")
        rna_seq = dna_seq.transcribe()
        self.assertEqual(str(rna_seq), "AGCU")

        rna_seq = RNASequence("AGCU")
        dna_seq = rna_seq.reverse_transcribe()
        self.assertEqual(str(dna_seq), "AGCT")

    def test_gc_content(self):
        dna_seq = DNASequence("AGCT")
        self.assertAlmostEqual(dna_seq.gc_content(), 0.5)

        rna_seq = RNASequence("AGCU")
        self.assertEqual(rna_seq.gc_content(), 0.5)

    def test_calculate_pi(self):
        aa_seq = AminoAcidSequence("LPTQQ")
        self.assertAlmostEqual(aa_seq.calculate_pi(), 5.885, places=3)


class TestModuleBioprocessor(unittest.TestCase):

    def test_convert_multiline_fasta_to_oneline(self, test_file_path='data/example_gbk.GBK'):
        output_file_path = "converted_oneline"

        convert_multiline_fasta_to_oneline(test_file_path, output_fasta=output_file_path)
        self.assertTrue(os.path.isfile(f'{output_file_path}.fasta'))

        os.remove(f'{output_file_path}.fasta')

    def test_select_genes_from_gbk_to_fasta(self, test_file_path='data/example_gbk.GBK'):
        genes_of_interest = ["sdhD"]
        output_file_path = "selected_genes"

        select_genes_from_gbk_to_fasta(test_file_path, genes_of_interest, output_fasta=output_file_path)
        self.assertTrue(os.path.isfile(f'{output_file_path}.fasta'))

        os.remove(f'{output_file_path}.fasta')

    def test_change_fasta_start_pos(self, test_file_path='data/example_fasta.fasta'):

        output_file_path = "changed_start_pos"

        change_fasta_start_pos(test_file_path, shift=1, output_fasta=output_file_path)
        self.assertTrue(os.path.isfile(f'{output_file_path}.fasta'))

        os.remove(f'{output_file_path}.fasta')


if __name__ == '__main__':
    unittest.main()

