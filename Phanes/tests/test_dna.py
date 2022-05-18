import unittest

from phanes.sequence import DNA, RNA


class TestDNA(unittest.TestCase):
    def setUp(self):
        self.sequence = DNA(
            "ATGAAATTATGAATGAGCCTCAGCTGAAGCATCGCGCATCAGACTACGCTCAGACTCAGACTCAGCATTATAGTGAATGTTAATAAATAAAATAA")

    def test_length(self):
        length = 95

        self.assertEqual(length, len(self.sequence))

    def test_reverse_complement(self):
        reverse_complement = DNA(
            "TTATTTTATTTATTAACATTCACTATAATGCTGAGTCTGAGTCTGAGCGTAGTCTGATGCGCGATGCTTCAGCTGAGGCTCATTCATAATTTCAT")

        self.assertEqual(reverse_complement, self.sequence.reverse_complement())

    def test_rna(self):
        rna = RNA("AUGAAAUUAUGAAUGAGCCUCAGCUGAAGCAUCGCGCAUCAGACUACGCUCAGACUCAGACUCAGCAUUAUAGUGAAUGUUAAUAAAUAAAAUAA")

        self.assertEqual(rna, self.sequence.transcribe())


if __name__ == "__main__":
    unittest.main()
