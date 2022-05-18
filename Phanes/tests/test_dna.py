import unittest

from Phanes.src.phanes import DNA


class TestDNA(unittest.TestCase):
    def setUp(self):
        self.sequence = DNA("ATGAAATTATGAATGAGCCTCAGCTGAAGCATCGCGCATCAGACTACGCTCAGACTCAGACTCAGCATTATAGTGAATGTTAATAAATAAAATAA")

    def test_reverse_complement(self):
        reverse_complement = DNA("TTATTTTATTTATTAACATTCACTATAATGCTGAGTCTGAGTCTGAGCGTAGTCTGATGCGCGATGCTTCAGCTGAGGCTCATTCATAATTTCAT")

        self.assertEqual(reverse_complement, self.sequence)

    # def test(self):
    #     print(f"Sequence: {dna}")
    #     print(f"Length: {len(dna)}")
    #     print(f"RNA: {dna.transcribe()}")
    #     print(f"Reverse complement: {dna.reverse_complement()}")
    #
    #     self.assertEqual((5), 6)


if __name__ == "__main__":
    unittest.main()
