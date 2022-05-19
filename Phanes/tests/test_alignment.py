import unittest

from phanes import NeedlemanWunsch, AlignedSequences, MultipleAlignment


class TestAlignment(unittest.TestCase):
    def setUp(self):
        pass

    def test_global_alignment(self):
        seq1 = "PHSWG"
        seq2 = "HGWAG"

        needleman_wunsch = NeedlemanWunsch(seq1, seq2, -8)
        needleman_wunsch.calculate()

        self.assertListEqual([-40, -24, -10, 3, 11, 9], needleman_wunsch.S[-1])

    def test_consensus(self):
        aligned = AlignedSequences(["ATAGC", "A-ACC"])
        consensus = "ATAGC"

        self.assertEqual(consensus, aligned.consensus())

    def test_multiple_alignment(self):
        sequences = ["ATAGC", "AACC", "ATGAC"]
        aligned = MultipleAlignment(sequences).calculate()

        self.assertEqual("AT-GAC", aligned[-1])


if __name__ == "__main__":
    unittest.main()
