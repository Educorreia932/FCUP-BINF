import unittest

from phanes import NeedlemanWunsch, AlignedSequences


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
        aligned = AlignedSequences(["ATAGC", "A_ACC"])
        consensus = "ATACC"

        self.assertEqual(consensus, aligned.consensus())


if __name__ == "__main__":
    unittest.main()
