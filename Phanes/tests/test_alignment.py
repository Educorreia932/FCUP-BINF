import unittest

from Phanes.src.phanes import DNA, AlignedSequences


class TestAlignment(unittest.TestCase):
    def setUp(self):
        pass

    def test_consensus(self):
        aligned = AlignedSequences(["ATAGC", "A_ACC"])
        consensus = "ATACC"

        self.assertEqual(consensus, aligned.consensus())


if __name__ == "__main__":
    unittest.main()
