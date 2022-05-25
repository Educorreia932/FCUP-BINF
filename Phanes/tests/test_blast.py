import unittest


class TestAlignment(unittest.TestCase):
    def test(self):
        mb = Blast("seqBlast.txt", 11)
        query2 = "cgacgacgacgacgaatgatg"
        r = mb.bestAlignment(query2)

        print(r)
