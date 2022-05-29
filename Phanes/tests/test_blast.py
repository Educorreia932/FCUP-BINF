import unittest

from phanes import Blast


class TestAlignment(unittest.TestCase):
    def test(self):
        blast = Blast("database.txt", 11)
        query = "cgacgacgacgacgaatgatg"

        self.assertTupleEqual(blast.calculate(query), (0, 0, 21, 21, 4))
