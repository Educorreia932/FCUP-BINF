import unittest

from phanes import UPGMA
from phanes import PairwiseAlignment, NeedlemanWunsch
from phanes import BinaryTree, Matrix, HierarchicalClustering


class TestClustering(unittest.TestCase):
    def test_tree_leaves(self):
        # Leaves
        a = BinaryTree(1)
        b = BinaryTree(2)
        c = BinaryTree(3)
        d = BinaryTree(4)

        # Internal nodes
        e = BinaryTree(-1, 2.0, b, c)
        f = BinaryTree(-1, 1.5, d, a)
        g = BinaryTree(-1, 4.5, e, f)

        leaves = [leaf.value for leaf in g.leaves()]

        self.assertListEqual(leaves, [2, 3, 4, 1])

    def test_clustering(self):
        matrix = Matrix([
            [0, 0, 0, 0, 0],
            [3, 0, 0, 0, 0],
            [4, 6, 0, 0, 0],
            [4, 6, 7, 0, 0],
            [2, 5, 7, 9, 0],
        ])

        clustering = HierarchicalClustering(matrix)
        result = clustering.calculate()

        self.assertEqual(result.distance, 3.25)

    def test_upgma(self):
        sequences = [
            "ATAG-C-",
            "AT-GAC-",
            "A-A--CG",
            "A-AT-CG"
        ]

        substitution_matrix = PairwiseAlignment.create_submatrix(1, -1, "ACGT")
        alignment = NeedlemanWunsch(sm=substitution_matrix, gap=-2)
        upgma = UPGMA(sequences, alignment)

        result = upgma.calculate()

        self.assertEqual(result.distance, 2)
