import unittest

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
