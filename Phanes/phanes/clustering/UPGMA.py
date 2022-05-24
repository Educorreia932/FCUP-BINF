from alignment import PairwiseAlignment
from . import HierarchicalClustering, Matrix, BinaryTree


class UPGMA:
    def __init__(self, sequences, alignment: PairwiseAlignment):
        self.sequences = sequences
        self.alignment = alignment

    def _calculate_distance_matrix(self) -> Matrix:
        distance_matrix = Matrix(len(self.sequences))

        for i, s1 in enumerate(self.sequences):
            for j, s2 in enumerate(self.sequences[i + 1:]):
                j += i + 1

                self.alignment.set_sequences(s1, s2)
                self.alignment.calculate()

                aligned = self.alignment.recover_align()
                distance = sum(1 if a != b else 0 for a, b in zip(*aligned))
                distance_matrix[i, j] = distance

        return distance_matrix

    def calculate(self) -> BinaryTree:
        return HierarchicalClustering(self._calculate_distance_matrix()).calculate()
