from copy import copy

from . import BinaryTree, Matrix


class HierarchicalClustering:
    def __init__(self, distance_matrix: Matrix):
        self.distance_matrix = copy(distance_matrix)

    def calculate(self) -> BinaryTree:
        """
        Returns a tree based on distance matrix
        """

        # Create a tree for each sequence
        trees = [BinaryTree(-1) for _ in range(self.distance_matrix.size()[0])]
        iterations = self.distance_matrix.size()[0]

        for k in range(iterations, 1, -1):
            minimum_distance = -1
            minimum_distance_indexes = (0, 0)

            # Find indices of the minimum distance in the matrix
            for i in range(1, self.distance_matrix.size()[0]):
                for j in range(i):
                    if self.distance_matrix[i, j] < minimum_distance or minimum_distance == -1:
                        minimum_distance_indexes = (i, j)
                        minimum_distance = self.distance_matrix[i, j]

            # Create a new tree to join the clusters with minimum distance
            i, j = minimum_distance_indexes
            tree = BinaryTree(k, distance=minimum_distance / 2, left=trees[i], right=trees[j])

            if k > 2:
                # Remove from the list of trees the joined branches
                tree_i = trees.pop(i)
                tree_j = trees.pop(j)

                distances = []

                # Calculate the distances for the new cluster
                for x in range(self.distance_matrix.size()[0]):
                    if x != i and x != j:
                        si = len(tree_i.leaves())  # |si|
                        sj = len(tree_j.leaves())  # |sj|

                        # Use the weighted average to calculate the distances between the clustersd
                        distance = (si * self.distance_matrix[i, x] + sj * self.distance_matrix[j, x]) / (si + sj)
                        distances.append(distance)

                # Remove column corresponding to j
                self.distance_matrix.remove_column(i)
                self.distance_matrix.remove_column(j)

                # Remove row corresponding to i
                self.distance_matrix.remove_row(i)
                self.distance_matrix.remove_row(j)

                # Add row with new distances
                self.distance_matrix.add_row(distances)

                # Add column with zero distances: of len (|dists| + 1)
                self.distance_matrix.add_column([0 for _ in range(len(distances) + 1)])

                trees.append(tree)

            else:
                return tree
