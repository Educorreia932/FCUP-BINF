class BinaryTree:
    def __init__(self, value, distance=0, left=None, right=None):
        self.value = value
        self.distance = distance
        self.left = left
        self.right = right

    def leaves(self, leaves=None):
        if leaves is None:
            leaves = []

        # Root node
        if self.left is None and self.right is None:
            leaves.append(self)

            return leaves.copy()

        else:
            return self.left.leaves(leaves.copy()) + self.right.leaves(leaves.copy())
