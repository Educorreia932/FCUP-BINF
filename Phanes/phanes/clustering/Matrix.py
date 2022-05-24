class Matrix:
    def __init__(self, rows):
        if type(rows) == int:
            self.matrix = [[0 for _ in range(rows)] for _ in range(rows)]

        else:
            self.matrix = [row for row in rows]

    def add_row(self, row: list):
        self.matrix.append(row)

    def add_column(self, column: list):
        for i, row in enumerate(self.matrix):
            row.append(column[i])

    def remove_row(self, index: int):
        self.matrix.pop(index)

    def remove_column(self, index: int):
        for row in self.matrix:
            row.pop(index)

    def size(self):
        return len(self.matrix[0]), len(self.matrix)

    def __getitem__(self, position: [int]):
        i, j = position

        if i > j:
            return self.matrix[i][j]

        else:
            return self.matrix[j][i]

    def __setitem__(self, position: [int], value: int):
        i, j = position

        if i > j:
            self.matrix[i][j] = value

        else:
            self.matrix[j][i] = value
