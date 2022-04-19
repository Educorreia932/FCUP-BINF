import matplotlib.pyplot as plt

def create_matrix(nrows, ncols):
    return [[0 for j in range(ncols)] for i in range(nrows)]


def dotplot(seq1, seq2):
    """Create a matrix based on the the input sequences and fill the cells that correspond to a match"""

    matrix = create_matrix(len(seq1), len(seq2))

    for i, c1 in enumerate(seq1):
        for j, c2 in enumerate(seq2):
            if c1 == c2:
                matrix[i][j] = 1

    return matrix


def extended_dotplot(seq1, seq2, window, stringency):
    mat = create_matrix(len(seq1), len(seq2))
    start = int(window / 2)

    for i in range(start, len(seq1) - start):
        for j in range(start, len(seq2) - start):
            matches = 0
            l = j - start

            for k in range(i - start, i + start + 1):
                if seq1[k] == seq2[l]:
                    matches += 1

                l += 1

                if matches >= stringency:
                    mat[i][j] = 1

    return mat


def print_dotplot(matrix, seq1, seq2):
    seq1 = " " + seq1
    seq2 = " " + seq2

    for i in range(len(matrix) + 1):
        for j in range(len(matrix[0]) + 1):
            if i == 0:
                print(seq2[j], end="")

            else:
                if j == 0:
                    print(seq1[i], end="")

                else:
                    if matrix[i - 1][j - 1] == 1:
                        print("*", end="")

                    else:
                        print(" ", end=" ")

        print()


def dotplot_chart(matrix):
    x = [i for i in range(len(matrix)) for j in range(len(matrix[0])) if matrix[i][j] == 1]  
    y = [j for i in range(len(matrix)) for j in range(len(matrix[0])) if matrix[i][j] == 1]  

    plt.scatter(x, y)
    plt.show()

def test_diagonal_length(mat, istart, jstart):
    # given the starting indices on the row and column
    # check along the diagonal that starts in istart and jstart
    # the longest sub-sequences of matches; return this value
    # ....
    pass


def test():
    HBA = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
    HBB = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"

    matrix = extended_dotplot(HBA, HBB, 10, 4)
    dotplot_chart(matrix)


if __name__ == "__main__":
    test()
