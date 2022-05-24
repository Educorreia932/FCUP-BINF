import abc
import os

from .AlignedSequences import AlignedSequences


class PairwiseAlignment:
    def __init__(self, seq1: str = None, seq2: str = None, gap: int = -1, sm=None, submat_file="blosum62.mat"):
        self.sm = None
        self.seq1 = seq1
        self.seq2 = seq2
        self.gap = gap
        self.S = None
        self.T = None
        self._score = None

        if sm is None:
            self.sm = self.read_submat_file(submat_file)

        else:
            self.sm = sm

    def score_position(self, c1, c2):
        """Score of a position (column)"""

        if c1 == "-" or c2 == "-":
            return self.gap

        else:
            return self.sm[c1 + c2]

    @staticmethod
    def read_submat_file(filename: str) -> dict:
        """Read substitution matrix from file"""

        sm = {}

        filepath = os.path.join(os.path.dirname(__file__), filename)

        with open(filepath, "r") as f:
            line = f.readline()
            tokens = line.split("\t")
            ns = len(tokens)
            alphabet = []

            for i in range(0, ns):
                alphabet.append(tokens[i][0])

            for i in range(0, ns):
                line = f.readline()
                tokens = line.split("\t")

                for j in range(0, len(tokens)):
                    k = alphabet[i] + alphabet[j]
                    sm[k] = int(tokens[j])

        return sm

    @staticmethod
    def create_submatrix(match: int, mismatch: int, alphabet: str) -> dict:
        sm = {}

        for c1 in alphabet:
            for c2 in alphabet:
                if c1 == c2:
                    sm[c1 + c2] = match

                else:
                    sm[c1 + c2] = mismatch

        return sm

    @staticmethod
    def max3t(v1, v2, v3):
        """Provides the integer to fill in T"""

        if v1 > v2:
            if v1 > v3:
                return 1

            else:
                return 3

        else:
            if v2 > v3:
                return 2

            else:
                return 3

    @property
    def alignment_score(self):
        return self._score

    def set_sequences(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2

    @abc.abstractmethod
    def recover_align(self) -> AlignedSequences:
        pass

    @abc.abstractmethod
    def calculate(self):
        pass


# Global alignment
class NeedlemanWunsch(PairwiseAlignment):
    def calculate(self):
        self.S = [[0]]
        self.T = [[0]]

        # Initialize gaps in rows
        for j in range(1, len(self.seq2) + 1):
            self.S[0].append(self.gap * j)
            self.T[0].append(3)  # horizontal move: 3

        # Initialize gaps in cols
        for i in range(1, len(self.seq1) + 1):
            self.S.append([self.gap * i])
            self.T.append([2])  # vertical move: 2

        # Apply the recurrence to fill the matrices
        for i in range(0, len(self.seq1)):
            for j in range(len(self.seq2)):
                s1 = self.S[i][j] + self.score_position(self.seq1[i], self.seq2[j])  # Diagonal
                s2 = self.S[i][j + 1] + self.gap  # Vertical
                s3 = self.S[i + 1][j] + self.gap  # Horizontal

                self.S[i + 1].append(max(s1, s2, s3))  # na matrix score add max value
                self.T[i + 1].append(self.max3t(s1, s2, s3))

        self._score = self.S[-1][-1]

    def recover_align(self) -> AlignedSequences:
        # alignment are two strings
        res = ["", ""]
        i = len(self.seq1)
        j = len(self.seq2)

        while i > 0 or j > 0:
            # Diagonal move
            if self.T[i][j] == 1:
                res[0] = self.seq1[i - 1] + res[0]  # add to align of seq1 a symbol from seq1(i-1)
                res[1] = self.seq2[j - 1] + res[1]  # add to align of seq2 a symbol from seq2(i-1)
                i -= 1
                j -= 1

            # Horizontal move
            elif self.T[i][j] == 3:
                res[0] = "-" + res[0]  # insert gap na seq 1
                res[1] = self.seq2[j - 1] + res[1]  # insert symbol from seq2
                j -= 1

            # Vertical move
            else:
                res[0] = self.seq1[i - 1] + res[0]  # insert symbol from seq1
                res[1] = "-" + res[1]  # insert gap na seq 2
                i -= 1

        return AlignedSequences(res)


# Local alignment
class SmithWaterman(PairwiseAlignment):
    def calculate(self):
        S = [[0]]
        T = [[0]]
        maxscore = 0

        # First row filled with zeros
        for j in range(1, len(self.seq2) + 1):
            S[0].append(0)
            T[0].append(0)

        # First column filled with zeros
        for i in range(1, len(self.seq1) + 1):
            S.append([0])
            T.append([0])

        for i in range(0, len(self.seq1)):
            for j in range(len(self.seq2)):
                s1 = S[i][j] + self.score_position(self.seq1[i], self.seq2[j])
                s2 = S[i][j + 1] + self.gap
                s3 = S[i + 1][j] + self.gap
                b = max(s1, s2, s3)

                if b <= 0:
                    S[i + 1].append(0)
                    T[i + 1].append(0)

                else:
                    S[i + 1].append(b)
                    T[i + 1].append(self.max3t(s1, s2, s3))

                    if b > maxscore:
                        maxscore = b
                        maxmat = (i, j)

        self.S = S
        self.T = T
        self._score = maxscore
        self.maxmat = maxmat

        return S, T, maxscore

    def recover_align(self):
        """Recover one of the optimal alignments"""
        res = ["", ""]

        # Determine the cell with max score
        i, j = self.maxmat

        # Terminates when finds a cell with zero
        while self.T[i][j] > 0:
            if self.T[i][j] == 1:
                res[0] = self.seq1[i - 1] + res[0]
                res[1] = self.seq2[j - 1] + res[1]
                i -= 1
                j -= 1

            elif self.T[i][j] == 3:
                res[0] = "-" + res[0]
                res[1] = self.seq2[j - 1] + res[1]
                j -= 1

            elif self.T[i][j] == 2:
                res[0] = self.seq1[i - 1] + res[0]
                res[1] = "-" + res[1]
                i -= 1

        return res
