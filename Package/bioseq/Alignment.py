
class Alignment:
    def __init__(self, seq1: str, seq2: str, gap: int, submat_file="blosum62.mat"):
        self.seq1 = seq1
        self.seq2 = seq2
        self.gap = gap
        self.read_submat_file(submat_file)

    def read_submat_file(self, filename):
        """Read substitution matrix from file """
        
        self.sm = {}

        with open(filename, "r") as f:
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
                    self.sm[k] = int(tokens[j])
        
    def score_position(self, c1, c2):
        """Score of a position (column)"""

        if c1 == "-" or c2 == "-":
            return self.gap

        else:
            return self.sm[c1 + c2]

    def max3t(self, v1, v2, v3):
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

    # Global alignment
    def needleman_wunsch(self):
        S = [[0]]
        T = [[0]]

        # Initialize gaps in rows
        for j in range(1, len(self.seq2) + 1):
            S[0].append(self.gap * j)
            T[0].append(3)  # horizontal move: 3
        
        # Initialize gaps in cols
        for i in range(1, len(self.seq1) + 1):
            S.append([self.gap * i])
            T.append([2])  # vertical move: 2
        
        # Apply the recurrence to fill the matrices
        for i in range(0, len(self.seq1)):
            for j in range(len(self.seq2)):
                s1 = S[i][j] + self.score_position(self.seq1[i], self.seq2[j])  # Diagonal
                s2 = S[i][j + 1] + self.gap                                                        # Vertical
                s3 = S[i + 1][j] + self.gap                                                        # Horizontal

                S[i + 1].append(max(s1, s2, s3))  # na matrix score add max value
                T[i + 1].append(self.max3t(s1, s2, s3))

        return (S, T)

    # Local alignment
    def smith_waterman(self):
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
        
        return (S, T, maxscore)

