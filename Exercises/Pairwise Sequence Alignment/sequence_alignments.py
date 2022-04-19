def create_submat(match, mismatch, alphabet):
    """ substitution matrix as dictionary """
    sm = {}

    for c1 in alphabet:
        for c2 in alphabet:
            if (c1 == c2):
                sm[c1 + c2] = match

            else:
                sm[c1 + c2] = mismatch

    return sm


def read_submat_file(filename):
    """read substitution matrix from file """
    sm = {}
    f = open(filename, "r")
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


def score_pos(c1, c2, sm, g):
    """Score of a position (column)"""

    if c1 == "-" or c2 == "-":
        return g

    else:
        return sm[c1 + c2]

# score of the whole alignment


def score_align(seq1, seq2, sm, g):
    """
    Score of the whole alignment; iterate through the two sequences
    sum the score of each position and return its sum; assume sequences are of equal length
    """

    score = 0

    for c1 in seq1:
        for c2 in seq2:
            score += score_pos(c1, c2, sm, g)
    
    return score


def score_affinegap(seq1, seq2, sm, g, r):
    ''' calculates the score of alignment based on : affine_gap(len) = g + r*len
    if the gap is open (first occurrence) sum value g; if gap continues sum r to each new gap position;
    if there is no gap use the substitution matrix for the score.
    '''
    res = 0
    ingap1 = False  # two f are true when inside gap sequences
    ingap2 = False
    for i in range(len(seq1)):
        if seq1[i] == "-":
            # gap is already open; add r
            if ingap1:
                res += r
            else:
                # gap is open for the first time; add g
                ingap1 = True
                res += g
        elif seq2[i] == "-":
            # gap is already open; add r
            if ingap2:
                res += r
            else:
                # gap is open for the first time; add g
                ingap2 = True
                res += g
        else:
            # no gaps; use substitution matrix
            if ingap1:
                ingap1 = False
            if ingap2:
                ingap2 = False
            res += sm[seq1[i] + seq2[i]]
    return res

## global alignment


def needleman_Wunsch(seq1, seq2, sm, g):
    """Global Alignment"""
    S = [[0]]
    T = [[0]]

    # Initialize gaps in rows
    for j in range(1, len(seq2) + 1):
        S[0].append(g * j)
        T[0].append(3)  # horizontal move: 3
    
    # Initialize gaps in cols
    for i in range(1, len(seq1) + 1):
        S.append([g * i])
        T.append([2])  # vertical move: 2
    
    # Apply the recurrence to fill the matrices
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos(seq1[i], seq2[j], sm, g)  # Diagonal
            s2 = S[i][j + 1] + g                               # Vertical
            s3 = S[i + 1][j] + g                               # Horizontal

            S[i + 1].append(max(s1, s2, s3))  # na matrix score add max value
            T[i + 1].append(max3t(s1, s2, s3))

    return (S, T)


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


def recover_align(T, seq1, seq2):
    # alignment are two strings
    res = ["", ""]
    i = len(seq1)
    j = len(seq2)

    while i > 0 or j > 0:
        # Diagonal move
        if T[i][j] == 1:    
            res[0] = seq1[i - 1] + res[0]  # add to align of seq1 a symbol from seq1(i-1)
            res[1] = seq2[j - 1] + res[1]  # add to align of seq2 a symbol from seq2(i-1)
            i -= 1
            j -= 1

        # Horizontal move
        elif T[i][j] == 3:  
            res[0] = "-" + res[0]   # insert gap na seq 1
            res[1] = seq2[j - 1] + res[1]  # insert symbol from seq2
            j -= 1

        # Vertical move
        else:               
            res[0] = seq1[i - 1] + res[0]  # insert symbol from seq1
            res[1] = "-" + res[1]  # insert gap na seq 2
            i -= 1
    return res

# local alignment


def smith_Waterman(seq1, seq2, sm, g):
    """Local alignment"""
    S = [[0]]
    T = [[0]]
    maxscore = 0
    # first row filled with zero
    for j in range(1, len(seq2) + 1):
        S[0].append(0)
        T[0].append(0)
    # first column filled with zero
    for i in range(1, len(seq1) + 1):
        S.append([0])
        T.append([0])
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos(seq1[i], seq2[j], sm, g)
            s2 = S[i][j + 1] + g
            s3 = S[i + 1][j] + g
            b = max(s1, s2, s3)
            if b <= 0:
                S[i + 1].append(0)
                T[i + 1].append(0)
            else:
                S[i + 1].append(b)
                T[i + 1].append(max3t(s1, s2, s3))
                if b > maxscore:
                    maxscore = b
    return (S, T, maxscore)


def recover_align_local(S, T, seq1, seq2):
    """recover one of the optimal alignments"""
    res = ["", ""]
    """determine the cell with max score"""
    i, j = max_mat(S)
    """terminates when finds a cell with zero"""
    while T[i][j] > 0:
        if T[i][j] == 1:
            res[0] = seq1[i - 1] + res[0]
            res[1] = seq2[j - 1] + res[1]
            i -= 1
            j -= 1
        elif T[i][j] == 3:
            res[0] = "-" + res[0]
            res[1] = seq2[j - 1] + res[1]
            j -= 1
        elif T[i][j] == 2:
            res[0] = seq1[i - 1] + res[0]
            res[1] = "-" + res[1]
            i -= 1
    return res


def max_mat(matrix):
    """Finds the max cell in the matrix"""
    result = tuple()

    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if len(result) == 0 or matrix[i][j] > matrix[result[0]][result[1]]:
                result = (i, j)

    return result


def identity(seq1, seq2, alphabet="ACGT"):
    """Calculate the identity score between two sequences"""

    score = 0

    for c1 in seq1:
        for c2 in seq2:
            if c1 == c2:
                sccore += 1

    return score


def print_mat(mat):
    for i in range(0, len(mat)):
        print(mat[i])


def test_DNA():
    sm = create_submat(2, -2, "ACGT")
    seq1 = "-CAGTGCATG-ACATA"
    seq2 = "TCAG-GC-TCTACAGA"
    g = -3

    print(score_align(seq1, seq2, sm, g))


def test_prot():
    # test the alignment of the two sequence
    # plot the score of alignment using the subtitution matrix blosum62.mat
    # plot the score of alignment using affine gap score with gap value.
    # complete here ...
    pass


def test_global_alig():
    sm = read_submat_file("blosum62.mat")
    seq1 = "PHSWG"      
    seq2 = "HGWAG"      
    res = needleman_Wunsch(seq1, seq2, sm, -2)  
    S = res[0]  # Scores matrix
    T = res[1]  # Trace-back matrix
    print("Score of optimal alignment:", S[len(seq1)][len(seq2)])  # cell lower right
    print_mat(S)
    print_mat(T)
    alig = recover_align(T, seq1, seq2)
    print(alig[0])
    print(alig[1])


def test_local_alig():
    sm = read_submat_file("blosum62.mat")
    seq1 = "PHSWG"
    seq2 = "HGWAG"
    res = smith_Waterman(seq1, seq2, sm, -8)
    S = res[0]
    T = res[1]
    print("Score of optimal alignment:", res[2])
    print_mat(S)
    print_mat(T)
    alinL = recover_align_local(S, T, seq1, seq2)
    print(alinL[0])
    print(alinL[1])
    i, j = max_mat(S)  # from Score matrix get the indices i,j from cell with max value
    best_score = S[i][j]  # best score of the alignment from cell i,j
    print("best score: " + str(best_score))


def test_DNA_GlobalAlign():
    # test function
    # Test sequences seq1 and seq2
    # create a substitution matrix with the match and mismatch values
    # solve the NW algorithm with gap p
    # obtain the score of the alignment: using matrix cells and score alignment function
    # recover the alignment and print the aligned sequences 1 and 2

    # complete here ...
    pass


def test_Prot_LocalAlign():
    # Test local alignment SW to sequences seq1 and seq2
    pass


def exam_local_alig():
    sm = create_submat(2, 0, "MCSNHL")
    seq1 = "MCSNH"
    seq2 = "SDHL"
    res = smith_Waterman(seq1, seq2, sm, -2)
    S = res[0]
    T = res[1]
    print("Score of optimal alignment:", res[2])
    print_mat(S)
    print_mat(T)
    alinL = recover_align_local(S, T, seq1, seq2)
    print(alinL[0])
    print(alinL[1])
    i, j = max_mat(S)  # from Score matrix get the indices i,j from cell with max value
    best_score = S[i][j]  # best score of the alignment from cell i,j
    print("best score: " + str(best_score))


if __name__ == "__main__":
    test_DNA()
    test_prot()
    test_global_alig()
    test_local_alig()
