"""
Bruno Vaz - up201705247
Eduardo Correia - up201806433
Filipe Justiça - up201606339 
"""

from bioseq.BioSequence import BioSequence
from bioseq.Alignment import NeedlemanWunsch, SmithWaterman

proteins = BioSequence.read_fasta("glycoproteinS.fas")

seq1 = proteins['YP_009724390.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]']
seq2 = proteins['QHO60594.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]']


def count_mismatches(s1, s2):
    mismatches = 0

    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            mismatches += 1

    return mismatches


def calculate_score_matrix(proteins, algorithm):
    score_matrix = []

    for i, p1 in enumerate(list(proteins.items())):  # Para alterar!!
        score_matrix.append([])

        # Add elements that were already calculated
        # We only calculate them once, since similarity between A and B is the same as B and A
        for j in range(0, i):
            score_matrix[-1].append(score_matrix[j][i])

        for p2 in list(proteins.items())[i:]:  # Para alterar
            alignment = algorithm(p1[1], p2[1], -8)
            alignment.calculate()

            score_matrix[-1].append((p2[0], alignment.alignment_score, alignment.recover_align()))

    return score_matrix


def calculate_mismatch_matrix(score_matrix):
    mismatch_matrix = []

    for i, _ in enumerate(score_matrix):
        mismatch_matrix.append([])

        for p2 in score_matrix[i]:
            mismatch_matrix[-1].append(count_mismatches(*p2[2]))

    return mismatch_matrix


def plot_matrix(matrix, protein_names):
    print("\t\t\t  {}".format(protein_names[0]), end="")
    print("%16s" * 9 % tuple(protein_names[1:]))

    for i, p1 in enumerate(matrix):
        print("{:>22}\t".format(protein_names[i]), end="")

        for j, p2 in enumerate(matrix[i]):
            print("%16s" % matrix[i][j], end="")

        print()

# 1)


score_matrix = calculate_score_matrix(proteins, NeedlemanWunsch)
scores_list = list(sorted(score_matrix[0], key=lambda item: item[1], reverse=True))

print(f"Exercise 1:")
print("\t<Protein>: <Score>")

for score in scores_list:
    print(f"\t{score[0]}: {score[1]}")

# 2)

# TODO:

# 3)

protein_names = [list(proteins.keys())[i].split(" ")[0] for i in range(len(proteins))]


print("\nExercise 3:")
print("\n\tScores matrix:\n")

plot_matrix([[p[1] for p in row] for row in score_matrix], protein_names)

print("\n\n\tMismatches matrix:\n")

mismatch_matrix = calculate_mismatch_matrix(score_matrix)
plot_matrix(mismatch_matrix, protein_names)

# 4)

# TODO:

# 5.1)

score_matrix = calculate_score_matrix(proteins, SmithWaterman)
scores_list = list(sorted(score_matrix[0], key=lambda item: item[1], reverse=True))

print(f"Exercise 5.1:")
print("\t<Protein>: <Score>")

for score in scores_list:
    print(f"\t{score[0]}: {score[1]}")

# 5.2)

# TODO:

# 5.3)

print("\nExercise 5.3:")
print("\n\tScores matrix:\n")

plot_matrix([[p[1] for p in row] for row in score_matrix], protein_names)

print("\n\n\tMismatches matrix:\n")

mismatch_matrix = calculate_mismatch_matrix(score_matrix)
plot_matrix(mismatch_matrix, protein_names)

# 5.4)

# TODO:
