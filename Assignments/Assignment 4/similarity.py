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

# 1)

score_matrix = []

for i, p1 in enumerate(list(proteins.items())):
    score_matrix.append([])

    # Add elements that were already calculated
    # We only calculate them once, since similarity between A and B is the same as B and A
    for j in range(0, i):
        score_matrix[-1].append(score_matrix[j][i])

    for p2 in list(proteins.items())[i:]:
        needleman_wunsch = NeedlemanWunsch(p1[1], p2[1], -8)
        needleman_wunsch.calculate()

        score_matrix[-1].append((p2[0], needleman_wunsch.alignment_score, needleman_wunsch.recover_align()))

scores_list = list(sorted(score_matrix[0], key=lambda item: item[1], reverse=True))

print(f"Exercise 1:")
print("\t<Protein>: <Score>")

for score in scores_list:
    print(f"\t{score[0]}: {score[1]}")

# 2)

# 3)

mismatch_matrix = []

for i, _ in enumerate(score_matrix):
    mismatch_matrix.append([])

    for p2 in score_matrix[i]:
        mismatch_matrix[-1].append(count_mismatches(*p2[2]))

def plot_matrix(matrix):
    print("{}".format("".join('\t' + str(i) for i in range(len(matrix[0])))))

    for i, p1 in enumerate(matrix):
        print(f"{i}\t", end="")

        for j, p2 in enumerate(matrix[i]):
            print(f"{matrix[i][j]}\t", end="")

        print()

print()
print("Exercise 3:")

print()
print("Scores matrix:")
print()

plot_matrix([[p[1] for p in row] for row in score_matrix])

print()
print("Mismatches matrix:")
print()

plot_matrix(mismatch_matrix)

# 4)

# 5)
