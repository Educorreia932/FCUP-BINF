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

needleman_wunsch = NeedlemanWunsch(seq1, seq2, -8)
needleman_wunsch.calculate()
s1, s2 = needleman_wunsch.recover_align()

print(f"Exercise 1:")
print(f"\tThe number of mismatches in global alignment is {count_mismatches(s1, s2)}.")
print()

# 2)

smith_waterman = SmithWaterman(seq1, seq2, -8)
smith_waterman.calculate()
s1, s2 = smith_waterman.recover_align()

print(f"Exercise 2:")
print(f"\tThe number of mismatches in local alignment is {count_mismatches(s1, s2)}.")
print()

# 3)
ref_protein = seq1
protas_score = {}
for prota in list(proteins.keys())[1:]:
    seq = proteins[prota]

    needleman_wunsch = NeedlemanWunsch(ref_protein, seq, -8)
    needleman_wunsch.calculate()
    protas_score[prota] = needleman_wunsch.alignment_score

scores_list = list(
    sorted(protas_score.items(), key=lambda item: item[1], reverse=True)
)

print(f"Exercise 3:")
print("\t Protein: Score")
for tuplo in scores_list:
    print(f"\t {tuplo[0]}: {tuplo[1]}")



# 4)

