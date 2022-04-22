from bioseq.BioSequence import BioSequence
from bioseq.Alignment import NeedlemanWunsch, SmithWaterman

proteins = list(sequence for sequence in BioSequence.read_fasta("glycoproteinS.fas").values())

seq1 = proteins[0] # YP_009724390.1
seq2 = proteins[1] # QHO60594.1

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
s1, s2 = needleman_wunsch.recover_align()

print(f"Exercise 2:")
print(f"\tThe number of mismatches in local alignment is {count_mismatches(s1, s2)}.")
print()

# 3)

# 4)

