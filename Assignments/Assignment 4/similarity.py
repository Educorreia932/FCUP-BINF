from bioseq.BioSequence import BioSequence
from bioseq.Alignment import NeedlemanWunsch

proteins = list(sequence for sequence in BioSequence.read_fasta("glycoproteinS.fas").values())

seq1 = proteins[0] # YP_009724390.1
seq2 = proteins[1] # QHO60594.1

# 1)
needleman_wunsch = NeedlemanWunsch(seq1, seq2, -8)
needleman_wunsch.calculate()
s1, s2 = needleman_wunsch.recover_align()
mismatches = 0

for c1, c2 in zip(s1, s2):
    if c1 != c2:
        mismatches += 1

print(f"The number of mismatches in global alignment is {mismatches}")

# 2)

# 3)

# 4)

