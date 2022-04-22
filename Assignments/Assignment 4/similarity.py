from bioseq.BioSequence import BioSequence
from bioseq.Alignment import Alignment

proteins = list(sequence for sequence in BioSequence.read_fasta("glycoproteinS.fas").values())

seq1 = proteins[0] # YP_009724390.1
seq2 = proteins[1] # QHO60594.1

alignment = Alignment(seq1, seq2, -8)

# 1)
S1, T1 = alignment.needleman_wunsch()
mismatches = 0
s1, s2 = alignment.recover_align(T1)

for c1, c2 in zip(s1, s2):
    if c1 != c2:
        mismatches += 1

print(f"The number of mismatches in global alignment is {mismatches}")

# 2)
# alignment.smith_waterman()

# 3)

# 4)

