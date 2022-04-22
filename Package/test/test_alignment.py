from bioseq.Alignment import Alignment
from rich import print

seq1 = "PHSWG"
seq2 = "HGWAG"

alignment = Alignment(seq1, seq2, -8)

S1, T1 = alignment.needleman_wunsch()
print(alignment.recover_align(T1))

# print(alignment.smith_waterman())