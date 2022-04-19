from bioseq.BioSequence import BioSequence
from bioseq.Alignment import Alignment

proteins = BioSequence.read_fasta("glycoproteinS.fas")

seq1 = proteins[0]
seq2 = proteins[1]

alignment = Alignment(seq1, seq2)

# 1)
alignment.needleman_wunsch()

# 2)
alignment.smith_waterman()

# 3)

# 4)

