from .AlignedSequences import AlignedSequences
from .PairwiseAlignment import NeedlemanWunsch


class MultipleAlignment:
    def __init__(self, sequences: list[str], gap):
        self.sequences = sequences
        self.gap = gap

    def add_alignment(self, aligned: AlignedSequences, sequence) -> AlignedSequences:
        result = ["" for _ in range(len(aligned) + 1)]

        # Create consensus from give alignments
        consensus = aligned.consensus()
        alignment = NeedlemanWunsch(consensus, sequence, self.gap)
        aligned2 = alignment.recover_align()
        origin = 0

        for i in range(len(aligned2[0])):
            # Check if is a gap
            if aligned2[0, i] == "_":
                for k in range(len(aligned.sequences)):
                    result[k] += "_"

            else:
                for k in range(len(aligned.sequences)):
                    result[k] += aligned[k, origin]

                origin += 1

        result[len(aligned.sequences)] = aligned2.sequences[1]

        return AlignedSequences(result)

    def align(self) -> AlignedSequences:
        alignment = NeedlemanWunsch(self.sequences[0], self.sequences[1], self.gap)
        result = alignment.recover_align()

        for sequence in self.sequences[2:]:
            result = self.add_alignment(result, sequence)

        return result
