from .AlignedSequences import AlignedSequences
from .PairwiseAlignment import NeedlemanWunsch


class MultipleAlignment:
    def __init__(self, sequences: list[str]):
        self.sequences = sequences
        self.alignment = NeedlemanWunsch(sm=NeedlemanWunsch.create_submatrix(1, -1, "ACGT"))

    def add_alignment(self, aligned: AlignedSequences, sequence) -> AlignedSequences:
        result = ["" for _ in range(len(aligned) + 1)]

        # Create consensus from give alignments
        consensus = aligned.consensus()
        self.alignment.set_sequences(consensus, sequence)
        self.alignment.calculate()

        aligned2 = self.alignment.recover_align()
        origin = 0

        for i in range(len(aligned2[0])):
            # Check if is a gap
            if aligned2[0, i] == "-":
                for k in range(len(aligned)):
                    result[k] += "-"

            else:
                for k in range(len(aligned.sequences)):
                    result[k] += aligned[k, origin]

                origin += 1

        result[len(aligned.sequences)] = aligned2.sequences[1]

        return AlignedSequences(result)

    def calculate(self) -> AlignedSequences:
        self.alignment.set_sequences(self.sequences[0], self.sequences[1])
        self.alignment.calculate()

        result = self.alignment.recover_align()

        for sequence in self.sequences[2:]:
            result = self.add_alignment(result, sequence)

        return result
