from ordered_set import OrderedSet


class AlignedSequences:
    def __init__(self, sequences):
        self.sequences = sequences

    def consensus(self) -> str:
        """Calculate the consensus of all sequences with the new sequence"""

        result = ""
        transposed_sequences = zip(*self.sequences)

        for column in transposed_sequences:
            frequencies = tuple((item, column.count(item)) for item in OrderedSet(column) if item != "-")
            highest_frequency = max(frequencies, key=lambda x: x[1])[1]
            most_frequent = tuple(x[0] for x in frequencies if x[1] == highest_frequency)

            result += most_frequent[0]

        return result

    def __len__(self) -> int:
        return len(self.sequences)

    def __getitem__(self, n):
        if type(n) is tuple and len(n) == 2:
            i, j = n

            return self.sequences[i][j]

        elif type(n) is int:
            return self.sequences[n]

        return None
