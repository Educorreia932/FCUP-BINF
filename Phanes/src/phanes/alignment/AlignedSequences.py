class AlignedSequences:
    def __init__(self, sequences):
        self.sequences = sequences

    def consensus(self) -> str:
        result = ""
        transposed_sequences = zip(*self.sequences)

        for column in transposed_sequences:
            frequencies = tuple((item, column.count(item)) for item in set(column))
            highest_frequency = max(frequencies, key=lambda x: x[1])[1]
            most_frequent = sorted(x[0] for x in frequencies if x[1] == highest_frequency)

            result += most_frequent[0]

        return result
