class AlignedSequences:
    def __init__(self, sequences):
        self.sequences = sequences

    def consensus(self) -> list[str]:
        result = []
        transposed_sequences = zip(*self.sequences)

        for column in transposed_sequences:
            most_common = max(((item, column.count(item)) for item in set(column)), key=lambda a: a[1])[0]

            result.append(most_common)

        return result
