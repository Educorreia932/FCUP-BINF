class Blast:
    def __init__(self, filename: str | None = None, word_length: int = 3):
        if filename is None:
            self.database = []

        else:
            self.database = self._read_database(filename)

        self.word_length = word_length

    @staticmethod
    def _read_database(filename: str) -> [str]:
        """From file with sequences line by line read the sequences to a list"""

        with open(filename, "r") as f:
            return [line for line in f.readlines()]

    def _build_map(self, query):
        result = {}

        for i in range(len(query) - self.word_length + 1):
            subsequence = query[i:i + self.word_length]

            if subsequence not in in result:
                result[subsequence] = [i]

            else:
                result[subsequence].append(i)

        return result

    def get_hits(self, sequence: str, query: str):
        result = []
        self.map = self._build_map(query)

        for i in range(len(sequence) - self.word_length + 1):
            subsequence = sequence[i:i + self.word_length]

            if subsequence in self.map:
                l = self.map[subsequence]

                for j in l:
                    result.append((j, i))

        return result
