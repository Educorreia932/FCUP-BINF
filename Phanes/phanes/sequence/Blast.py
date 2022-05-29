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

            if subsequence not in result:
                result[subsequence] = [i]

            else:
                result[subsequence].append(i)

        return result

    def get_hits(self, sequence: str, query: str):
        result = []
        hashmap = self._build_map(query)

        for i in range(len(sequence) - self.word_length + 1):
            subsequence = sequence[i:i + self.word_length]

            if subsequence in hashmap:
                for j in hashmap[subsequence]:
                    result.append((j, i))

        return result

    def extend_hit(self, sequence, hit, query):
        stq, sts = hit[0], hit[1]
        matfw = 0
        k = 0
        bestk = 0

        while 2 * matfw >= k and stq + self.word_length + k < len(query) and sts + self.word_length + k < len(sequence):
            if query[stq + self.word_length + k] == sequence[sts + self.word_length + k]:
                matfw += 1
                bestk = k + 1

            k += 1

        size = self.word_length + bestk

        k = 0
        matbw = 0
        bestk = 0

        while 2 * matbw >= k and stq > k and sts > k:
            if query[stq - k - 1] == sequence[sts - k - 1]:
                matbw += 1
                bestk = k + 1
            k += 1

        size += bestk

        return stq - bestk, sts - bestk, size, self.word_length + matfw + matbw

    def hit_best_score(self, sequence, query):
        hits = self.get_hits(sequence, query)
        best_score = -1
        best = ()

        for hit in hits:
            extension = self.extend_hit(sequence, hit, query)
            score = extension[3]

            if score > best_score or (score == best_score and extension[2] < best[2]):
                best_score = score
                best = extension

        return best

    def calculate(self, query: str):
        """Calculates best alignment"""

        best_score = -1
        result = (0, 0, 0, 0, 0)

        for k, sequence in enumerate(self.database):
            best_hit = self.hit_best_score(sequence, query)

            if len(best_hit) != 0:
                score = best_hit[3]

                if score > best_score or (score == best_score and best_hit[2] < result[2]):
                    best_score = score
                    result = (best_hit[0], best_hit[1], best_hit[2], best_hit[3], k)

        if best_score < 0:
            result = tuple()

        return result
