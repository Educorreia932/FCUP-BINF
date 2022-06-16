def dna_ohe(sequence):
    encoding = {x:"".join("1" if i == j else "0" for j in range(4)) for i, x in enumerate("ACGT")}

    return "".join(encoding[x] for x in sequence)

def word_to_kmer(word, k):
    subsequences = set(word[i:i + k] for i in range(len(word) - k + 1))
    frequencies = {subsequence:word.count(subsequence) for subsequence in subsequences}

    return frequencies

def file_to_kmer_table(file_name):
    pass


if __name__ == "__main__":
    sequence = "GTAGAGCTGT"

    print(dna_ohe(sequence))
    print(word_to_kmer(sequence, 2))
