from rich import print
from read_fasta import read_fasta


def search_first_occ(seq, pattern):
    found = False
    i = 0
    while i <= len(seq) - len(pattern) and not found:
        j = 0
        while j < len(pattern) and pattern[j] == seq[i + j]:
            j = j + 1
        if j == len(pattern):
            found = True
        else:
            i += 1
    if found:
        return i
    else:
        return -1


def search_all_occurrences(seq, pattern):
    ''' returns a list of the indices where the pattern occurrs within seq'''
    res = []

    for i in range(len(seq) - len(pattern)):
        if seq[i:i + len(pattern)] == pattern:
            res.append(i)

    return res


def search_all_occurrences_mismatch(seq, pattern, m):
    '''Find all occurrences of pattern in seq with maximum hamming distance of m'''
    res = []

    for i in range(len(seq) - len(pattern)):
        if hamming_distance(seq[i:i + len(pattern)], pattern) <= m:
            res.append(i)

    return res


def hamming_distance(s, p):
    distance = 0

    if len(s) != len(p):
        return -1

    else:
        for i in range(0, len(s)):
            if s[i] == p[i]:
                distance += 1

    return distance


def simpleBooyerMoore(seq, pattern):
    """Very simplified version of Booyer-Moore with the BCR rule"""
    alphabet = "".join(set(list(seq)))
    # process Bad-character rule
    occ = {}
    for symb in alphabet:
        occ[symb] = -1
    for j in range(len(pattern)):
        c = pattern[j]
        occ[c] = j
    # search pattern
    res = []
    i = 0
    while i <= len(seq) - len(pattern):
        j = len(pattern) - 1
        print(i)
        while j >= 0 and pattern[j] == seq[j + i]:
            j -= 1
        if (j < 0):
            res.append(i)
            i += 1
        else:
            c = seq[j + i]
            i += max(1, j - occ[c])
    return res


def repeated_subsequences_frequency(dna_seq, k=10):
    '''Write a function that, given a DNA sequence, allows to detect if there are repeated sequences of size k
    The result should be a dictionary with sub-sequences as keys, and their frequency as values.'''

    subsequences = [dna_seq[i:i + k] for i in range(len(dna_seq) - k + 1)]

    return {x: subsequences.count(x) for x in subsequences if subsequences.count(x) > 1}


def finds_patterns_frequency(seqA, seqB, low, high):
    patterns = set()

    for k in range(low, high + 1):
        patterns |= set(seqA[i:i + k] for i in range(len(seqA) - k + 1))

    occurrences = {pattern: count for pattern in patterns if (count := len(search_all_occurrences(seqB, pattern))) > 0}

    return dict(sorted(occurrences.items(), key=lambda x: x[1], reverse=True))


def intron_boundaries(dna_seq):
    sequences = {}

    for i in search_all_occurrences(dna_seq, "GTG"):
        for j in search_all_occurrences(dna_seq, "CAG"):
            # Add intron sequence if a CAG occurs after a GTG signal
            # If a sequence already exists for a given CAG signal, update it if the new one is shorter
            if j > i and (j not in sequences or sequences[j] < i):
                sequences[j] = i

                # Stop at first occurrence of CAG signal
                break

    return [(item[1], item[0] + 3) for item in sequences.items()]


def validate_dna_re(seq):
    from re import search

    if search("[^ACTGactg]", seq) != None:
        return False

    else:
        return True


def translate_codon_re(cod):
    import re
    if re.search("GC.", cod):
        aa = "A"
    elif re.search("TG[TC]", cod):
        aa = "C"
    elif re.search("GA[TC]", cod):
        aa = "D"
    elif re.search("GA[AG]", cod):
        aa = "E"
    elif re.search("TT[TC]", cod):
        aa = "F"
    elif re.search("GG.", cod):
        aa = "G"
    elif re.search("CA[TC]", cod):
        aa = "H"
    elif re.search("AT[TCA]", cod):
        aa = "I"
    elif re.search("AA[AG]", cod):
        aa = "K"
    elif re.search("TT[AG]|CT.", cod):
        aa = "L"
    elif re.search("ATG", cod):
        aa = "M"
    elif re.search("AA[TC]", cod):
        aa = "N"
    elif re.search("CC.", cod):
        aa = "P"
    elif re.search("CA[AG]", cod):
        aa = "Q"
    elif re.search("CG.|AG[AG]", cod):
        aa = "R"
    elif re.search("TC.|AG[TC]", cod):
        aa = "S"
    elif re.search("AC.", cod):
        aa = "T"
    elif re.search("GT.", cod):
        aa = "V"
    elif re.search("TGG", cod):
        aa = "W"
    elif re.search("TA[TC]", cod):
        aa = "Y"
    elif re.search("TA[AG]|TGA", cod):
        aa = "_"
    else:
        aa = ""
    return aa


def find_pattern_re(seq, pat):
    from re import search
    mo = search(pat, seq)
    if (mo != None):
        return mo.span()[0]
    else:
        return -1


def find_all_occurrences_re(seq, pat):
    from re import finditer
    mos = finditer(pat, seq)
    res = []
    for x in mos:
        res.append(x.span()[0])
    return res


def find_all_overlap(seq, pat):
    return find_all_occurrences_re(seq, "(?=" + pat + ")")


def test_RE():
    seq = input("Input sequence:")
    pat = input("Input pattern (as a regular expression):")

    res = find_pattern_re(seq, pat)
    if res >= 0:
        print("Pattern found in position: ", res)
    else:
        print("Pattern not found")

    all_res = find_all_occurrences_re(seq, pat)
    if len(all_res) > 0:
        print("Pattern found in positions: ", all_res)
    else:
        print("Pattern not found")

    all_ov = find_all_overlap(seq, pat)
    if len(all_ov) > 0:
        print("Pattern found in positions (overlap): ", all_ov)
    else:
        print("Pattern not found")


def test():
    fasta = read_fasta("test_files/HBA1.DNA.fasta")
    sequence = fasta["HBA1"]

    print(f"All occurrences: {search_all_occurrences(sequence, 'AAT')}")
    print(f"Patterns frequency: {finds_patterns_frequency('AAATG', 'CAAATC', 2, 3)}")

    boundaries = intron_boundaries(sequence.upper())
    introns = [sequence[i:j] for (i, j) in boundaries if len(sequence[i:j]) % 3 == 0]

    print("Introns:")
    print(introns)


if __name__ == "__main__":
    test()
