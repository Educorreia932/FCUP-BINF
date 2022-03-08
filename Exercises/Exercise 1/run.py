def count_occurrences(s):
    return " ".join(f"{item[0]}:{item[1]}" for item in sorted({x:s.count(x) for x in set(s)}.items()))

print(count_occurrences("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"))
