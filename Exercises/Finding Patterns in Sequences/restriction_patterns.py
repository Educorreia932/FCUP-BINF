from read_fasta import read_fasta
from re import finditer


def iub_to_RE(iub):
    """converts IUB ambiguity code into regular expression"""
    dic = {"A": "A", "C": "C", "G": "G", "T": "T",
           "R": "[GA]", "Y": "[CT]", "M": "[AC]", "K": "[GT]", "S": "[GC]", "W": "[AT]",
           "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}

    site = iub.replace("^", "")
    regexp = ""

    for c in site:
        regexp += dic[c]

    return regexp


def cut_positions(enzyme, sequence):
    """Returns cut position of a restriction enzyme (in IUB code) in a sequence"""

    cutpos = enzyme.find("^")
    regexp = iub_to_RE(enzyme)

    matches = finditer(regexp, sequence)
    locs = []

    for m in matches:
        locs.append(m.start() + cutpos)

    return locs

# determines subsequences resulting from a sequence cut in a list of positions


def cut_subsequences(locs, sequence):
    """Determines subsequences resulting from a sequence cut in a list of positions
"""
    res = []
    positions = locs
    positions.insert(0, 0)
    positions.append(len(sequence))

    for i in range(len(positions) - 1):
        res.append(sequence[positions[i]:positions[i + 1]])

    return res


def test():
    print(iub_to_RE("G^AMTV"))
    pos = cut_positions("G^ATTC", "GTAGAAGATTCTGAGATCGATTC")
    print(pos)
    print(cut_subsequences(pos, "GTAGAAGATTCTGAGATCGATTC"))


