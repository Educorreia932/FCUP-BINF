from re import search
from read_fasta import read_fasta


def find_zync_finger(seq):
    from re import search
    regexp = "C.H.[LIVMFY]C.{2}C[LIVMYA]"
    mo = search(regexp, seq)
    if (mo != None):
        return mo.span()[0]
    else:
        return -1


def find_prosite(seq, profile):
    regexp = profile.replace("-", "")
    regexp = regexp.replace("x", ".")
    regexp = regexp.replace("(", "{")
    regexp = regexp.replace(")", "}")
    mo = search(regexp, seq)

    if (mo != None):
        return mo.span()

    else:
        return -1


def test():
    seq = "HKMMLASCKHLLCLKCIVKLG"
    print(find_zync_finger(seq))
    print(find_prosite(seq, "C-x-H-x-[LIVMFY]-C-x(2)-C-[LIVMYA]"))


# Task 2
def taks_2():
    prosite = "C-x-H-x-[LIVMFY]-C-x(2)-C-[LIVMYA]"
    sequence = read_fasta("test_files/Q8RXD4.fasta")["sp|Q8RXD4|BRCA1_ARATH"]
    i, j = find_prosite(sequence, prosite)

    print(sequence[:i].lower() + sequence[i:j].upper() + sequence[j:].lower())


if __name__ == "__main__":
    taks_2()
