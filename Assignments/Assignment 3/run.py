def frequency(seq, percentage=False):
    """ Calculates the frequency of each symbol in the sequence. Returns a dictionary. """

    if percentage:
        return {x: seq.count(x) / len(seq) for x in set(seq)}

    else:
        return {x: seq.count(x) for x in set(seq)}


def gc_content(dna_seq):
    """ Returns the percentage of G and C nucleotides in a DNA sequence. """
    gc_count = 0

    for s in dna_seq:
        if s in "GCgc":
            gc_count += 1

    return gc_count / len(dna_seq)


def reverse_complement(dna_seq):
    """ Computes the reverse complement of the inputted DNA sequence. """

    return dna_seq[::-1] \
        .replace("A", "t") \
        .replace("C", "g") \
        .replace("G", "C") \
        .replace("T", "A") \
        .upper()


def translate_codon(cod):
    """Translates a codon into an aminoacid using an internal dictionary with the standard genetic code."""

    tc = {
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TGT": "C", "TGC": "C",
        "GAT": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",
        "TTT": "F", "TTC": "F",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "CAT": "H", "CAC": "H",
        "ATA": "I", "ATT": "I", "ATC": "I",
        "AAA": "K", "AAG": "K",
        "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATG": "M", "AAT": "N", "AAC": "N",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TGG": "W",
        "TAT": "Y", "TAC": "Y",
        "TAA": "_", "TAG": "_", "TGA": "_"
    }

    if cod in tc:
        return tc[cod]

    else:
        return None


def translate_seq(dna_seq, ini_pos=0):
    """Translates a DNA sequence into an aminoacid sequence."""

    seqm = dna_seq.upper()

    return "".join(translate_codon(seqm[i:i + 3]) for i in range(ini_pos, len(seqm) - 2, 3))


def reading_frames(dna_seq):
    """Computes the six reading frames of a DNA sequence (including the reverse complement."""

    res = []

    for i in range(3):
        res.append(translate_seq(dna_seq, i))

    rc = reverse_complement(dna_seq)

    for i in range(3):
        res.append(translate_seq(rc, i))

    return res


def all_proteins_rf(aa_seq, coordinates=False):
    """Computes all possible proteins in an aminoacid sequence."""

    aa_seq = aa_seq.upper()
    current_prot = ""
    proteins = []

    for i, aa in enumerate(aa_seq):
        if aa == "_" and current_prot != "":
            if coordinates:
                proteins.append((i - len(current_prot), i, current_prot))

            else:
                proteins.append(current_prot)

            current_prot = ""

        else:
            if aa == "M" or current_prot != "":
                current_prot += aa

    return proteins


def all_orfs(dna_seq, minsize=1, sorted=True, coordinates=False, ):
    """Computes all possible proteins for all open reading frames."""
    res = []

    # If coordinates are enabled, the ORF is the last element of the tuple (Start1,	End1, ORF1). Otherwise, it is just a string
    def orf(x): return x[2] if coordinates else x

    for frame in reading_frames(dna_seq):
        res += all_proteins_rf(frame, coordinates=coordinates)

    # Sort ORFs
    if sorted:
        res.sort(key=lambda x: len(orf(x)))

    # Filter ORFs by minimum size
    res = filter(lambda x: len(orf(x)) >= minsize, res)

    return res


def read_fasta_2dictionary(filename):
    """Reads a FASTA file to a dictionary"""

    result = {}
    sequence = ""

    with open(filename) as fasta:
        for line in fasta.readlines():
            if line[0] == ">":
                if sequence != "":
                    result[identifier] = sequence

                identifier = line.split()[0][1:]

                sequence = ""

            else:
                sequence += line.strip()

        result[identifier] = sequence

    return result


if __name__ == "__main__":
    dna_seq = read_fasta_2dictionary("sequence.fasta")["NC_045512.2"]

    # 1. Length of the sequence
    print(len(dna_seq))

    # 2. Frequency (in %) of A, C, G, T
    print(frequency(dna_seq, percentage=True))

    # 3. GC content
    print(gc_content(dna_seq))

    # 4. Number of Start (AUG) codons found
    aa_seq = translate_seq(dna_seq)

    print(aa_seq.count("M"))

    # 5. Number of Stop Codons (UAA, UAG, UGA)
    print(aa_seq.count("_"))

    # 6. Most and less frequent codons
    print(max(frequency(aa_seq).items(), key=lambda x: x[1]))
    print(min(frequency(aa_seq).items(), key=lambda x: x[1]))

    # 7. Output protein sequences to file
    orfs = all_orfs(dna_seq, minsize=40)

    with open("all_potential_proteins.txt", "w") as outfile:
        outfile.writelines(orf + "\n" for orf in orfs)

    # 8. Output genomic	coordinates to file
    orfs_coordinates = all_orfs(dna_seq, minsize=40, coordinates=True)

    with open("orf_coordinates.txt", "w") as outfile:
        outfile.writelines(f"{', '.join(str(x) for x in orf)}\n" for orf in orfs_coordinates)
