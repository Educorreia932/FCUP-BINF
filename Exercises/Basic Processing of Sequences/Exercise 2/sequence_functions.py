from rich import print

def read_seq_from_file(filename):
    """ Reads a sequence from a multi-line text file. """

    fh = open(filename, "r")
    lines = fh.readlines()
    seq = "".join(line.replace("\n", "") for line in lines)
    
    fh.close()
    
    return seq

def write_seq_to_file(seq, filename):
    """ Writes a sequence to file. """

    with open(filename, "w") as outfile:
        outfile.writelines([seq[i:i + 60] + "\n" for i in range(0, len(seq), 60)])

def read_genetic_code_from_file(filename):
    """ Reads the genetic code to a dictionary from a multi-line text file. """

    with open(filename, "r") as f:
        return {line[0]:line[1] for line in map(lambda x: x.replace("\"", "").split(), f.readlines())} 

def validate_dna (dna_seq):
    """ Checks if DNA sequence is valid. Returns True is sequence is valid, or False otherwise. """

    seqm = dna_seq.upper()

    return all(x in {"A", "C", "G", "T"} for x in seqm)

def frequency (seq):
    """ Calculates the frequency of each symbol in the sequence. Returns a dictionary. """

    return {x:seq.count(x) for x in set(seq)}

def gc_content (dna_seq):
    """ Returns the percentage of G and C nucleotides in a DNA sequence. """
    gc_count = 0
    for s in dna_seq:
        if s in "GCgc": gc_count += 1
    return gc_count / len(dna_seq)

def gc_content_subseq (dna_seq, k=100):
    """ Returns GC content of non-overlapping sub-sequences of size k. """
    # complete
    # ...

def transcription (dna_seq):
    """ Function that computes the RNA corresponding to the transcription of the DNA sequence provided. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"

    return dna_seq.upper().replace("T", "U")

def reverse_complement (dna_seq):
    """ Computes the reverse complement of the inputted DNA sequence. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
 
    return dna_seq[::-1] \
        .replace("A", "t") \
        .replace("C", "g") \
        .replace("G", "C") \
        .replace("T", "A") \
        .upper()

def translate_codon (cod):
    """Translates a codon into an aminoacid using an internal dictionary with the standard genetic code."""

    tc = {
        "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
        "TGT":"C", "TGC":"C",
        "GAT":"D", "GAC":"D",
        "GAA":"E", "GAG":"E",
        "TTT":"F", "TTC":"F",
        "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
        "CAT":"H", "CAC":"H",
        "ATA":"I", "ATT":"I", "ATC":"I",
        "AAA":"K", "AAG":"K",
        "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
        "ATG":"M", "AAT":"N", "AAC":"N",
        "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
        "CAA":"Q", "CAG":"Q",
        "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
        "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
        "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
        "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
        "TGG":"W",
        "TAT":"Y", "TAC":"Y",
        "TAA":"_", "TAG":"_", "TGA":"_"
    }
    
    if cod in tc: 
        return tc[cod]
    
    else: 
        return None

def translate_seq (dna_seq, ini_pos = 0):
    """Translates a DNA sequence into an aminoacid sequence."""

    assert validate_dna(dna_seq), "Invalid DNA sequence"

    seqm = dna_seq.upper()

    return "".join(translate_codon(seqm[i:i + 3]) for i in range(ini_pos, len(seqm) - 2, 3))

def codon_usage(dna_seq, aa):
    """Provides the frequency of each codon encoding a given aminoacid, in a DNA sequence ."""

    assert validate_dna(dna_seq), "Invalid DNA sequence"
    
    seqm = dna_seq.upper()

    aa_codons = [seqm[i:i + 3] for i in range(0, len(seqm) - 2, 3) if codon == aa]
    
    return {codon:aa_codons.count() / len(aa_codons) for codon in set(aa_codons)}

def reading_frames (dna_seq):
    """Computes the six reading frames of a DNA sequence (including the reverse complement."""

    assert validate_dna(dna_seq), "Invalid DNA sequence"

    res = []

    for i in range(3):
        res.append(translate_seq(dna_seq, i))

    rc = reverse_complement(dna_seq)
    
    for i in range(3):
        res.append(translate_seq(rc, i))
    
    return res

def all_proteins_rf(aa_seq):
    """Computes all posible proteins in an aminoacid sequence."""
    aa_seq = aa_seq.upper()
    current_prot = []
    proteins = []

    for aa in aa_seq:
        if aa == "_" and current_prot:
            for p in current_prot:
                proteins.append(p)

            current_prot = []
        
        else:
            if aa == "M":
                current_prot.append("")
            
            for i in range(len(current_prot)):
                current_prot[i] += aa

    return proteins

def all_orfs(dna_seq):
    """Computes all possible proteins for all open reading frames."""
    res = []

    for frame in reading_frames(dna_seq):
        res += all_proteins_rf(frame)

    return res

def all_orfs_ord(dna_seq, minsize = 0):
    """Computes all possible proteins for all open reading frames. Returns ordered list of proteins with minimum size."""
    
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    
    return filter(lambda x: len(x) > minsize, sorted(all_orfs(dna_seq), key=lambda x: len(x)))

def insert_prot_ord (prot, list_prots):
    ''' inserts prot in list_prots in a sorted way '''

    i = 0

    while i < len(list_prots) and len(prot) < len(list_prots[i]):
        i += 1
    
    list_prots.insert(i, prot)

def test_frequency():
    seq_aa = input("Protein sequence:")
    freq_aa = frequency(seq_aa)
    list_f = sorted(freq_aa.items(), key=lambda x: x[1], reverse = True)

    for (k,v) in list_f:
        print("Aminoacid:", k, ":", v)

def test_all(seq):
    if validate_dna (seq):
        assert validate_dna(seq), "Invalid DNA sequence"

        print("Valid sequence\n")
        print(f"[bold]Transcription:[/bold] {transcription(seq)}\n")
        print(f"[bold]Reverse complement:[/bold] {reverse_complement(seq)}\n")
        print(f"[bold]GC Content:[/bold] {gc_content(seq) * 100}%\n")
        print(f"[bold]Direct translation:[/bold] {translate_seq(seq)}")

        orfs = all_orfs_ord(seq)

        for i, orf in enumerate(orfs):
            write_seq_to_file(orf, f"orf/orf-{i + 1}.txt")
    
    else: 
        print("DNA sequence is not valid")

def test_files():
    fname = input("Insert input filename:")
    test_sequence_from_file(fname)

def test_sequence_from_file(filename):
    seq = read_seq_from_file(filename)

    test_all(seq)

if __name__ == "__main__":
    write_seq_to_file("ATGAGCGACAT" * 10, "seq.txt")
    test_sequence_from_file("example_Hinfluenzae.txt")
