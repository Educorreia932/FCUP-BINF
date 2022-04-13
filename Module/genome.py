class BioSequence():
    def __init__(self, sequence: str):
        if type(self) == BioSequence:
            raise Exception("BioSequence must be subclassed.")

        self.sequence = sequence

    @staticmethod
    def translate_codon(codon: str):
        translation = {
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

        if codon in translation:
            return translation[codon]

        else:
            return None

    def translate(k=0):
        aminoacid_sequence = "".join(self.translate_codon(self.sequence[i:i + 3]) for i in range(k, len(self.sequence, len(self.sequence - k + 2), 3)))

        return Protein(aminoacid_sequence)

    def all_proteins():
        pass

    def __len__(self) -> int:
        return len(self.sequence)

    def __str__(self) -> int:
        return self.sequence

class DNA(BioSequence):
    def __init__(self, sequence):
        BioSequence.__init__(self, sequence)

    def reverse_complement(self) -> str:
        complement = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C"
        }

        return "".join(complement[nucleotide] for nucleotide in self.sequence[::-1])

    def transcribe(self):
        return RNA(self.sequence.replace("T", "U"))


class RNA(BioSequence):
    def __init__(self, sequence: str):
        BioSequence.__init__(self, sequence)

    def transcribe(self):
        return DNA(self.sequence.replace("U", "T"))


class Protein():
    def __init__(self, sequence: str):
        self.sequence = sequence

# TODO: Read from file
# TODO: Validate sequence