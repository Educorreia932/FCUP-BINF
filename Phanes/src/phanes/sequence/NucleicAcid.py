from __future__ import annotations

from .ORF import ORF
from .BioSequence import BioSequence

import re


class NucleicAcid(BioSequence):
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

    def translate(self, k) -> str:
        a = [self.translate_codon(self.sequence[i:i + 3]) for i in range(k, len(self.sequence) - 2, 3)]

        return "".join(a)

    def compute_ORFs(self, k) -> [ORF]:
        result = []
        aminoacid_sequence = self.translate(k)
        current_orf = ""

        for aminoacid in aminoacid_sequence:
            # Stop codon
            if aminoacid == "_" and current_orf != "":
                result.append(ORF(current_orf))

                current_orf = ""

            # Start codon or ongoing sequence
            elif aminoacid == "M" or current_orf != "":
                current_orf += aminoacid

        return result

    def all_ORFs(self) -> [ORF]:
        result = []

        # Frames starting in positions 0, 1 and 2
        for i in range(3):
            result += self.compute_ORFs(i)

        # Compute in both ways
        if type(self) == DNA:
            reverse_complement = self.reverse_complement()

            for i in range(3):
                result += reverse_complement.compute_ORFs(i)

        return result


class DNA(NucleicAcid):
    def __init__(self, sequence):
        BioSequence.__init__(self, sequence)

    @staticmethod
    def validate(sequence: str) -> bool:
        return bool(re.match(r"^[ACGT]+$", sequence))

    def reverse_complement(self) -> DNA:
        complement = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C"
        }

        return DNA("".join(complement[nucleotide] for nucleotide in self.sequence[::-1]))

    def transcribe(self) -> RNA:
        return RNA(self.sequence.replace("T", "U"))


class RNA(NucleicAcid):
    def __init__(self, sequence: str):
        BioSequence.__init__(self, sequence)

    @staticmethod
    def validate(sequence: str) -> bool:
        return bool(re.match(r"^[ACGU]+$", sequence))

    def transcribe(self) -> DNA:
        return DNA(self.sequence.replace("U", "T"))
