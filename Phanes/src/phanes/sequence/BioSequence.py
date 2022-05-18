from __future__ import annotations

class BioSequence:
    def __init__(self, sequence: str):
        if type(self) == BioSequence:
            raise Exception("BioSequence must be subclassed.")

        elif not self.validate(sequence):
            raise Exception("Sequence is not valid")

        self.sequence = sequence

    @staticmethod
    def read_fasta(filename):
        sequences = {}
        current_sequence = ""

        with open(filename) as fasta:
            for line in fasta.readlines():
                if line[0] == ">":
                    if current_sequence != "":
                        sequences[identifier] = current_sequence
                    identifier = line[1:].strip()
                    current_sequence = ""

                else:
                    current_sequence += line.strip()

        sequences[identifier] = current_sequence

        return sequences

    def __len__(self) -> int:
        return len(self.sequence)

    def __str__(self) -> str:
        return self.sequence

    def __eq__(self, other):
        return self.sequence == other.sequence
