from .BioSequence import BioSequence

import re

class ORF(BioSequence):
    def __init__(self, sequence: str):
        BioSequence.__init__(self, sequence)

    @staticmethod
    def validate(sequence):
        return bool(re.match(r"^[A-IK-NP-Z\*-]+$", sequence))