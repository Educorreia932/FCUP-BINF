import sys
from read_fasta import read_fasta
from prosite_re import *

if __name__ == "__main__":
    filename = sys.argv[1]
    fasta = read_fasta(f"test_files/{filename}")
    motif = "N-x-G-x-R-[LIVM]-D-[LIVMFYH]-x-[LV]-x-S"

    for identifier, sequence in fasta.items():
        matches = find_prosite(sequence, motif) != -1

        print(f"{identifier} {'MATCH' if matches else 'NOT_MATCH'}")

