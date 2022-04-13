from genome import *

if __name__ == "__main__":
    sequence = "ATGAAATTATGAATGAGCCTCAGCTGAAGCATCGCGCATCAGACTACGCTCAGACTCAGACTCAGCATTATAGTGAATGTTAATAAATAAAATAA"

    dna = DNA(sequence)

    print(f"Sequence: {dna}")
    print(f"Length: {len(dna)}")
    print(f"RNA: {dna.transcribe()}")
    print(f"Reverse complement: {dna.reverse_complement()}")
    print("ORFs:")

    for orf in dna.all_ORFs():
        print(f"\t{orf}")
