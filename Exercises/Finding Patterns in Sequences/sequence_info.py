import sys

from read_fasta import read_fasta

def sequence_identifiers(fasta):
    return fasta.keys()

if __name__ == "__main__":
    filename = sys.argv[1]
    fasta = read_fasta(f"test_files/{filename}")

    print("Sequence identifiers:")
    
    for identifier in sequence_identifiers(fasta):
        print(f"\t{identifier}") 

    print()
    
     