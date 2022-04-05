import re

def pubmed_identifiers(gb):
    result = []
    lines = [line.strip() for line in gb.split("\n")]

    for line in lines:
        if line.find("PUBMED") != -1:
            result.append(line)

    return result

def list_genes(gb):
    genes = set()

    for gene in re.findall(r"CDS.*\n.*\/gene=\"(.+)\"", gb):
        genes.add(gene)
    
    return genes

def list_proteins(gb):
    proteins = []

    for protein in re.findall(r"CDS.*\n(?:.*\n)*?.*\/protein_id=\"(.+)\"", gb):
        proteins.append(protein)

    return proteins


def translated_protein_sequences(gb):
    CDS\s*.*\n\s*/gene=\"(\w+)\"\n(?:.*\n)*?\s*?/translation=\"(\w+\n\s+)
    

if __name__ == "__main__":
    with open("test_files/sequence.gb") as gb:
        gb = gb.read()

        # a)
        print("PUBMED Identifiers:")

        for identifier in pubmed_identifiers(gb):
            print(f"\t{identifier}")

        print()

        # b)
        print("List of genes (ORFs):")

        for gene in list_genes(gb):
            print(f"\t{gene}")

        print()

        # c)
        print("Protein identifiers:")

        for protein in list_proteins(gb):
            print(f"\t{protein}")
        