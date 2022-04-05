def read_fasta(filename):
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
    