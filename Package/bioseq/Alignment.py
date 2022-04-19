class Alignment:
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2

    def read_submat_file(filename):
        """Read substitution matrix from file """
        
        self.sm = {}
        f = open(filename, "r")
        line = f.readline()
        tokens = line.split("\t")
        ns = len(tokens)
        alphabet = []

        for i in range(0, ns):
            alphabet.append(tokens[i][0])
        
        for i in range(0, ns):
            line = f.readline()
            tokens = line.split("\t")
            
            for j in range(0, len(tokens)):
                k = alphabet[i] + alphabet[j]
                self.sm[k] = int(tokens[j])
        
    # Global alignment
    def needleman_wunsch():
        pass

    # Local alignment
    def smith_waterman():
        pass
