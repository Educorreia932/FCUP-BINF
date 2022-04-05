"""
Class with methods to read, transcribe, translate and get features from DNA or RNA genomes

Bruno Vaz - up201705247
Eduardo Correia - up201806433
Filipe Justiça - up201606339 

"""
from typing import Union
from collections import Counter

class genome:
    def __init__(self):
        """
        Genome sequence class

        :param sequence: string composed by sequence of DNA or RNA
        """

        self.isDNA = False
        self.isRNA = False
        self.seq = None

    def read_from_str(self, sequence: str) -> None:
        """
        Read sequence from string, asserts if DNA or RNA
        
        :param sequence: string composed by sequence of DNA or RNA
        """
        if (not isinstance(sequence, str)): raise TypeError("Sequence must be a string")
        if len(sequence) < 3: raise ValueError("Sequence must have at least one codon")

        self.isDNA = self.validate_dna(sequence)
        self.isRNA = False

        if not self.isDNA:
            self.isRNA = self.validate_rna(sequence)

        if not self.isDNA and not self.isRNA:
            raise ValueError("Sequence is not valid") 

        self.seq = sequence.upper()

    def read_from_file(self, filename):
        """ 
        Reads a sequence from a multi-line text file. 
        """
        fh = open(filename, "r")
        lines = fh.readlines()
        seq = ""
        for l in lines:
            seq += l.replace("\n","")
        fh.close()
        
        self.read_from_str(seq)

    def read_from_fasta(self, filename):
        """
        Reads a FASTA file to a dictionary
        """
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

        self.read_from_str(sequence)
        

        
    @staticmethod 
    def validate_dna(dna_seq):
        """ 
        Checks if DNA sequence is valid. Returns True is sequence is valid, or False otherwise. 

        :param dna_seq: string composed by sequence of DNA
        """
        seqm = dna_seq.upper()
        valid = seqm.count("A") + seqm.count("C") + seqm.count("G") + seqm.count("T")
        return valid == len(seqm)


    @staticmethod
    def validate_rna(rna_seq):
        """ 
        Checks if RNA sequence is valid. Returns True is sequence is valid, or False otherwise. 

        :param rna_seq: string composed by sequence of RNA
        """
        seqm = rna_seq.upper()
        valid = seqm.count("A") + seqm.count("C") + seqm.count("G") + seqm.count("U")
        return valid == len(seqm)

    def nucleotides_frequency(self, percentage: bool):
        """ 
        Calculates the frequency of each symbol in the sequence. Returns a dictionary. 

        :param percentage: Calculates relative frequency in percentage if true
        """
        if percentage:
            return {x: 100*self.seq.count(x) / len(self.seq) for x in set(self.seq)}

        else:
            return {x: self.seq.count(x) for x in set(self.seq)}

    def transcribe(self):
        """ 
        Function that computes the RNA corresponding to the transcription of the DNA sequence provided. 
        """
        if self.isRNA: raise ValueError("RNA can't be trancribed")
        return self.__class__.read_from_str(self.seq.replace("T","U"))

    def back_transcribe(self):
        """ 
        Function that computes the RNA corresponding to the transcription of the DNA sequence provided. 
        """
        if self.isDNA: raise ValueError("DNA can't be back trancribed")
        return self.__class__.read_from_str(self.seq.replace("U","T"))

    def reverse_complement(self):
        """ 
        Computes the reverse complement of the inputted DNA sequence. 
        """
        if self.isRNA: raise ValueError("There is no reverse complement of RNA sequence")

        rev = self.seq[::-1]
        complement_transfer = { 
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C"
        }
        comp_rev = ""
        for letter in rev:
            comp_rev += complement_transfer[letter.upper()]
            
        return comp_rev

    def gc_content(self):
        """ 
        Returns the percentage of G and C nucleotides in a DNA sequence. 
        """
        gc_count = 0
        for s in self.seq:
            if s in "GCgc": gc_count += 1
        return gc_count / len(self.seq)

    def gc_content_subseq(self, k=100):
        """ 
        Returns GC content of non-overlapping sub-sequences of size k. 
        """
        return [self.gc_content(self.seq[i:i+k]) for i in range(0, len(self.seq) - (k-1), k)]


    @staticmethod
    def translate_codon(cod: str):
        """
        Translates a codon into an aminoacid using an internal dictionary with the standard genetic code.

        :param cod: string of codon to be translated to AA
        """

        tc = {
            "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
            "TGT":"C", "TGC":"C",
            "GAT":"D", "GAC":"D",
            "GAA":"E", "GAG":"E",
            "TTT":"F", "TTC":"F",
            "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
            "CAT":"H", "CAC":"H",
            "ATA":"I", "ATT":"I", "ATC":"I",
            "AAA":"K", "AAG":"K",
            "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
            "ATG":"M", "AAT":"N", "AAC":"N",
            "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
            "CAA":"Q", "CAG":"Q",
            "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
            "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
            "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
            "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
            "TGG":"W",
            "TAT":"Y", "TAC":"Y",
            "TAA":"_", "TAG":"_", "TGA":"_"
        }
        if cod in tc: return tc[cod]
        else: return None

    def translate_seq(self, reverse_comp: bool = False, ini_pos: int = 0) -> list():
        """ 
        Translates a DNA sequence into an aminoacid sequence. 
        :param reverse_comp: bool that determines if ORF should be found in the negative direction
        :param ini_pos: ini_pos = 0 -> frame 1
                        ini_pos = 1 -> frame 2
                        ini_pos = 2 -> frame 3

        :returns: list of strings with aminoacid sequence for reading frame selected
        """
        dna_seq = self.seq if self.isDNA else self.back_transcribe().seq

        if reverse_comp: 
            dna_seq = self.reverse_complement()

        seq_aa = ""
        for start_triple in range(ini_pos, len(dna_seq)-(ini_pos+2),3):
            seq_aa += self.translate_codon(dna_seq[start_triple:start_triple+3])
        return seq_aa

    def codon_frequency(self,  reverse_comp: bool = False, ini_pos: int = 0) -> dict():
        """
        """
        seq = self.seq
        if reverse_comp: 
            seq = self.reverse_complement()
        freq = dict()

        for start_triple in range(ini_pos, len(seq)-(ini_pos+2),3):
            codon = seq[start_triple:start_triple+3]
            if codon in freq.keys():  freq[codon] += 1
            else: freq[codon] = 1

        return freq

    def min_max_codon_freq(self, ini_pos=Union[str, list]) -> dict():
        """
        Calculate minimum and maximum frquency codons for the reading frames selected by ini_pos
        :param ini_pos: str or list of str representing possible reading frames

        :returns: dictionary with min and max codon and their frequency for each rf as well as for
        the sum of all reading frames
        """
        ORF_freqs = {}
        total_freq = Counter({})
        for rf in ini_pos:
            cod_freq = self.codon_frequency(ini_pos=rf)
            max_codon = max(cod_freq , key = cod_freq.get)
            min_codon = min(cod_freq , key = cod_freq.get)
            ORF_freqs[rf] = {
                "max_codon": {max_codon: cod_freq[max_codon]}, 
                "min_codon": {min_codon: cod_freq[min_codon]}
            }
            total_freq += Counter(cod_freq)

        max_codon = max(total_freq, key = total_freq.get)
        min_codon = min(total_freq, key = total_freq.get)
        ORF_freqs["total"] = {
            "max_codon": {max_codon: total_freq[max_codon]}, 
            "min_codon": {min_codon: total_freq[min_codon]}
        }
        return ORF_freqs

    def codon_usage(self, aa: str):
        """
        Provides the frequency of each codon encoding a given aminoacid, in a DNA sequence.

        :param aa: aminoacid to find
        """
        seqm = self.seq if self.isDNA else self.back_transcribe().seq
        dic = {}
        total = 0
        for i in range(0, len(seqm)-2, 3):
            cod = seqm[i:i+3]
            if self.translate_codon(cod) == aa:
                if cod in dic:
                    dic[cod] += 1
                else: dic[cod] = 1
                total += 1
        if total >0:
            for k in dic:
                dic[k] /= total
        return dic

    def reading_frames(self):
        """
        Computes the six reading frames for a DNA sequence (including the reverse complement)
        And the 3 for a RNA sequence.
        """

        res = []
        res.append(self.translate_seq(ini_pos=0))
        res.append(self.translate_seq(ini_pos=1))
        res.append(self.translate_seq(ini_pos=2))

        if self.isDNA:
            res.append(self.translate_seq(reverse_comp = True, ini_pos=0))
            res.append(self.translate_seq(reverse_comp = True, ini_pos=1))
            res.append(self.translate_seq(reverse_comp = True, ini_pos=2))
        return res

    @staticmethod
    def all_proteins_rf(aa_seq: str, overlapping: bool = True, ini_pos: int = 0):
        """
        Computes all posible proteins in an aminoacid sequence.

        :params aa_seq: string representing aa sequence
        :param overlapping: boolean indicating if overlapping proteins in the same
            region should be computed, if false only keeps the biggest protein
            of each region
        """
        aa_seq = aa_seq.upper()
        current_prot = []
        proteins = []
        # Lista com coordenadas (Start, End)
        coords = []
        # Variável auxiliar p/ guardar coordenadas
        # (vai corresponder à coordenada de início)
        x = ini_pos
        # Counter para saber em que coordenada estamos
        counter = ini_pos
        # Flag para saber se uma nova coordenada já foi
        # atribuida a x. Isto porque, se a sequência de AA
        # for MAMFT, x primeiro ficaria com a coordenada do
        # primerio M, mas depois teria a coordenada do segundo M
        # Nós só queremos a coordenada do primeiro M que
        # corresponderá à maior sequência
        flag = 0
        for aa in aa_seq:
            if aa == "_":
                if current_prot:
                    if overlapping:
                        proteins.append(p for p in current_prot)
                    else: 
                        proteins.append(max(current_prot, key=len))
                    coords.append((x, counter + 3))
                    current_prot = []
                    flag = 0
            else:
                if aa == "M":
                    current_prot.append("")
                    if flag == 0:  # Se estivermos num aa inicial (M)
                        # e ainda não tivermos atribuido valor a x
                        x = counter  # Então x passa a ter o valor da coordenada atual
                        flag = 1  # E a flag fica ativa
                for i in range(len(current_prot)):
                    current_prot[i] += aa
            counter += 3
        return proteins, coords

    def orf(
        self, 
        ini_pos: int = 0,  
        reverse_comp: bool = False, 
        overlapping: bool = True,
        minsize: int = 0
    ) -> list():
        """
        Find proteins in specific reading frame

        :param reverse_comp: bool that determines if ORF should be found in the negative direction
        :param ini_pos: ini_pos = 0 -> frame 1
                        ini_pos = 1 -> frame 2
                        ini_pos = 2 -> frame 3
        :param overlapping: boolean indicating if overlapping proteins in the same
            region should be computed, if false only keeps the biggest protein
            of each region
        :param minsize: minimum size of reading of ORF

        :returns: list of strings with aminoacid sequence for each protein found
        """
        rf = self.translate_seq(reverse_comp = reverse_comp, ini_pos=ini_pos)
        proteins, coords = self.all_proteins_rf(rf, overlapping, ini_pos)
        res = []
        # Counter para saber em que orf estou
        orf_counter = 0
        # Lista de listas com coordenadas e orf
        coords_all = []
        for p, c in zip(proteins, coords):
            if len(p) >= minsize:
                orf_counter += 1
                #self.insert_prot_ord(p, res)
                res.append(p)
                coords_all.append((c[0], c[1], orf_counter))

        return res, coords_all

    def all_orfs(self, overlapping: bool = True, minsize: int = 0):
        """
        Computes all possible proteins for all open reading frames.

        :param overlapping: boolean indicating if overlapping proteins in the same
            region should be computed, if false only keeps the biggest protein
            of each region
        """
        rfs = self.reading_frames()
        res = []
        # Counter para saber em que orf estou
        orf_counter = 0
        # Counter para saber qual o reading frame
        rf_counter = 0
        # Lista de listas com coordenadas e orf
        coords_all = []
        for rf in rfs:
            prots, coords = self.all_proteins_rf(rf, overlapping, rf_counter)
            rf_counter += 1
            for p, c in zip(prots, coords):
                if len(p) >= minsize:
                    orf_counter += 1
                    #self.insert_prot_ord(p, res)
                    res.append(p)
                    coords_all.append((c[0], c[1], orf_counter))

        return res, coords_all

    def all_orfs_ord(self, overlapping: bool= True, minsize = 0):
        """
        Computes all possible proteins for all open reading frames. 
        Returns ordered list of proteins with minimum size.

        :param overlapping: boolean indicating if overlapping proteins in the same
            region should be computed, if false only keeps the biggest protein
            of each region
        """
        orfs, coords = self.all_orfs(overlapping, minsize)
        for i in range(len(orfs)):
            orfs[i].sort(key = len, reverse = True)

        return orfs

    def __str__(self):
        return self.seq

    def __repr__(self):
        return self.__str__()


if __name__ == "__main__":
    import sys
    import csv

    if len(sys.argv) != 2:
        print("Usage: python3 sequence.py <filename>")

        exit()
        
    seq_file = sys.argv[1]

    covid_sequence = genome()
    covid_sequence.read_from_fasta(
        seq_file
    )

    print(
        "Exercise 1. \n Sequence Length:",
        len(covid_sequence.seq)
    )

    print(
        "Exercise 2. \n Nucleotides Frequency:",
        covid_sequence.nucleotides_frequency(percentage=True)
    )

    print(
        "Exercise 3. \n GC Content:",
        covid_sequence.gc_content()
    )

    print(
        "Exercise 4. \n Nr Start Codons in Poistive Reading Frames:",
        sum( [covid_sequence.translate_seq(ini_pos = pos).count("M") for pos in range(0,3)] )
    )

    print(
        "Exercise 5. \n Nr Stop Codons in Poistive Reading Frames:",
        sum( [covid_sequence.translate_seq(ini_pos = pos).count("_") for pos in range(0,3)] )
    )

    ORF_freqs = covid_sequence.min_max_codon_freq(ini_pos=[0,1,2])
    print(
        "Exercise 6. Codons Frequency:",
        f""" 
            ORF | Most Frequent | Less Frequent
            0 | {ORF_freqs[0]['max_codon']} | {ORF_freqs[0]['min_codon']}
            1 | {ORF_freqs[1]['max_codon']} | {ORF_freqs[1]['min_codon']}
            2 | {ORF_freqs[2]['max_codon']} | {ORF_freqs[2]['min_codon']}
            Total | {ORF_freqs['total']['max_codon']} | {ORF_freqs['total']['min_codon']}
        """
    )
    rfs = []
    coords = []
    for ini_pos in [0,1,2]:
        rf, coord = covid_sequence.orf(ini_pos, overlapping=False, minsize=40)
        rfs += rf
        coords += coord

    proteins = []
    with open("all_potential_proteins.txt", "w") as outfile:
        outfile.writelines(protein + "\n" for protein in rfs)
    print(
        "\nExercise 7.",
        f"""
        {len(rfs)} proteins with {40} or more aminoacids were found
        The file all_potential_proteins.txt was written
        """
    )
    
    c = 0
    with open("orf_coordinates.txt", "w") as outfile:
        for x in coords:
            c += 1
            outfile.writelines(str(x[0]) + ', ' + str(x[1]) + ', ' + str(c) + "\n")
    print(
        "\nExercise 8.",
        """
        The file orf_coordinates.txt was written
        """
    )

    # opening the file using "with" 
    # statement
    overlaps = []
    with open("proteins_86693_757732.csv", 'r') as data:
        for protein in csv.DictReader(data):
            max_overlap = 0
            for coord in coords:
                if (int(protein["Start"]) <= coord[1]) and (int(protein["Stop"]) >= coord[0]):
                    overlap = ( min(int(protein["Stop"]), coord[1]) - max(int(protein["Start"]), coord[0]) ) // 3 / int(protein["Length"])
                    if overlap > max_overlap: max_overlap = overlap
            overlaps.append((protein["Locus"], f"{int(round(max_overlap*100, 0))}%"))

    # with open("orf_overlaps.txt", "w") as outfile:
    #     for overlap in overlaps:
    #         outfile.writelines(f"{overlap[0], overlap[1]}\n")


    print('\nExercise 9')
    for o in overlaps:
        print('', o)
