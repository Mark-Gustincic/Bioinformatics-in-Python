#!/usr/bin/env python3
# Name: Mark Gustincic
# Group Members: Thomas Richards (tarichar)


class NucParams (str):
    """
    The NucParams class is designed to take in a string of nucleotides and
    returns certain properties associated with the sequence. Such properties
    include the number of amino acids, the number of nucleotides, the number of
    codons, and the protein's amino acid composition expressed as percentages.
    
    The below section of code initializes dictionaries and value names as class
    attributes for use in methods within the ProteinParam class.
    """
    def __init__(self, dna):
        """
        The rnaCodonTable is a dictionary used to relate codons with their
        single letter code. The dnaCodonTable is created by replacing all
        Us with Ts.
        """
        self.rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
        }
        self.dnaCodonTable = {key.replace('U','T'):value for key, value in self.rnaCodonTable.items()}
        """
        nucString, dnaString, and rnaString are all initialized here. dnaString
        and rnaString are the same, except with the replacement of T for U in
        rnaString.
        """
        self.nucString = dna.upper()
        self.dnaString = self.nucString.replace('U', 'T')
        self.rnaString = self.nucString.replace('T', 'U')
        """
        The list of appropriate nucleotides is initialized here, along with the
        protein, codon, and nucleotide dictionaries. The protein & codon dictionaries
        are keyed with the values & keys of the rnaCodonTable, respectively. The
        nucleotide dictionay is keyed with the components of the nucleotide list.
        """

        nucList = ['A','C','G','T','U','N']
        self.proteinDict = {key:0 for key in self.rnaCodonTable.values()}
        self.codonDict = {key:0 for key in self.rnaCodonTable.keys()}
        self.ntDict = {key:0 for key in nucList}
               
    def addSequence (self, thisSequence):
        """
        The addSequence method is responsible for creating codons out of the
        inputted sequence (thisSequence). First, it cuts the sequence into 3
        letter blocks, then converts the sequence from DNA to RNA.

        The codon list is then iterated twice times, the first time counting the number
        of codons in sequence, and adding the counts to the codon dictionary. The
        second time, the number of amino acids is counted and added to the protein
        dictionary. The inputted sequence (thisSequence) is then iterated and each
        nucleotide is counted and the counts are adde to the nucleotide dictionary.
        """
        codonList = [thisSequence[i:i+3] for i in range(0, len(thisSequence), 3)]
        for codon in codonList:
            codon = codon.replace('T', 'U')
            if codon in self.rnaCodonTable.keys():
                self.codonDict[codon] += 1
                
            aa = self.rnaCodonTable[codon]
            if aa in self.proteinDict.keys():
                self.proteinDict[aa] += 1

        for nt in thisSequence:
             self.ntDict[nt] += 1

    def aaComposition(self):
        """
        aaComposition returns the dictionary of amino acids.
        """
        return self.proteinDict
            
    def nucComposition(self):
        """
        nucComposition returns the dictionary of nucleotides.
        """
        return self.ntDict
    
    def codonComposition(self):
        """
        codonComposition returns the dictionary of codons.
        """
        return self.codonDict
        
    def nucCount(self):
        """
        nucCount returns the counts of all nucleotides in the rnaString
        and returns the sum of all the values.
        """

        return sum(self.ntDict.values())
        
        
# Please do not modify any of the following.  This will produce a standard output that can be parsed



class FastAreader :
    '''
    Class to provide reading of a file containing one or more FASTA
    formatted sequences:
    object instantiation:
    FastAreader(<file name>):
 
    object attributes:
    fname: the initial file name
 
    methods:
    readFasta() : returns header and sequence as strings.
    Author: David Bernick
    Date: April 19, 2013
    '''
    def __init__ (self, fname):
        '''contructor: saves attribute fname '''
        self.fname = fname
 
    def readFasta (self):
        '''
        using filename given in init, returns each included FastA record
        as 2 strings - header and sequence.
        whitespace is removed, no adjustment is made to sequence contents.
        The initial '>' is removed from the header.
        '''
        header = ''
        sequence = ''
        
        with open(self.fname) as fileH:
            # initialize return containers
            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            # header is saved, get the rest of the sequence
            # up until the next header is found
            # then yield the results and wait for the next call.
            # next call will resume at the yield point
            # which is where we have the next header
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        # final header and sequence will be seen with an end of file
        # with clause will terminate, so we do the final yield of the data
        yield header,sequence
 
# presumed object instantiation and example usage
# myReader = FastAreader ('testTiny.fa');
# for head, seq in myReader.readFasta() :
#     print (head,seq)
