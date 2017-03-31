#!/usr/bin/env python3
# Name: Mark Gustincic
# Group Members: Thomas Richards (tarichar)

class SeqCleaner:
    """
    Returns new string in upper case.
    """
    def __init__(self,s):
        self.DNAstring = s.upper()
        self.compDict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    """
    Counts number of "N"s in inputted string, and splits the
    string into sections i, j, and k. Section containing "N"
    is replaced by {num_N}, and three parts of string are joined
    together.s
    """
    def clean(self):
        ntList = list(self.DNAstring)
        count = 0
        for base in ntList:
            if not base in self.compDict.keys():
                ntList[count] = '-'
            count += 1
        cleanStr = ''.join(ntList)
        return(cleanStr.lower())
    """
    Prompts user to input dna sequence, and instantiates dnaStr as an
    instance of SeqCleaner using argument dDna. dnaStr is printed to output
    after being run through the method clean.
    """
dDna = input("What is your DNA sequence?\n")
dnaStr = SeqCleaner(dDna)
cleanDNA = dnaStr.clean()
print(cleanDNA)
    
    
