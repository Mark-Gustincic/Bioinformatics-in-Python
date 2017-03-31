#!/usr/bin/env python3
# Name: Mark Gustincic
# Group Members: Thomas Richards (tarichar)

class Converter (str):
    def __new__(self,name):
        return str.__new__(self,name.upper())

AA = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
"""
Flips keys and values of dictionary AA and creates new dictionary OneToThree.
"""
OneToThree = {value:key for key,value in AA.items()}
RNA_codon_table = {
# Second Base
# U             C             A             G
#U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
#C 
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
#A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
#G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
}
dnaCodonTable = {key.replace('U','T'):value for key, value in RNA_codon_table.items()}

"""
Asks for string input and calls it name, then instantiates name1 as instance
of converter class.
"""
name = input("Please input a single string: ")
name1 = Converter(name)

if name1 in OneToThree:
    """
    If input string is located in dictionary OneToThree, prints value associated
    with key (obtained from input).
    """
    print (name1,"=",OneToThree[name1])
elif name1 in AA:
    """
    If input string is located in dictionary AA, prints value associated
    with key (obtained from input).
    """
    print (name1,"=",AA[name1])
elif name1 in RNA_codon_table:
    """
    If input string is located in dictionary RNA_codon_table, prints value
    associated with key (obtained from input).
    """
    print (name1,"=",RNA_codon_table[name1])
elif name1 in dnaCodonTable:
    """
    If input string is located in dictionary dnaCodonTable, prints value associated
    with key (obtained from input).
    """
    print (name1,"=",dnaCodonTable[name1])
else:
    """
    If input string is not found in any included dictionaries, returns the
    string "unknown".
    """
    print ("unknown")


    
    
