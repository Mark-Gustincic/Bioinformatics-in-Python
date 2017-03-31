#!/usr/bin/env python3
# Name: Mark Gustincic
# Group Members: Thomas Richards (tarichar)
"""
genomeAnalyzer organizes the output for the sequenceAnalysis module. It is
responsible for importing sequenceAnalysis and reading the file 'testGenome.fa'.
It also finds the relative usage of each codon, along with the gc content of
sequence. testGenome is iterated, and each codon is counted, along with the
length of the sequence in mega bases.
"""
import sequenceAnalysis
myReader = sequenceAnalysis.FastAreader('testGenome.fa')
seq = ''
nuc = sequenceAnalysis.NucParams(seq)
for head, seq in myReader.readFasta():
    nuc.addSequence(seq)

length = nuc.nucCount()
mbLength = length/1000000
g = nuc.ntDict['G']
c = nuc.ntDict['C']
gcContent = ((g+c)/length)*100

print ('sequence length = {:.2f} Mb'.format(mbLength))
print (' ')
print ('GC content = {:.1f}%'.format(gcContent))
print (' ')
codonList = []
for codon in nuc.rnaCodonTable.keys():
    codonList += [codon[i:i+3] for i in range(0, len(codon),3)]
codonList.sort()
for codon in codonList:
    if codon in nuc.rnaCodonTable.keys():
        aa = nuc.rnaCodonTable[codon]
        w = nuc.codonDict[codon]#codon count
        b = nuc.proteinDict[aa]#aa count
        z = (w/b)*100
        print ('{0} : {1} {2:5.1f} ({3:6d})'.format(codon, aa, z, w))




