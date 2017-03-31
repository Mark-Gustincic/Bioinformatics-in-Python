#!/usr/bin/env python3 
# Name: Mark Gustincic
# Group Members: Thomas Richards (tarichar)
'''
Pseudo Code for 

class OrfFinder():
    def __init__(self, seq):
        initialize: sequence(upper case), list of start codons, list of stop codons,
                    dictionary of nucleotide complements. List of found start
                    codons with an initial value of 0. Set of all found ORFs in sequence.

        set bases equal to a list of the values present in self.seq.
        set self.reverseBases to the nucleotide comeplement of self.seq and join list values
        together into a sequence called self.reverseSeq.
      
    def findOrfs (self):
        loop through frames 0, 1, and 2
            set self.startList to contain only 0
            loop through all points within reading frame
                set codon to the value of the point to the point + 3
                conditional for codon in self.start (list of starts)
                    add start positions to self.startList
                    set mimimum start to the start w/ smallest start coordinate + 1
                conditional for codon in self.stop (list of stops)
                    loop through starts in self.startList
                        save frame, start, stop and length
                    blank list of found start codons
                conditional for position >= end of sequence - 2
                    conditional if start list contains any values
                        same frame, start, stop, and length
                conditional for anything else
                    pass

    def findRevOrfs (self):
        loop through frames 0, 1, and 2
            set self.startList to contain only 0
            loop through all points within reading frame
                set codon to the value of the point to the point + 3
                conditional for codon in self.start (list of starts)
                    add start positions to self.startList
                    set mimimum start to the start w/ smallest start coordinate + 1
                conditional for codon in self.stop (list of stops)
                    loop through starts in self.startList
                        save frame, start, stop and length
                    blank list of found start codons
                conditional for position >= end of sequence - 2
                    conditional if start list contains any values
                        same frame, start, stop, and length
                conditional for anything else
                    pass

    def saveORFs (self):
        conditional if on positive strand
            add frame, start, stop & frame to set of ORFs

        condtionatal if on negative strand
            set a to value of start
            set b to value of stop
            add frame, start, stop & frame to set of ORFs

    def ORFsorter (self):
        set ORFs to values in a list of ORFs
        sort ORFs by their 2nd value, which is starting position
        sort ORFs by their 4th value, which is ORF length
        reverse order of ORFs so that they are listed from longest to shortest   
'''

import sys

class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        self.parser.add_argument('outFile', action = 'store', help='output file name') 
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=True, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= range(0, 2000), action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', nargs='?', help='start Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

class OrfFinder():
    '''
    The __init__ method initializes lists & dictionaries to be used later in the code. Start & stop lists
    contain the start codons to be used when searching for ORFs. The comp dict is used to find the nucleotide
    complement in order to produce the complementary starnd. Startlist is set to 0 to handle dangling start codons.
    ORFlist is actually a set to account for duplicate ORFs. Bottom section of method creates the complementary strand.
    '''
    def __init__(self, seq):
        self.seq = seq.upper()
        self.start = ['ATG']
        self.stop = ['TAA','TGA','TAG']
        self.compDict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
        self.startList = [] #Instantiates a list of all start codon coordinates present before stop codon.
        self.ORFlist = set()

        bases = list(self.seq)
        self.reverseBases = [self.compDict[nt] for nt in reversed(bases)]
        self.reverseSeq = ''.join(self.reverseBases)    

    def findOrfs (self):
        '''
        The findORFs method finds ORFs in the top strand of an inputted sequence. Once the ORFs are found, they are saved in
        a list of ORFs to be sorted in ORFsorter. This method handles both the dangling start and stop problems.
        '''
        for frame in (0,1,2): 
            self.startList = []
            for p in range(frame, len(self.seq), 3):
                codon = self.seq[p : p + 3]
                if codon in self.start:
                    self.startList.append(p)#adds start codon coordinate to a list
                    minStart = min(self.startList) + 1 
                elif codon in self.stop:
                    for thisStart in self.startList:
                        self.saveORF(frame+1, minStart, p + 3, p - minStart + 4)
                    self.startList = []
                elif p >= len(self.seq)-2:
                    if self.startList: #dangling start
                        self.saveORF(frame+1, self.minStart, len(self.seq), len(self.seq)- (minStart + 4))
                #else:
                    #pass
                
                                                
    def findRevOrfs (self):
        '''
        The findORFs method finds ORFs in the complementary strand of an inputted sequence. Once the ORFs are found, they are saved in
        a list of ORFs to be sorted in ORFsorter. This method handles both the dangling start and stop problems.
        '''
        
        for frame in (0,1,2):
            self.startList = []
            for q in range(frame, len(self.reverseSeq), 3):
                codon = self.reverseSeq[q : q + 3]
                if codon in self.start:
                    self.startList.append(q)#adds start codon coordinate to a list
                    minStart = min(self.startList) + 1 
                elif codon in self.stop:
                    for thisStart in self.startList:
                        self.saveORF(-(frame+1), minStart, q + 3, q - minStart + 4)
                    self.startList = []    
                elif q >= len(self.seq)-2:
                    if self.startList: #dangling stop
                        self.saveORF(-(frame+1), minStart, len(self.reverseSeq), len(self.reverseSeq)- (minStart + 4))
                #else:
                    #pass
                                
                
                    
    def saveORF (self, frame, start, stop, length):#takes ORF and puts in list of many ORFs to be sorted
        '''
        The saveORF method saves ORFs found in the sequence to a list of ORFs.
        '''
        if frame > 0: #Top strand ORFs
            self.ORFlist.add((frame, start, stop, length))

        if frame < 0: #Bottom strand ORFs
            a = len(self.reverseSeq) - stop #finds start position according to top strand coordinates
            b = len(self.reverseSeq) - start #finds stop position according to bottom strand 
            self.ORFlist.add((frame, a + 1, b + 1, length))

            
    def ORFsorter (self):
        '''
        The orfSorter method takes in a set of ORFs, and sorts them according to start position and length.
        First, the contents of self.ORFlist(set) are made into an list. ORFs are then sorted by their start
        position. An optional, commented out line of code gives the programmer the ability to reverse the order
        of the sorted list. Next, the ORFs are sorted by their length, and then the list is reversed to sort the
        ORFs from longest to shortest.
        '''
        ORFs = list(self.ORFlist)
        ORFs.sort(key=lambda entry: entry[3])
        ORFs.reverse()
        return ORFs

   
    
                
class FastAreader :
	
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        
        self.fname = fname
            
    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
 
    def readFasta (self):
		
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
			
            header = ''
            sequence = ''
 
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()
 
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
						
        yield header,sequence


def main(myCommandLine=None):
    '''
    Implements the Usage exception handler that can be raised from anywhere in process.
    If - else conditional below denotes whether or not to use the command line.
    If command line is empty, the program is to use the default parameters.
    Else the program will listen to the command line. 

    myReader is initialized to run the inFile through FastAreader.
    Loop through head and sequence in file. Set genome to instance of OrfFinder
    with input file as input. The program then searches for ORFs on the top and
    bottom strand. The set 'seen' gets initialized as an empty set. Program then
    loop through ORFs in genome.ORFsorter. If the ORF length is greater than 100,
    print frame, start, stop, and length to an output file.
    '''
    
    if myCommandLine is None:
        myCommandLine = CommandLine([ 'tass2.fa',
                                      'tass2ORFdata-ATG-100.txt',
                                      '--longestGene',
                                      '--start=ATG',
                                      '--minGene=100'])
    else :
        myCommandLine = CommandLine(myCommandLine)

    myReader = FastAreader (myCommandLine.args.inFile)
    outputFile = open(myCommandLine.args.outFile, 'w')
    for head, seq in myReader.readFasta():
        genome = OrfFinder(seq)
        genome.start = myCommandLine.args.start
        genome.findOrfs()
        genome.findRevOrfs()
        print (head, file = outputFile)
        for ORF in genome.ORFsorter():
            if ORF[3] >= myCommandLine.args.minGene:
                print ('{0:+d} {1:>5d}..{2:>5d} {3:5d}'.format(ORF[0], ORF[1], ORF[2], ORF[3]), file=outputFile)
    outputFile.close()

if __name__ == "__main__":
    main()



    
    

    
    

        
        
        
            
                       
        
            
                             
                

        
    
