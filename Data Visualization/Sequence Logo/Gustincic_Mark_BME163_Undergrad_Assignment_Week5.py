import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import sys
import matplotlib.image as mpimg


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



plt.style.use('BME163.mplstyle')#if throws error, remove '.mplstyle'

plt.figure(figsize=(8,3))

#panel=plt.axes([left, bottom, width, height])
panel1=plt.axes([.075,.2,.42,.6])
panel2=plt.axes([.52,.2,.42,.6])

yBins = np.arange(0, 2.5, .5)
xBins = np.arange(-10, 15, 5)

panel1.set_xticks(xBins)
panel1.set_yticks(yBins)
panel2.set_xticks(xBins)
panel2.set_yticks(yBins)

logoA = mpimg.imread('A.png')
logoT = mpimg.imread('T.png')
logoG = mpimg.imread('G.png')
logoC = mpimg.imread('C.png')

three = []
five = []
freqP3 = []
freqP5 = []

if sys.argv[1]:
    inf = sys.argv[1]
else:
    print('No input file',file = sys.stderr)

myReader = FastAreader (inf)

'''
The below loop iterates through the input file by each header/sequence pair.
The header is stripped & split on ' to obtain direction of read(5 or 3).
After direction is identified, the list(sequence) is appended to its appropriate list(three or five).
Liblist contains the 3 and 5 lists(list of lists).
'''

for head, seq in myReader.readFasta():
    split_head=head.strip().split("'")
    direction = split_head[0]
    seqL = list(seq)
    if direction == '3':
        three.append(seqL)
    else:
        five.append(seqL)


libList = [three, five]

'''
The below loop iterates through each library(3 & 5).
List of 20 lists is initialized to hold all the nucleotides together from the same
position in sequence.
Position is also initialized at 0 for each library.


The first nested loop(1) iterates through each sequence(list) present in the library.
Each nucleotide is isolated and appended to its appropriate position list within the
posCountLists list.
The conditional is included to prevent the position from iterating past the
length of the sequence.

frequencyList is itintialized as a list of 20 lists each containing [0,0,0,0]
and will be used to store the counts of each nucleotide at each position.

The next nexted loop(1) iterates through each position(list containing nucleotides)
held in posCountLists.
A,T,G,C are counted at each position and value is added to nt's respective locations in
position list contained by frequencyList

'''
#Values acquired are about 1/20th the expected size.
#Sum of values at each position should match number of sequences in library
'''

'''
libCount = 0
for library in libList:
    en = ((1/np.log(2))*(3/(2*len(library))))
    posCountLists = [[] for i in range(20)]
    for seqList in library:
        count = 0
        for base in seqList:
            posCountLists[count].append(base)
            count += 1

    pos2 = 0
    frequencyList = [[0,0,0,0] for i in range(20)]
    
    for position in posCountLists:
        A = position.count('A')
        T = position.count('T')
        G = position.count('G')
        C = position.count('C')

        frequencyList[pos2][0] += A
        frequencyList[pos2][1] += T
        frequencyList[pos2][2] += G
        frequencyList[pos2][3] += C
        pos2 += 1
        
    #print(frequencyList)

    for position in frequencyList:
        freqA = (position[0])/len(library)
        freqT = (position[1])/len(library)
        freqG = (position[2])/len(library)
        freqC = (position[3])/len(library)
        baseList = [freqA, freqT, freqG, freqC]
        if libCount == 0:
            freqP3.append(baseList)
        else:
            freqP5.append(baseList)
    libCount += 1

posCount = 0
for freqList in freqP3:
    baseCount = 0
    unc = -1*((freqList[0]*np.log2(freqList[0]))+(freqList[1]*np.log2(freqList[1]))+(freqList[2]*np.log2(freqList[2]))+(freqList[3]*np.log2(freqList[3])))
    Ri = np.log2(4) - (unc + en)
    x = 0
    for myFreq in freqList:
        height = myFreq*Ri
        freqList[x] = height
        x += 1

posCount = 0
for freqList in freqP5:
    baseCount = 0
    unc = -1*((freqList[0]*np.log2(freqList[0]))+(freqList[1]*np.log2(freqList[1]))+(freqList[2]*np.log2(freqList[2]))+(freqList[3]*np.log2(freqList[3])))
    Ri = np.log2(4) - (unc + en)
    y = 0
    for myFreq in freqList:
        height = myFreq*Ri
        freqList[y] = height       
        y += 1
    
counter = 0
posList3 =[[[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]]]

posList5 = [[[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]],
           [[logoA,0],[logoT,0],[logoG,0],[logoC,0]]]

p = 0
for freqList in freqP3:
    b = 0
    for height in freqList:
        posList3[p][b][1] = height
        b += 1
        if b > 3:
            continue
    p += 1
    if p > 19:
        break
    
p = 0
for freqList in freqP5:
    b = 0
    for height in freqList:
        posList5[p][b][1] = height
        b += 1
        if b > 3:
            continue
    p += 1
    if p > 19:
        break
        
counter3 = 0        
for position in posList3:
    previous_height=0        
    for base in sorted(position, key = lambda x:x[1]):
        panel2.imshow(base[0], extent=[counter3-10, counter3-9, previous_height, base[1]+previous_height], aspect='auto', alpha=1, zorder=1000) #left, right, bottom, top
        previous_height += base[1]
    counter3 += 1

counter5 = 0
for position in posList5:
    previous_height=0
    for base in sorted(position, key = lambda x:x[1]):
        panel1.imshow(base[0], extent=[counter5-10, counter5-9, previous_height, base[1]+previous_height], aspect='auto', alpha=1, zorder=1000)
        previous_height += base[1]
    counter5 += 1

panel1.set_xlabel('Distance to \nSplice Site', fontsize=10)
panel1.set_ylabel('Bits', fontsize=10)
panel2.set_xlabel('Distance to \nSplice Site', fontsize=10)


panel1.plot([0, 0], [0, 2], color='black', linestyle='-', linewidth=.5)
panel2.plot([0, 0], [0, 2], color='black', linestyle='-', linewidth=.5)

panel1.set_title("5'SS")
panel2.set_title("3'SS")

panel1.tick_params(axis='both',which='both',\
                   bottom='on', labelbottom='on',\
                   left='on', labelleft='on',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')
panel2.tick_params(axis='both',which='both',\
                   bottom='on', labelbottom='on',\
                   left='off', labelleft='off',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')

panel1.set_xlim(-10,10)
panel1.set_ylim(0,2)
panel2.set_xlim(-10,10)
panel2.set_ylim(0,2)
    

plt.savefig('Gustincic_Mark_BME163_Undergrad_Assignment_Week5.pdf')






    
