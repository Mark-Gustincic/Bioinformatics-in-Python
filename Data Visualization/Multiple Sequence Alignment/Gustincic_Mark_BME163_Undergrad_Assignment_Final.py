import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mplpatches

plt.style.use('BME163.mplstyle')
fig_1 = plt.figure(figsize=(10,5))
panel1 = plt.axes([0, .0625, 1, .25],frameon=True)
panel2 = plt.axes([0, .3625, 1, .25],frameon=True)
panel3 = plt.axes([0, .6625, 1, .25],frameon=True) 

left_bound=45232945
right_bound=45240000

ill_read_list=[]
for line in open('BME163_Input_data5.psl'):
    a=line.strip().split('\t')
    chromosome=a[13]
    if chromosome=='chr7':
        start=int(a[15])
        end=int(a[16])
        blockstarts=np.array(a[20].split(',')[:-1],dtype=int)
        blocksizes=np.array(a[18].split(',')[:-1],dtype=int)
        if start > left_bound and end < right_bound:
            ill_read_list.append([start,end,blockstarts,blocksizes,'not plotted'])
    

illSort = sorted(ill_read_list,key=lambda x: x[1])#length = 336

lineListIll = [] #lineList will be list containing all readLists
readPos = 1
for read in illSort:
    if read[4] == 'not plotted':
        readList = [] #readList will contain a list of all reads on same line
        lineEnd = read[1]
        readList.append(read)
        read[4] = 'plotted'

        for nextRead in illSort[readPos:]:
            if nextRead[4] == 'not plotted' and nextRead[0] >= lineEnd:
                lineEnd = nextRead[1]
                readList.append(nextRead)
                nextRead[4] = 'plotted'
        lineListIll.append(readList)
    readPos += 1

y_pos = 0
for readList in lineListIll:
    y_pos += 1
    for read in readList:
        blockstarts=read[2]
        blocksizes=read[3]
        rectangle1=mplpatches.Rectangle((read[0],y_pos),read[1]-read[0],0.125,facecolor='black',linewidth=0)
        panel2.add_patch(rectangle1)
        for pos in range(0,len(blockstarts),1):
            rectangle1=mplpatches.Rectangle((blockstarts[pos],y_pos-.2),blocksizes[pos],0.525,facecolor='black',linewidth=0)
            panel2.add_patch(rectangle1)

nano_read_list=[]
for line in open('BME163_Input_data4.psl'):
    a=line.strip().split('\t')
    chromosome=a[13]
    if chromosome=='chr7':
        start=int(a[15])
        end=int(a[16])
        blockstarts=np.array(a[20].split(',')[:-1],dtype=int)
        blocksizes=np.array(a[18].split(',')[:-1],dtype=int)
        if start > left_bound and end < right_bound:
            nano_read_list.append([start,end,blockstarts,blocksizes,'not plotted'])

nanoSort = sorted(nano_read_list,key=lambda x: x[1])
lineList = [] #lineList will be list containing all readLists
readPos = 1
for read in nanoSort:
    if read[4] == 'not plotted':
        readList = [] #readList will contain a list of all reads on same line
        lineEnd = read[1]
        readList.append(read)
        read[4] = 'plotted'

        for nextRead in nanoSort[readPos:]:
            if nextRead[4] == 'not plotted' and nextRead[0] >= lineEnd:
                lineEnd = nextRead[1]
                readList.append(nextRead)
                nextRead[4] = 'plotted'
        lineList.append(readList)
    readPos += 1

y_pos = 0
for readList in lineList:
    y_pos += 1
    for read in readList:
        blockstarts=read[2]
        blocksizes=read[3]
        rectangle1=mplpatches.Rectangle((read[0],y_pos),read[1]-read[0],0.125,facecolor='black',linewidth=0)
        panel1.add_patch(rectangle1)
        for pos in range(0,len(blockstarts),1):
            rectangle1=mplpatches.Rectangle((blockstarts[pos],y_pos-.2),blocksizes[pos],0.525,facecolor='black',linewidth=0)
            panel1.add_patch(rectangle1)

exon_read_list = []
transcript_id_list = []
counter = 0
for line in open('gencode.vM12.annotation.gtf.txt'):
    if counter < 5:
        counter += 1
        continue
    a = line.strip().split('\t')
    if a[0] == 'chr7' and a[2] == 'exon':
        start = int(a[3])
        end = int(a[4])
        exon_length = end - start
        if left_bound<start<right_bound and end<right_bound:                
            attributes = a[8].strip().split('; ')
            x = attributes[1].split('"')
            transcript_id = x[1]
            if transcript_id not in transcript_id_list:
                transcript_id_list.append(transcript_id)
            exon_read_list.append((start,end,exon_length,transcript_id))

transcriptDict = {key: [] for key in transcript_id_list}
lengthDict = {key: [] for key in transcript_id_list}
tsEndList = [] #list containing the endpoint of each transcript
tsStartList = []

for transID in transcriptDict.keys():
    startList = []
    endList = []
    for entry in exon_read_list:
        if entry[3] == transID:
            transcriptDict[transID].append(entry)
    for exon in transcriptDict[transID]:
        startList.append(exon[0])
        endList.append(exon[1])
        

    tsStart = min(startList)
    tsStartList.append(tsStart)
    tsEnd = max(endList)
    tsEndList.append(tsEnd)
    transcript_length = tsEnd - tsStart
    lengthDict[transID] = [tsStart, tsEnd, transcript_length]

tsEndSorted = sorted(tsEndList)


for ID, exons in transcriptDict.items(): #exons is a list of lists that contain exon parameters
    for exon_params in exons:
        exStart = exon_params[0]
        exEnd = exon_params[1]
        if exStart in tsStartList:
            tsStart = exStart
        if exEnd in tsEndList:
            tsEnd = exEnd
        height = tsEndSorted.index(tsEnd) + 1
        if height > 3:
            height = tsEndSorted.index(tsEnd)

        exon = mplpatches.Rectangle((exStart, height-.2), exEnd - exStart, 0.5, facecolor='black',linewidth=.6)
        if len(exons) == 2:
            exon = mplpatches.Rectangle((exStart, .8), exEnd - exStart, 0.5, facecolor='black',linewidth=.6)
        panel3.add_patch(exon)


    intron=mplpatches.Rectangle((tsStart,height),tsEnd - tsStart,0.1,facecolor='black',linewidth=.1)
    if exStart == max(tsStartList):
        intron=mplpatches.Rectangle((tsStart,1),tsEnd - tsStart,0.1,facecolor='black',linewidth=.1)
    panel3.add_patch(intron)

panel1.set_xlim(left_bound,right_bound)
panel1.set_ylim(0,y_pos+50)
panel2.set_xlim(left_bound,right_bound)
panel2.set_ylim(0,len(lineListIll)+5)
panel3.set_xlim(left_bound,right_bound)
panel3.set_ylim(0,len(tsEndSorted))

panel1.tick_params(axis='both',which='both',\
                   bottom='off', labelbottom='off',\
                   left='off', labelleft='off',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')

panel1.xaxis.set_visible(False)

panel2.tick_params(axis='both',which='both',\
                   bottom='off', labelbottom='off',\
                   left='off', labelleft='off',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')

panel2.xaxis.set_visible(False)

panel3.tick_params(axis='both',which='both',\
                   bottom='off', labelbottom='off',\
                   left='off', labelleft='off',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')

panel3.xaxis.set_visible(False)

plt.savefig('Gustincic_Mark_BME163_Undergrad_Assignment_Final.png',dpi=2000)

