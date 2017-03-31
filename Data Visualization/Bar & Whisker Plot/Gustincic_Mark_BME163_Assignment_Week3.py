import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np

plt.style.use('BME163.mplstyle')#if throws an error, remove '.mplstyle' from file name

plt.figure(figsize=(5,2))

#bins are a list. create dictionary w/ key as bin value, value is y/x ratio

panel1=plt.axes([0.1,0.3,4/5,1/2])

x_list=[]
y_list=[]
gene_length_list=[]
lengthRatioList = []
newLengthList = []
uniqueLengthList =[]
ratioList = []
normalizedLength = []

for line in open('BME163_Input_data_1.txt'):
    split_line=line.strip().split('\t')
    gene_name=split_line[0]
    gene_length=int(split_line[1])
    x_value=int(split_line[2])
    y_value=int(split_line[3])
    if x_value and y_value > 0:
        lengthRatio = [gene_length, y_value/x_value]
        gene_length_list.append(gene_length) #contains unaltered lengths
        lengthRatioList.append(lengthRatio)

for length in gene_length_list:
    newLength = int(length/300)*300
    newLengthList.append(newLength)
dirtyList = list(set(newLengthList))

bins = np.arange(0,6300,300)
for length in dirtyList:
    if length in bins:
        uniqueLengthList.append(length)
        
binDict = {key:[] for key in uniqueLengthList}
       
    

for values in lengthRatioList: #for each list within lengthRatioList
    values[0] = int(values[0]/300)*300 #lengthRatioList now contains lists that contain normalized gene lengths & y/x ratio
    if values[0] in binDict.keys(): #if normalized gene length is between 0 and 6300
        binDict[values[0]].append(values[1])

#binDict now contains normalized gene length: list of ratios as key:value pair
        

for key,value in binDict.items():
    normalizedLength.append(key+300)
    ratioList.append(value)


bp=panel1.boxplot(ratioList, positions= normalizedLength, patch_artist=True, widths=125)

for box in bp['boxes']:
    box.set(edgecolor='black',facecolor='white',linewidth=.8)
for whisker in bp['whiskers']:
    whisker.set(color='black', linestyle='-',linewidth=.8)
for median in bp['medians']:
    median.set(color='black', linestyle='-',linewidth=.8)
for flier in bp['fliers']:
    flier.set(markersize=0)

labelList = []
for step in np.arange(0,6000,600):
    myString = ('{lowerLimit}\n to \n {upperLimit}'.format(lowerLimit = step, upperLimit = step+300))
    labelList.append(myString)

panel1.set_xticks(np.arange(300,6600,600))
panel1.set_xticklabels(labelList)

panel1.set_xlabel('Gene length')
panel1.set_ylabel('Expression level ratio')

panel1.tick_params(axis='both',which='both',\
                   bottom='on', labelbottom='on',\
                   left='on', labelleft='on',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')




panel1.set_xlim(0,6000)
panel1.set_ylim(0,3)


#plt.show()
plt.savefig('Gustincic_Mark_BME163_Assignment_Week3.pdf')
