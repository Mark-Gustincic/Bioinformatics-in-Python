import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import scipy.stats as stats
import math

plt.style.use('BME163.mplstyle')#if throws error, remove '.mplstyle'

plt.figure(figsize=(2,2))

panel1=plt.axes([0.25,0.18,.7,.7])

bins = np.arange(1,5,1)
thisDict = {key:[] for key in bins}
indexDict = {}
normExpList = []
neList = []
indexNameList = []
inList = []
countList = []




for line in open('BME163_Input_data2.psl'):
    split_line=line.strip().split('\t')
    tScript=split_line[13]
    #read = split_line[9]
    index=int(split_line[9][0])
    indexNameList = []
    if index not in indexDict.keys():
        indexDict[index] = {}
    elif tScript not in indexDict[index].keys():
        indexDict[index][tScript] = 1
    else:
        count = indexDict[index][tScript]
        count += 1
        indexDict[index][tScript] = count

for infoDict in indexDict.values():
    denominator = 0
    collectList = []
    for count in infoDict.values():
        denominator += count
    for count in infoDict.values():
        normExp = np.log2(((count/denominator)*10000)+1)
        collectList.append(normExp)
    neList.append(collectList)

p_value1 = stats.mannwhitneyu(neList[0], neList[1],alternative='two-sided')
p_value2 = stats.mannwhitneyu(neList[0], neList[2],alternative='two-sided')
p_value3 = stats.mannwhitneyu(neList[2], neList[3],alternative='two-sided')


plt.text(.9,12.1,'p={}'.format(round(p_value1[1],3)),fontsize=6)
plt.plot([.9, 2], [11.9, 11.9], color='black', linestyle='-', linewidth=.75)

plt.text(1.5,13.8,'p={}'.format(round(p_value2[1],3)),fontsize=6)
plt.plot([1, 3], [13.6, 13.6], color='black', linestyle='-', linewidth=.75)

plt.text(3,12.1,'p={}'.format(round(p_value3[1],3)),fontsize=6)
plt.plot([3, 4], [11.9, 11.9], color='black', linestyle='-', linewidth=.75)

bp=panel1.boxplot(neList, positions= [1,2,3,4], patch_artist=True, widths=.5)


for box in bp['boxes']:
    box.set(edgecolor='black',facecolor='white',linewidth=1)
for whisker in bp['whiskers']:
    whisker.set(color='black', linestyle='-',linewidth=1)
for median in bp['medians']:
    median.set(color='black', linestyle='-',linewidth=1)
for flier in bp['fliers']:
    flier.set(markersize=0)

step_size=0.01
cutoff=0.05

for num,data in enumerate(neList):
    center = num +1
    left_bound=center - .75
    right_bound=center + 1.25
    placed_points=[]
    for y_value in data:
        if len(placed_points)==0:
            placed_points.append((center, y_value))
        else:
            potential_x_position=[]
            for x_position in np.arange(left_bound,right_bound, step_size):
                distances=[]
                for placed_point in placed_points:
                    distance=((x_position-placed_point[0])**2+((y_value-placed_point[1])/3)**2)**0.5
                    distances.append(distance)
                if min(distances)>cutoff:
                    potential_x_position.append(x_position)
            if len(potential_x_position)>0:
                 best_x_position=sorted(potential_x_position,key=lambda x: np.absolute(x-center))[0]
                 placed_points.append((best_x_position,y_value)) 
            else:
                 print('point not placed: ',y_value)
    for point in placed_points:
        panel1.plot(point[0],point[1],marker='o',ms=1.5,mfc='black', mew=0,linewidth=0)



for box in bp['boxes']:
    box.set(edgecolor='black',facecolor='white',linewidth=.7)
for whisker in bp['whiskers']:
    whisker.set(color='black', linestyle='-',linewidth=.7)
for median in bp['medians']:
    median.set(color='black', linestyle='-',linewidth=.7)
for flier in bp['fliers']:
    flier.set(markersize=0)


panel1.set_xticks(bins)
panel1.set_yticks(np.arange(0,16,2))

panel1.set_xlabel('Library Index', fontsize=6)
panel1.set_ylabel('Normalized Expression (log$\mathregular{_2}$)', fontsize=6)

panel1.tick_params(axis='both',which='both',\
                   bottom='on', labelbottom='on',\
                   left='on', labelleft='on',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')


panel1.set_xlim(0,5)
panel1.set_ylim(0,15)

plt.savefig('Gustincic_Mark_BME163_Assignment_Week4.pdf')



        





    
    
