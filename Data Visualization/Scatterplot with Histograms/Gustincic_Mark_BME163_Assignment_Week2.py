import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np

plt.style.use('BME163.mplstyle') #if error occurs, remove '.mplstyle' from file name

plt.figure(figsize=(2,2))

#panel=plt.axes([left, bottom, width, height])
panel1=plt.axes([0.165,0.2,0.125,0.5]) #left
panel2=plt.axes([0.33,0.2,0.5,0.5]) #right
panel3=plt.axes([0.33,0.74,0.5,0.125]) #top

x_list=[]
y_list=[]
gene_length_list=[]

for line in open('BME163_Input_data_1.txt'):
    split_line=line.strip().split('\t')
    gene_name=split_line[0]
    gene_length=int(split_line[1])
    x_value=int(split_line[2])
    y_value=int(split_line[3])
    x_list.append(x_value)
    y_list.append(y_value)
    gene_length_list.append(gene_length)

#turn lists into arrays in order to be able to do math on the lists
x_array=np.array(x_list)
y_array=np.array(y_list)

x_array_log=np.log2(x_array+1)
y_array_log=np.log2(y_array+1)

alphaVals = np.linspace(.3,.38,len(gene_length_list))
count = len(gene_length_list)
rgbColors = np.zeros((count,4))
rgbColors[:,3] = alphaVals

panel2.scatter(x_array_log,y_array_log, \
               s=2,\
               color = rgbColors,\
               linewidth=0.0)

bins=np.arange(0,15,.5)

x_array_log_hist,bins1=np.histogram(x_array_log,bins)
y_array_log_hist,bins1=np.histogram(y_array_log,bins)

left = 0

for step in np.arange(0,len(x_array_log_hist),1):
    bottom = 0
    width=bins[step+1]-bins[step]
    height=x_array_log_hist[step]
    
    rectangle1=mplpatches.Rectangle((left,bottom),1,height,\
                                    linewidth=0.1,\
                                    facecolor=(0.5,0.5,0.5),\
                                    edgecolor=(0,0,0))

    panel3.add_patch(rectangle1)

    left += 1
    
for step in np.arange(0,len(y_array_log_hist),1):   
    left1=0
    bottom1=bins[step]
    height1=bins[step+1]-bins[step]
    width1=y_array_log_hist[step]

    rectangle2=mplpatches.Rectangle((left1,bottom1),width1,height1,\
                                linewidth=0.1,\
                                facecolor=(0.5,0.5,0.5),\
                                edgecolor=(0,0,0))

    panel1.add_patch(rectangle2)

panel1.set_xlim([400,0])
panel1.set_ylim([0,14])

panel2.set_xlim([0,14])
panel2.set_ylim([0, max(y_array_log)*1.3])

panel3.set_xlim([0,(len(x_array_log_hist))])
panel3.set_ylim([0,400])

panel1.set_xticks(np.arange(0, 800, 400))
panel1.set_yticks(np.arange(0, 16, 2))

panel3.set_yticks(np.arange(0, 800, 400))


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

panel3.tick_params(axis='both',which='both',\
                   bottom='off', labelbottom='off',\
                   left='on', labelleft='on',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')


plt.savefig('Gustincic_Mark_BME163_Assignment_Week2.png')


