import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.stats as stats
import matplotlib.image as mpimg
import matplotlib.patches as mplpatches

plt.style.use('BME163.mplstyle')

fig_1 = plt.figure(figsize=(5,3)) #width,height
panel1 = plt.axes([0.5/3 , 0.12 , 0.75/5 , 2.5/3],frameon=True) 
panel2 = plt.axes([.43 , 0.12 , 2.5/5 , 2.5/3],frameon=False)

panel1.set_xticks(np.arange(.5,9.5,1))

xLabels = [0,'',6,'',12,'',18,'']
panel1.set_xticklabels(xLabels)


R=np.linspace(255/255,56/255,101)
G=np.linspace(225/255,66/255,101)
B=np.linspace(40/255,157/255,101)


geneList = []
ppList = []

lineCount = 0

for line in open('BME163_Input_data3.txt'):
    a = line.strip().split('\t')
    if lineCount == 0:
        lineCount += 1
        continue
    CT0 = int(a[4])
    CT3 = int(a[5])
    CT6 = int(a[6])
    CT9 = int(a[7])
    CT12 = int(a[8])
    CT15 = int(a[9])
    CT18 = int(a[10])
    CT21 = int(a[11])
    peakPhase = float(a[13])
    ppList.append(peakPhase)
    geneInfo = [peakPhase,CT0,CT3,CT6,CT9,CT12,CT15,CT18,CT21]
    geneList.append(geneInfo)

bins = np.arange(0,26,2)
ppHist,bins1=np.histogram(ppList,bins)

#geneList contains 1262 elements

normalized = []
    
for gene in geneList:
    newGene = [gene[0]]
    count = 1
    for data in gene[1:]:
        normalized_data=((data-min(gene[1:]))/(max(gene[1:])-min(gene[1:])))*100
        newGene.append(int(normalized_data))
        count += 1
    normalized.append(newGene)

y_pos = 0
for gene in sorted(normalized, key = lambda x:x[0], reverse = True):
    x_pos = 0
    for data in gene[1:]:
        rectangle=mplpatches.Rectangle([x_pos,y_pos],1,1,facecolor=(R[data],G[data],B[data]),linewidth=0)
        panel1.add_patch(rectangle)
        x_pos +=1
    y_pos += 1



x_listNeg=[]
y_listNeg=[]

for x_value in np.linspace(-1.2,0,200):
    y_value = np.sqrt(1.44-x_value**2)
    x_listNeg.append(x_value*100)
    y_listNeg.append(y_value*100)
    x_listNeg.append(x_value*100)
    y_listNeg.append(-y_value*100)

panel2.plot(x_listNeg, y_listNeg, color = 'black', zorder = 250)

x_listPos=[]
y_listPos=[]

for x_value in np.linspace(-1,1,200):
    y_value = np.sqrt(1-x_value**2)
    x_listPos.append(x_value*100)
    y_listPos.append(y_value*100)
    x_listPos.append(x_value*100)
    y_listPos.append(-y_value*100)

panel2.plot(x_listPos, y_listPos, color = 'white',zorder=500)

circle1 = plt.Circle((0,0), 120, edgecolor='black', facecolor = 'white', fill=True, zorder = 250, lw = .8)
panel2.add_patch(circle1)

circle2 = plt.Circle((0,0), 100, edgecolor='black', facecolor = 'white', fill=False, zorder = 1000, lw = .8)
panel2.add_patch(circle2)


start = np.pi/2
end = np.pi/3
xListBlack = []
yListBlack = []
xCirc1list = []
yCirc1list = []
xCirc2list = []
yCirc2list = []
xCirc3list = []
yCirc3list = []

for step in ppHist:
    #print(step)
    height = step
    xList = []
    yList = []
    for radian in np.linspace(start,end,100):
        for radius in np.linspace(120,height+120,height):
            if radius == height + 120:
                x_pos=np.cos(radian)*radius
                xListBlack.append(x_pos)
                y_pos=np.sin(radian)*radius
                yListBlack.append(y_pos)

            x_pos=np.cos(radian)*radius
            xList.append(x_pos)
            y_pos=np.sin(radian)*radius
            yList.append(y_pos)
                
        xCirc1 = np.cos(radian)*220
        xCirc1list.append(xCirc1)
        yCirc1 = np.sin(radian)*220
        yCirc1list.append(yCirc1)

        xCirc2 = np.cos(radian)*320
        xCirc2list.append(xCirc2)
        yCirc2 = np.sin(radian)*320
        yCirc2list.append(yCirc2)

        xCirc3 = np.cos(radian)*420
        xCirc3list.append(xCirc3)
        yCirc3 = np.sin(radian)*420
        yCirc3list.append(yCirc3)



    panel2.plot(xList,yList,marker='o', mfc='black', mew=0, markersize=0,linewidth=0.7,color='grey')
    panel2.plot(xListBlack,yListBlack,marker='o', mfc='black', mew=0, markersize=0,linewidth=0.5,color='black')
    

    start -= (np.pi)/6
    end -= (np.pi)/6


panel2.plot([120,231],[0,0], marker='o', mfc='black', mew=0, markersize=0,linewidth=0.7,color='black')
panel2.plot([-120,-187],[0,0], marker='o', mfc='black', mew=0, markersize=0,linewidth=0.7,color='black')
panel2.plot([0,0],[120,276], marker='o', mfc='black', mew=0, markersize=0,linewidth=0.7,color='black')
panel2.plot([0,0],[-120,-336], marker='o', mfc='black', mew=0, markersize=0,linewidth=0.7,color='black')
panel2.plot([0,-137],[0,237], marker='o', mfc='black', mew=0, markersize=0,linewidth=0.7,color='black')    
panel2.plot([0,-154.15],[0,89], marker='o', mfc='black', mew=0, markersize=0,linewidth=0.7,color='black')
panel2.plot([0,-176.67],[0,-102], marker='o', mfc='black', mew=0, markersize=0,linewidth=0.7,color='black')
panel2.plot([0,-109.5],[0,-189.66], marker='o', mfc='black', mew=0, markersize=0,linewidth=0.7,color='black')
panel2.plot([0,110],[0,190.525], marker='o', mfc='black', mew=0, markersize=0,linewidth=0.7,color='black')
panel2.plot([0,156.75],[0,90.5], marker='o', mfc='black', mew=0, markersize=0,linewidth=0.7,color='black')
panel2.plot([0,217.37],[0,-125.5], marker='o', mfc='black', mew=0, markersize=0,linewidth=0.7,color='black')
panel2.plot([0,121.5],[0,-210.44], marker='o', mfc='black', mew=0, markersize=0,linewidth=0.7,color='black')

panel2.plot(xCirc1list,yCirc1list,marker='o', mfc='black', dashes = (1,1,2,1), mew=0, markersize=0,linewidth=0.3,color='black',zorder=3000)
panel2.plot(xCirc2list,yCirc2list,marker='o', mfc='black', mew=0, dashes = (1,1,2,1), markersize=0,linewidth=0.3,color='black')
panel2.plot(xCirc3list,yCirc3list,marker='o', mfc='black', mew=0, dashes = (1,1,2,1), markersize=0,linewidth=0.3,color='black')



panel2.text(-16,-7.5,r'CT',zorder = 1500, fontsize = 6)
panel2.text(-8,66,r'0',zorder = 1500, fontsize = 6)
panel2.text(-16,-80,r'12',zorder = 1500, fontsize = 6)
panel2.text(47.5,35,r'4',zorder = 1500, fontsize = 6)
panel2.text(47.5,-45,r'8',zorder = 1500, fontsize = 6)
panel2.text(-67.5,-45,r'16',zorder = 1500, fontsize = 6)
panel2.text(-67.5,35,r'20',zorder = 1500, fontsize = 6)

panel2.text(-275,-7.5,r'100',zorder = 1500, fontsize = 6)
panel2.text(-375,-7.5,r'200',zorder = 1500, fontsize = 6)
panel2.text(-475,-7.5,r'300',zorder = 1500, fontsize = 6)


panel1.set_xlim(0,8)
panel1.set_ylim(0,len(geneList))

panel2.set_ylim(-450,450)
panel2.set_xlim(-450,450)

panel1.set_xlabel('CT', fontsize=10)
panel1.set_ylabel('Number of genes', fontsize=10)


panel1.tick_params(axis='both',which='both',\
                   bottom='on', labelbottom='on',\
                   left='on', labelleft='on',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')

panel2.tick_params(axis='both',which='both',\
                   bottom='off', labelbottom='off',\
                   left='off', labelleft='off',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')



plt.savefig('Gustincic_Mark_BME163_Undergrad_Assignment_Week6.png',dpi=600)
