import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np

plt.style.use('BME163.mplstyle') #if error occurs, remove '.mplstyle' from file name

plt.figure(figsize=(3.42,2))

panel1=plt.axes([0.1,0.2,1/3.42,1/2])
panel2=plt.axes([0.55,0.2,1/3.42,1/2])


panel1.set_xlim([0.0, 1.0])
panel1.set_ylim([0.1, 1.0])

panel1.set_xticks(np.arange(0.0, 1.1, 0.2))
panel1.set_yticks(np.arange(0.1, 1.1, 0.1))

panel1.tick_params(direction='out',\
                   axis='both',which='both',\
                   bottom='on', labelbottom='on',\
                   left='on', labelleft='on', \
                   right='off', labelright='off',\
                   top='off', labeltop='off')


panel2.set_xticks(np.arange(0.1, 1.0, 0.1))
panel2.set_yticks(np.arange(0.1, 1.0, 0.1))

x = np.arange(0, 1, .01)
y = np.sqrt(1-x**2)

for p in range(100):
    circ = mplpatches.Circle((x[p],y[p]),0.013,
                              linewidth=0,\
                              linestyle='-',\
                              color=(y[p],y[p],y[p]))
    panel1.add_patch(circ)

panel2.tick_params(direction='in',\
                   axis='both',which='both',\
                   bottom='off', labelbottom='off',\
                   left='off', labelleft='off', \
                   right='off', labelright='off',\
                   top='off', labeltop='off')

for y_pos in range(10):
    for x_pos in range(10):
        rectangle1=mplpatches.Rectangle((x_pos/10,y_pos/10),1/10,1/10,\
                                        facecolor=(x_pos/10, y_pos/10, 1),\
                                        linewidth=1,\
                                        edgecolor = 'black')
        panel2.add_patch(rectangle1)

plt.savefig('Gustincic_Mark_BME163_Assignment_Week1.png')


