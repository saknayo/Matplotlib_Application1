#rs_mapping2.py
import numpy
import matplotlib.pyplot as plt
import matplotlib
import os
plt.style.use('ggplot')
N=20

class RowData:
    def __init__(self,fp):
        self.fp=fp
        self.data_import()
        self.x=numpy.arange( self.rowdata.shape[1] )

    def data_import(self):
        tlist=[]
        with open(self.fp,'r',encoding='utf-8')as f:
            for line in f:
                tlist.append(eval(line.strip()))
        self.rowdata=numpy.array(tlist)
    def getX(self):
        return self.x
    def getReptag(self):
        return numpy.array(list( self.rowdata[:,i,0].mean() for i in range(N) ))
    def getRep(self):
        return numpy.array(list( self.rowdata[:,i,1].mean() for i in range(N) ))
    def getOther(self):
        return numpy.array(list( self.rowdata[:,i,2].mean() for i in range(N) ))

def autolabel(rects,ax):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height),
                ha='center', va='bottom')       
def meanPlot(ax,x,y):
    #print(x,y)
    width=1
    mbar=ax.bar(x+width,y,width=width)#,width=10*self.length/self.dots)
    ax.set_xticks(x+1.5*width)
    ax.set_xticklabels(x+1)
    #ax1.set_xlabel('row',size=11)
    ax.set_ylabel('Numbers(avg)',size=11)
    autolabel(mbar,ax)

data=RowData('row_sta.txt')
fig=plt.figure(figsize=[15,10])
fig.subplots_adjust(hspace=0.3, wspace=0.3)

ax1=fig.add_subplot(3,1,1)
meanPlot(ax1,data.getX(),data.getReptag())
ax1.set_xlabel('Reptag+c',size=11)
ax2=fig.add_subplot(3,1,2)
meanPlot(ax2,data.getX(),data.getRep())
ax2.set_xlabel('Rep',size=11)
ax3=fig.add_subplot(3,1,3)
meanPlot(ax3,data.getX(),data.getOther())
ax3.set_xlabel('Other',size=11)

plt.savefig('row_{}.png'.format( os.path.basename(os.getcwd()) ) ,dpi=300)
plt.show()