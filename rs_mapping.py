import matplotlib.pyplot as plt
import matplotlib
import numpy 
import re
from functools import lru_cache
from collections import defaultdict
import time
plt.style.use('ggplot')

_debug=False

class DataManager:
    def __init__(self,fp,incell=True,outcell=True):
        self._id=time.time()
        self.records=self.getRecordsList(fp)
        self.record_num=len(self.records)

        self.step_p=re.compile('=>(?P<step>\d+)')
        self.steps=numpy.zeros(self.record_num,dtype=int)

        self.incell=incell
        self.outcell=outcell

        self.rbz_in={}
        self.rbz_out={}
        self.others={}

        self.rbz_in_patterns={}
        self.rbz_out_patterns={}
        self.others_patterns={}

        self.rep_str='Rep'
        self.ar_str='Ar'
        self.nr_str='Nr'
    def setRep(self,s):
        self.rep_str=s
    def setAr(self,s):
        self.ar_str=s
    def setNr(self,s):
        self.nr_str=s
    def getId(self):
        return self._id

    def getRecordNum(self):
        return self.record_num

    def getX(self,plotint):
        return self.steps[::plotint] 

    def getRecord(self,tag,_type,plotint):
        if tag in self.others :
            return self.others[tag][::plotint]
        elif _type == 'incell' :
            try :
                return self.rbz_in[tag][::plotint]
            except KeyError :
                print(" {} is not prased in db{}!!!".format(tag,self.getId()))
                raise
        elif _type == 'outcell' :
            try :
                return self.rbz_out[tag][::plotint]
            except KeyError :
                print(" {} is not prased in db{}!!!".format(tag,self.getId()))
                raise

    def getRecordsList(self,fp):
        with open(fp,'r',encoding='utf-8')as f:
            return self.splitRecordsByStep( f.read() )
    def splitRecordsByStep(self,filecontent):
        return re.split('Step',filecontent)[1:]

    def rbzDBGen(self,tag_str):
        l=['{}',  '{}_s',  '{}_d',  '{}tagc_s',  '{}tag',  '{}tag_s',  '{}tag_d',  '{}tagc',  '{}tagc_d']
        for i in l :
            self.register( i.format(tag_str), 'incell'  )
            self.register( i.format(tag_str), 'outcell' )

    def rbzInit(self):
        if self.rep_str :
            self.rbzDBGen(self.rep_str)
        if self.ar_str :
            self.rbzDBGen(self.ar_str)
        if self.nr_str :
            self.rbzDBGen(self.nr_str)
        for i in ['Tag', 'Tag_s', 'Tag_d', 'DTag', 'DTag_s', 'DTag_d'] :
            self.register( i   , 'incell'  )
        for i in ['Tag', 'Tag_s', 'Tag_d', 'DTag', 'DTag_s', 'DTag_d'] :
            self.register( i   , 'outcell' )

    def othersInit(self):
        for i in [ 'RNA',   'Prna', 'AM',   'Pam',  'Chain_3', 'Cell', 'C_Rep', 'C_Nr', 'C_Ar', 'C_Re_Nr',  'C_Re_Ar',  'C_Nr_Ar',  'C_Re_Nr_Ar' ,'Health_Cell' ] :
            self.register( i,'other' )


    def addArray(self,tag,adic,_dtype=int):
        adic.update({ tag : numpy.zeros(self.record_num,dtype=_dtype) })

    def patternA(self,_type,s):
        ''' default pattern generate based on tag '''
        if   _type == 'incell' :
            return re.compile( '{}:[(](?P<{}>\d+)'.format(s,s) )
        elif _type == 'outcell' :
            return re.compile( '{}:[(]0  [)](?P<{}>\d+)'.format(s,s) ) 
        elif _type == 'other' :
            return re.compile( '{}:(?P<{}>\d+)'.format(s,s) )

    def register(self,tag,_type,pattern=None,dtype=int):
        ''' tag=name, _type=incell|outcell|other , pattern=re.compile(), dtype=int|float '''
        if pattern == None :
            pattern=self.patternA(_type,tag)

        if _type == 'incell' :
            self.rbz_in_patterns[tag]=pattern
            if tag not in self.rbz_in :
                self.addArray(tag, self.rbz_in, _dtype=dtype)
        elif _type == 'outcell' :
            self.rbz_out_patterns[tag]=pattern
            if tag not in self.rbz_out :
                self.addArray(tag, self.rbz_out, _dtype=dtype)
        else :
            self.others_patterns[tag]=pattern
            if tag not in self.others :
                self.addArray(tag, self.others, _dtype=dtype)
        if _debug :
            print( "#! {} registed in {}.".format(tag,_type) )

    def defaultRegister(self):
        ''' register all enzyme(incell&outcell) and other info '''
        self.rbzInit()
        self.othersInit()

    def recordPraser(self,s_record,n): 
        self.steps[n]=eval( self.step_p.search(s_record).group('step') )
        if self.incell :
            for i in self.rbz_in_patterns :
                self.rbz_in[i][n]    = eval( self.rbz_in_patterns[i].search(s_record).group(i) )
        if self.outcell :
            for i in self.rbz_out_patterns :
                self.rbz_out[i][n]   = eval( self.rbz_out_patterns[i].search(s_record).group(i) )
        for i in self.others_patterns :
            self.others[i][n]        = eval( self.others_patterns[i].search(s_record).group(i) )

    def praseAll(self):
        n=0
        for r in self.records :
            self.recordPraser(r,n)
            n+=1

    def __eq__(self,ba):
        return self.getId() == ba

class DrawManager:
    def __init__(self,maindb,plotint=1,dpi=600):
        self.main_db=maindb
        self.sub_db={}
        self.plotint=plotint
        self.dpi=dpi
        self.dots=int( self.main_db.getRecordNum()/self.plotint )
        self.length= ( 2*self.main_db.getRecordNum())/self.dpi
        self.to_draw=dict() # colum : {tag:[db1,db2]}
        self.marks=dict()
        self.standard_colors={  
              'Tn':'c' ,          'RNA':'c' ,         'Prna':'c' ,           'Ta':'c' ,           'AM':'c' ,         'Pam':'c' ,      'Chain_3':'r' ,    
            'Cell':'m' ,        'C_Rep':'m' ,         'C_Nr':'m' ,         'C_Ar':'m' ,      'C_Re_Nr':'m' ,     'C_Re_Ar':'m' ,      'C_Nr_Ar':'m' ,   'C_Re_Nr_Ar':'m' ,  'Health_Cell':'m' ,    
             'Tag':'k' ,        'Tag_s':'k' ,        'Tag_d':'k' ,         'DTag':'k' ,       'DTag_s':'k' ,      'DTag_d':'k' , 
             'Rep':'r' ,        'Rep_s':'r' ,        'Rep_d':'r' ,       'Reptag':'r' ,     'Reptag_s':'r' ,    'Reptag_d':'r' ,      'Reptagc':'r' ,    'Reptagc_s':'r' ,    'Reptagc_d':'r' , 
              'Nr':'g' ,         'Nr_s':'g' ,         'Nr_d':'g' ,        'Nrtag':'g' ,      'Nrtag_s':'g' ,     'Nrtag_d':'g' ,       'Nrtagc':'g' ,     'Nrtagc_s':'g' ,     'Nrtagc_d':'g' , 
              'Ar':'b' ,         'Ar_s':'b' ,         'Ar_d':'b' ,        'Artag':'b' ,      'Artag_s':'b' ,     'Artag_d':'b' ,       'Artagc':'b' ,     'Artagc_s':'b' ,     'Artagc_d':'b' ,}
        self.extra_colors=[
                    '#ADFF2F'     ,'#7FFF00'     ,'#7CFC00'     ,'#00FF00'     ,'#32CD32'     ,'#98FB98'     ,'#90EE90'     ,'#00FA9A'     ,'#00FF7F'     ,'#3CB371'     ,'#2E8B57'     ,
                    '#228B22'     ,'#008000'     ,'#006400'     ,'#9ACD32'     ,'#6B8E23'     ,'#556B2F'     ,'#66CDAA'     ,'#8FBC8F'     ,'#20B2AA'     ,'#008B8B'     ,'#008080'     ]
        self.standard_markers={  
              'Tn':'s' ,          'RNA':'x' ,         'Prna':'+' ,           'Ta':'o' ,           'AM':'^' ,         'Pam':'*' ,      'Chain_3':'D' ,    
            'Cell':'s' ,        'C_Rep':'x' ,         'C_Nr':'+' ,         'C_Ar':'o' ,      'C_Re_Nr':'^' ,     'C_Re_Ar':'*' ,      'C_Nr_Ar':'D' ,   'C_Re_Nr_Ar':'.' ,  'Health_Cell':',' ,    
             'Tag':'s' ,        'Tag_s':'x' ,        'Tag_d':'+' ,         'DTag':'o' ,       'DTag_s':'^' ,      'DTag_d':'*' , 
             'Rep':'s' ,        'Rep_s':'x' ,        'Rep_d':'+' ,       'Reptag':'o' ,     'Reptag_s':'^' ,    'Reptag_d':'*' ,      'Reptagc':'D' ,    'Reptagc_s':'.' ,    'Reptagc_d':',' , 
              'Nr':'s' ,         'Nr_s':'x' ,         'Nr_d':'+' ,        'Nrtag':'o' ,      'Nrtag_s':'^' ,     'Nrtag_d':'*' ,       'Nrtagc':'D' ,     'Nrtagc_s':'.' ,     'Nrtagc_d':',' , 
              'Ar':'s' ,         'Ar_s':'x' ,         'Ar_d':'+' ,        'Artag':'o' ,      'Artag_s':'^' ,     'Artag_d':'*' ,       'Artagc':'D' ,     'Artagc_s':'.' ,     'Artagc_d':',' ,}
        self.extra_markers=['1','2','3','4','8','h','H','x','p','d','<','>','v']

    def getMarker(self,tag):
        if tag in self.standard_markers :
            return self.standard_markers[tag]
        else :
            return self.extra_markers[hash(tag)%len(self.extra_markers)]
    def getFaceColor(self,tag):
        if tag in self.standard_colors  :
            return self.standard_colors[tag]
        else :
            return self.extra_colors[hash(tag)%len(self.extra_colors)]
    def getEdgeColor(self,tag,_type):
        if _type == 'outcell' :
            return self.getFaceColor(tag)
        return 'k'
    def getLinestyle(self,_db):
        if _db == self.main_db :
            return '-'
        else :
            return [':', '-.', None][hash(_db.getId())%3]


    def addPlot(self,tag,_type,colum,_db = None):
        ''' tag = Reptag , _type = incell|outcell|other , _db = instance of DataManager '''
        if colum not in self.to_draw :
            self.to_draw[colum]=defaultdict(list)
        if _db : # comparation view
            self.to_draw[colum][tag].append([_db,_type])
            #self.marks[str(colum)+tag+_type+str(_db.getId)]
        else :
            self.to_draw[colum][tag].append([self.main_db,_type])
    def defaultDraw(self):
        self.addPlot('Rep','outcell',1)
        self.addPlot('Reptag','outcell',1)
        self.addPlot('Reptagc','outcell',1)
        self.addPlot('Rep','outcell',2)
        self.addPlot('Reptag','outcell',2)
        self.addPlot('Reptagc','outcell',2)
        self.addPlot('Rep','outcell',3)
        self.addPlot('Reptag','outcell',4)
        self.addPlot('Reptagc','outcell',5)
        self.addPlot('RNA','other',5)

        self.draw()

    def outcellDraw(self):
        self.addPlot('Reptag','outcell',1)
        self.addPlot('Reptagc','outcell',1)
        self.addPlot('Chain_3','other',1)
        self.addPlot('Tag','outcell',2)
        self.addPlot('DTag','outcell',2)
        self.addPlot('RNA','other',3)
        self.addPlot('Prna','other',3)

        self.draw()

    def incellDraw(self):
        self.addPlot('Reptag','incell',1)
        self.addPlot('Artag','incell',1)
        self.addPlot('Nrtag','incell',1)

        self.addPlot('Tag','incell',2)
        self.addPlot('DTag','incell',2)

        self.addPlot('Rep','incell',3)
        self.addPlot('Reptag','incell',3)
        self.addPlot('Reptagc','incell',3)
        self.addPlot('Chain_3','other',3)

        self.addPlot('Ar','incell',4)
        self.addPlot('Artag','incell',4)
        self.addPlot('Artagc','incell',4)

        self.addPlot('Nr','incell',5)
        self.addPlot('Nrtag','incell',5)
        self.addPlot('Nrtagc','incell',5)

        self.addPlot('RNA','other',6)
        self.addPlot('Prna','other',6)

        self.draw()

    def draw(self):
        self.width=len(self.to_draw)*3
        fig=plt.figure(figsize=[self.length,self.width])
        for c in self.to_draw :
            ax=fig.add_subplot( max(self.to_draw), 1, c )
            for tag in self.to_draw[c] :
                for db in self.to_draw[c][tag] :
                    _type=db[1]
                    #print(tag,repr(db[0].getRecord(tag,_type,self.plotint)))
                    ax.plot(db[0].getX(self.plotint),
                            db[0].getRecord(tag,_type,self.plotint),
                            marker=self.getMarker(tag),
                            linestyle=self.getLinestyle(db[0]),
                            linewidth=0.5,
                            mfc =self.getFaceColor(tag) ,
                            mec =self.getEdgeColor(tag,_type),
                            mew=1,  # mew:markeredgewidth ,mec:markeredgecolor , mfc:markerfacecolor
                            markersize=1, 
                            label=tag,
                            )
                    plt.legend()
        plt.legend()
        plt.savefig('test_{}.png'.format(self.plotint),dpi=self.dpi)
        plt.show()





if __name__ == '__main__' :
    # outcell draw test
    '''testdraw=DataManager("Cell_Rep_Nsr.txt",incell=False)
    testdraw.setNr('')
    testdraw.setAr('')
    testdraw.defaultRegister()
    #testdraw.register('dsn','other',re.compile(">4:\d+,(?P<dsn>\d+),\d+"))
    testdraw.praseAll()'''

    testdraw=DataManager("Cell_Rep_Nsr.txt")
    testdraw.defaultRegister()
    testdraw.praseAll()

    drawManager=DrawManager(testdraw,plotint=10)
    #drawManager.addPlot('dsn','other',4)
    #drawManager.addPlot('Rep','outcell',3,testdraw2)
    #drawManager.defaultDraw()
    drawManager.incellDraw()










