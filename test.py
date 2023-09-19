from shapely import geometry
import math
import numpy as np
import pandas as pd
from ttim import *
from pylab import *
import matplotlib
from matplotlib import cm
import fiona
import os,json
from descartes.patch import PolygonPatch
import shutil
import PySimpleGUI as sg
import ctypes
import operator

pd.set_option('display.max_columns', None) 
pd.options.display.width=None
pd.set_option('display.max_rows', 5) 
pd.set_option('display.expand_frame_repr', True)
pd.options.display.float_format = '{:,.2f}'.format

params = {'font.family': 'sans-serif',
          'font.sans-serif': 'arial',
          'axes.labelsize': 10,
          'axes.facecolor': '#ffffff', 
          'axes.labelcolor': 'black',
          'legend.fontsize': 8,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'lines.linewidth': 1,
          'grid.color': 'grey',
          'grid.linestyle': 'dashed',
          'grid.linewidth': 0.5,
          'text.usetex': False,
          'font.style': 'normal',
          'font.variant':'normal',
          'figure.facecolor': 'white',
          'font.size':10,
          'figure.autolayout': True,
          'figure.figsize': (10,8),
          'figure.dpi': 100,
          }
plt.rcParams.update(params)

class GEF:
    def __init__(self):
        self._data_seperator = ' '
        self._columns = {}
        self.x = 0.
        self.y = 0.
        self.z = 0.
        self.dz = []
        self.qc = []
        self.pw = []
        self.wg = []
        self.c  = []
        self.kb = []
        self.dist =[]
        self.npor=[]
        
    def readFile(self, filename):
        lines = open(filename, 'r').readlines()
        for line in lines:
            reading_header = True
        for line in lines:   
            if reading_header:
                self._parseHeaderLine(line)
            else:
                self._parseDataLine(line)
            if line.find('#EOH') > -1:
                if self._check_header():
                    reading_header = False
                else:
                    print(filename,'bestaat al')
                    return
            
    def _check_header(self):
        if not 1 in self._columns:
            return False
        if not 2 in self._columns:
            return False
        return True

    def _parseHeaderLine(self, line):
        for xe in ['#COMMENT', 'Peil=', 'uitvoerder', 'materieel','WATERSTAND',
                    'opmerkingen#MEASUREMENTTEXT','==','= NAP','antropogeen']:
            if xe in line:
                return          
        if len(line.split()) == 0:
            return
        
        keyword, argline = line.split('=')         
        keyword = keyword.strip()
        argline = argline.strip()
        args = argline.split(',')
      
        if '#XYID' in line:
            argline = argline.replace('.','')        
        args = argline.split(',')

        if 'Waterspanning'  in line:
            self.u = float(args[3])
        if 'Waterdruk'  in line:
            self.u = float(args[3]) 
        try:
            if 'waterspanning'  in line:
                self.u = int(args[3])
        except ValueError:
            return

        if keyword=='#XYID':
            if float(args[1]) < 1e5:
                args[1] = args[1]
            else:
                args[1]=args[1].replace('.','')
            self.x = float(args[1])
            args[2]=args[2].replace('.','')
            self.y = float(args[2])
            if (len(str(int(self.x))))>5:
                 self.x=int(self.x/pow(10,len(str(int(self.x)))-6))
            if (len(str(int(self.y))))>5:
                 self.y=int(self.y/pow(10,len(str(int(self.y)))-6))
            if self.x > 3e5:
                self.x=self.x/10

        elif keyword=='#ZID':
            self.z = round(float(args[1]),3)
           
        elif keyword=='#COLUMNINFO':
            column = int(args[0])
            dtype = int(args[-1])
            if dtype==11:
                dtype = 10
            self._columns[dtype] = column - 1    
       
    def _parseDataLine(self, line):
        line=line.strip()
        line = line.replace('|',' ')
        line = line.replace(';',' ')
        line = line.replace('!',' ')
        line = line.replace(':',' ')
        args=line.split()
        for n, i in enumerate(args):
            if i in ['9.9990e+003','-9999.9900','-1.0000e+007','-99999','-99999.0',
                  '-99999.00','-99999.000','-9.9990e+003','999.000', '-9999.99', '999.999',
                  '99','9.999','99.999', '999.9']:
                args[n] = '0.1'       
        if len(line.split()) == 0:
            return
 
        zz  = round(abs(float(args[self._columns[1]])),4)
        dz = round(self.z - zz,4)
        qc = round(float(args[self._columns[2]]),4)  

        slope    =  0.0104
        intercept=  0.0190
                         
        try:
            pw = float(args[self._columns[3]]) 
            if pw<-10:
                pw=0.1
            elif pw>5:
                pw=slope*qc+intercept             
                
        except KeyError: #Deze regel maakt van een tweekolommer een driekolommer
            pw=slope*qc+intercept

        self.dz.append(dz)
        self.qc.append(qc)
        self.pw.append(pw)
        if qc<=0.001:
            qc=0.1
            self.wg.append(10.)
        else:
            wg = abs((pw / qc) * 100.)
        wg = abs((pw / qc) * 100.)

##############################################K-waarde        
        if wg>=0.0:
            if wg >5: wg=15
            ke=math.exp(wg)
        if ke <=0:  ke=1
        else:
            kb  = (qc / ke)*2
            self.kb.append(kb)
        
        if kb <=0.1:
            npor=0.05
        elif 0.1<kb<=1:
            npor=0.1
        elif 1<kb<=30:
            npor=0.35
        elif kb>30:
            npor=0.25
        self.npor.append(npor)
            
    def asNumpy(self):
        return np.transpose(np.array([self.dz, self.npor, self.kb]))

    def asDataFrame(self):
        a = self.asNumpy()
        return pd.DataFrame(data=a, columns=['depth','npor', 'k'])
        
    def plt(self, filename):
        df = self.asDataFrame()
        df = df.sort_values('depth', ascending=False)

        if df.empty:
            return df
        
        df = df.rolling(50).mean() 
        df = df.iloc[:: 50]
        df=df.dropna()
        df = df.reset_index(drop=True)
        
        df['kzkh']= 2/(33.653*np.exp(-0.066*(df['k'])))
        df['kzkh']= np.where(df['kzkh']>1, 1, df['kzkh'])

        dzend=-35
        print('Aangenomen diepte aq = ', dzend, '[m NAP]')
        dfd=df['depth']
        dfd.loc['ld']= dzend
        dfn = df.iloc[: , [1]].copy()   
        dfn = pd.concat([df, pd.DataFrame.from_records([{'depth':dzend},])], ignore_index=True)
        dfnp= pd.concat([df, pd.DataFrame.from_records([{'npor': 0.25},])],  ignore_index=True)

        hstar = 7.05
        print('hstar = ', hstar, '  [m NAP]')
        
        df['k'] =df['k']*0.475
        
        print('KD = ',int(df['k'].sum()+(abs(dzend)-abs(df.iloc[-1,0]))*df.iloc[-1,2]),'  [m2/dag]')
        
        """
        Op basis van de vorm moet er een wat hogere S zijn dan bij de pompproef.
        df['S'] = np.where(df['depth']>  6, 0.2, 1e-5)   
        df['S'] = np.where(df['depth']<  1, 5e-6, df['S'])   
        df['S'] = np.where(df['depth']< -9, 5e-6, df['S'])    
        df['S'] = np.where(df['depth']<-14, 1e-5, df['S'])      
        Verder zou bij 7 de k waarde wat lager moeten zijn (meer verlaging immers)
        Deze dus teruggezet naar .875 ipv .95.
        11/9 Het lijkt eerder of de laag = 1 meter wat hoger middelt dan de uitwerking
        van de PP. In beide gevallen ligt de KD rond 830 [m2/dag] (UITV) / 818 (PP)
        Maar, als d efiletrs inderdaad 10 meter lang waren moet de KD beperkt zijn tot 
        451 [m2/dag], meer conform de sondeerinterpretaties. PP doorrekenen mbv deze waarden.
        """
        df['S'] = np.where(df['depth']>  6, 0.2, 0.025)   
        df['S'] = np.where(df['depth']<  1, 5e-3, df['S'])   
        df['S'] = np.where(df['depth']< -9, 5e-6, df['S'])    
        df['S'] = np.where(df['depth']<-14, 1e-5, df['S'])      
        
        df.to_excel(filename+'_data.xlsx')

        ml=Model3D(kaq=df['k'], z=dfn[ 'depth'], Saq= df['S'], kzoverkh=df['kzkh'], tmin=1e-6, 
                   tmax=200, M=10)

        mat=pd.read_excel('G:/Dick/Digitale boringen/50/Project Rijen/Controle bemaling/loc_II_Voor_u.xlsx', engine='openpyxl')  
        mat=mat.replace(',','.')
        
        well    = []
        ret     = []
        pxy     = []
        damwand = []
        wellad  = []
        wellvs  = []
        
        well      = mat[mat['Naam'].str.match ('w')]        
        ret       = mat[mat['Naam'].str.match ('r')]        
        pxy       = mat[mat['Naam'].str.match ('p')]   
        damwand   = mat[mat['Naam'].str.match ('dw')]  
        damwand2  = mat[mat['Naam'].str.match ('dx')]  

        x1=pxy.iloc[0,1]
        y1=pxy.iloc[0,2]
        
        ##Onttrekkingsfilter
        tf =  0
        bf = -5
        # nwells= len(well)
        # nret  = len(ret)
        # print('Aantal bronnen ', nwells)
        # print('Aantal retourbronnen ', nret)
        
        tfind = df['depth'].sub(tf).abs().values.argmin()
        bfind = df['depth'].sub(bf).abs().values.argmin() 
        well1 = arange(tfind, bfind, 1)

        ## Retourfilters
        vstf = -10
        vsbf = -12
        vstfind = df['depth'].sub(vstf).abs().values.argmin()
        vsbfind = df['depth'].sub(vsbf).abs().values.argmin() 
        well2 = arange(vstfind, vsbfind, 1)
        
        debreg=pd.read_excel('G:/Dick/Digitale boringen/50/Project Rijen/Controle bemaling/Debieten_Py_u.xlsx', index_col=None, skiprows=0, engine='openpyxl')  

        debtijd = list(zip(debreg['dagnr']-1, debreg['totaal']/19*24)) #dagnr-1 omdat de registratie 1 dag achterloopt op de werkelijk aanzettijd pomp
        debret  = list(zip(debreg['dagnr']-1, debreg['retour']/45*-24)) 
        

        for rnum in well.index:
            punt = well['Naam'][rnum]
            x    = well['x'][rnum] 
            y    = well['y'][rnum]
            punt = Well(ml, x,y, tsandQ=debtijd,  rw=0.225, layers=well1)
 
        for rnum in ret.index:
            puntr = ret['Naam'][rnum]
            x     = ret['x'][rnum] 
            y     = ret['y'][rnum]
            fact  = ret['rd'][rnum]
            # voor  =  [ x for x in debret]
            dbrt = [(t[0], t[1]*fact) for t in debret]
            puntr = Well(ml, x,y, tsandQ=dbrt, rw=0.125, layers=well2)

####################################
        ml.solve()
# ###################################
        pb= [1,2,3,4,5,6]


        M1d = 1.2
        M2d = 1.5
        M3d = 1.4
        M4d = 1.1
        M5d = 1.5
        M6d = 1.0
        M8d = 3.5

        M4 = df['depth'].sub(M4d).abs().values.argmin()        
        M1 = df['depth'].sub(M1d).abs().values.argmin()        
        M2 = df['depth'].sub(M2d).abs().values.argmin()        
        M3 = df['depth'].sub(M3d).abs().values.argmin()        
        M5 = df['depth'].sub(M5d).abs().values.argmin()     
        M6 = df['depth'].sub(M6d).abs().values.argmin()        
        M8 = df['depth'].sub(M8d).abs().values.argmin()        

        for p_b in pb: #laagnummers: zie debieten_py xlsx
            t=np.linspace(0,40,480)

            s1=ml.head(25,2,t)   #PB_M1 start op 18-7 op +6.98
            plt.plot(t,s1[M1]+6.64, lw=3, alpha=0.5, color='purple')

            s2=ml.head(38,-2,t)   #PB_M2 start op 18-7 op +6.98
            plt.plot(t,s2[M2]+6.64, lw=3, alpha=0.5, color='blue')

            s3=ml.head(54,-5, t) #PB_M3 start op 18.7 op +7.10
            plt.plot(t,s3[M3]+6.63, lw=3, alpha=0.5, color='grey')

            s4=ml.head(91,-15,t)   #PB_M4 start op 18-7 op +7.19 hierop zijn de berkende lijnen vastgezet verschil = retour oid
            plt.plot(t,s4[M4]+6.71, lw=3, alpha=0.5, color='red')

            s5=ml.head(2,55, t) #PB_M5 start op 18.7 op +7.02
            plt.plot(t,s5[M5]+6.96, lw=3, alpha=0.5, color='orange')
           
            s6=ml.head(-56,53, t) #PB_M5 start op 18.7 op +7.06
            plt.plot(t,s6[M6]+6.92, lw=3, alpha=0.5, color='black')
 
            s8=ml.head(-34,-36, t) #PB_M8 start op 18.7 op +7.10
            plt.plot(t,s8[M8]+6.91, lw=3, alpha=0.5, color='green')
         
        pbuitv=pd.read_excel('G:/Dick/Digitale boringen/50/Project Rijen/Controle bemaling/M4_py_u.xlsx', index_col=None, skiprows=0, engine='openpyxl')  
        x  = pbuitv['dagnr']
        y1 = pbuitv['pb_M1_c']
        y1 = np.where(y1<1.7,np.nan, y1)
        y2 = pbuitv['pb_M2_c']
        y3 = pbuitv['pb_M3_c']
        y4 = pbuitv['pb_M4_c']
        y5 = pbuitv['pb_M5_c']
        y6 = pbuitv['pb_M6_c']
        y8 = pbuitv['pb_M8_c']
        
        plt.plot(x,y1, '+', color='purple', markersize=2)
        plt.plot(x,y2, '+', color='blue',   markersize=2)
        plt.plot(x,y3, '+', color='grey',   markersize=2)
        plt.plot(x,y4, '+', color='red',    markersize=2)
        plt.plot(x,y5, '+', color='orange', markersize=2)
        plt.plot(x,y6, '+', color='black',  markersize=2)
        plt.plot(x,y8, '+', color='green',  markersize=2)
        
        plt.legend(['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M8'])
        plt.grid(axis='both')
        plt.xlim(0,40)
        plt.xticks(np.arange(0,40,2))
        plt.ylim(1,8)
        plt.xlabel('tijd [dagen]    0 = 3 aug : 20 = 23 aug : 40 = 12 sept')
        plt.ylabel('GWS [m NAP]')
        plt.title(filename)

        plt.savefig(filename + '_u_Tijd_dh.png', bbox_inches='tight')
        plt.show()
        plt.close()
        
for filename in os.listdir(os.getcwd()):
    if filename.endswith ('.GEF') or filename.endswith ('.gef'):
        if __name__=="__main__":
            g=GEF()
            g.readFile(filename)
            g.plt(filename)
