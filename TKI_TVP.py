import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.offsetbox as offsetbox
import os
from timml import*

params = {'font.family': 'sans-serif',
          'font.sans-serif': 'arial',
          'axes.labelsize': 10,
          'axes.facecolor': '#ffffff', 
          'axes.labelcolor': 'black',
          'legend.fontsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'lines.linewidth': 1,
          'grid.color': 'grey',
          'grid.linestyle': 'dashed',
          'grid.linewidth': 0.5,
          'text.usetex': False,
          'font.style': 'normal',
          'font.variant':'normal',
          'figure.facecolor': 'white',
          'font.size':8,
          'figure.autolayout': True,
          'figure.figsize': (8,8),
          'figure.dpi': 100,
          }
plt.rcParams.update(params)

"""
Hieronder wordt aan de hand van een sondering de input voor de TIM berekening gegenereerd.
Dit is natuurlijk optioneel, maar als de uitdraai van deze bewerking wordt vergeleken
met de opgave van REGIS zijn de verschillen te billijken. 
Deze sondeerbewerking is op basis van een script van Rob van Putten van een aantal
jaar geleden.

Deze bewerking:
    depth  npor     k  kzkh  
0    9.68  0.32  6.59  0.09   Sterksel
1    8.68  0.35 11.93  0.14   
2    7.68  0.24  5.98  0.09
3    6.68  0.14  0.68  0.06   Waalre / Tegelen 
4    5.68  0.28  2.35  0.07
5    4.68  0.35  6.37  0.09   
6    3.68  0.35  9.00  0.11
7    2.68  0.34  8.92  0.11
8    1.68  0.35 10.35  0.12    
9    0.68  0.22  2.11  0.07
10  -0.32  0.28  5.92  0.09
11  -1.32  0.35 12.86  0.15   
12  -2.32  0.35  4.33  0.08
13  -3.32  0.35  2.03  0.07
14  -4.32  0.33  4.97  0.08
15  -5.32  0.31  7.16  0.10
16  -6.32  0.35 10.22  0.12
17  -7.32  0.35  8.55  0.11
18  -8.32  0.28  9.62  0.12
19  -9.32  0.30 18.46  0.21    
20 -10.32  0.35 14.85  0.17
21 -11.32  0.35 10.91  0.13
22 -12.32  0.35  6.98  0.10
23 -13.32  0.33  6.03  0.09
24 -14.32  0.33 17.79  0.20
25 -15.32  0.30 27.04  0.39
26 -16.32  0.35 20.14  0.24
27 -17.32  0.35 16.06  0.18
28 -18.32  0.31 25.11  0.34
29 -19.32  0.26 36.61  0.76   Top Maassluis
30 -35                        Hydrologische basis
"""

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
                         
        pw = float(args[self._columns[3]]) 
        if pw<-10:
            pw=0.1

        self.dz.append(dz)
        self.qc.append(qc)
        self.pw.append(pw)
        if qc<=0.001:
            qc=0.1
            self.wg.append(10.)
        else:
            wg = abs((pw / qc) * 100.)
        wg = abs((pw / qc) * 100.)

##############################################   K-waarde berekeing     
        if wg>=0.0:
            if wg >5: wg=15
            ke=np.exp(wg)
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
        return pd.DataFrame(data=a, columns=['depth', 'npor', 'k'])
        
    def plt(self, filename):
        df = self.asDataFrame()
        df = df.sort_values('depth', ascending=False)

        if df.empty:
            return df
        
        LD = 50 # laagdikte in model, gelijk aan LD*sondeerinterval, bij 2 cm is LD dus 1 m
        
        df = df.rolling(LD).mean()  # hier wordt de sondering vervormd van laagjes van 2 cm naar een laag van 1 meter
        df = df.iloc[:: LD]
        df = df.dropna()
        df = df.reset_index(drop=True)
        df['k'] = df['k']*0.95 # correctie o.b.v. pompproef
        
        df['kzkh']= 2/(33.653*np.exp(-0.066*(df['k']))) # op basis van literatuurgegevens is dit algoritme gemaakt
        df['kzkh']= np.where(df['kzkh']>1, 1, df['kzkh'])

        dzend=-35
        print('Aangenomen diepte aq = ', dzend, '[m NAP]')
        
        """
        Hieronder wordt het dataframe verwerkt tot een input voor TTIM
        """
        
        dfd=df['depth']
        dfd.loc['ld']= dzend
        dfn = df.iloc[: , [1]].copy()   
        dfn = pd.concat([df, pd.DataFrame.from_records([{'depth':dzend},])], ignore_index=True)
        dfnp= pd.concat([df, pd.DataFrame.from_records([{'npor': 0.25},])],  ignore_index=True)

        hstar = 7.5 # DE GHG op deze locatie
        print('hstar = ', hstar, '  [m NAP]')
        
       
        """
        Omdat hier sprake is van een langdurige bemaling (> 10 dagen) en naar een worst case 
        scenario wordt gevraagd wordt hier stationair gerkend (TIMml)
        
        """
        
        ml=Model3D(kaq=df['k'], z=dfn['depth'], kzoverkh=df['kzkh'], npor=dfnp['npor'], 
                    topboundary='semi', topres=1, topthick=0.01, hstar=hstar) 
        
        """
        De locatie van de onttrekkingsfilters en de retourfilters wordt opgehaald uit een 
        XLSX spreadsheet, niet omdat het moet maar het invoegen van nieuwe / verwijderen van
        overbodige locaties is dan eenvoudiger
        """
       
        mat=pd.read_excel('loc_.xlsx', engine='openpyxl')  
        mat=mat.replace(',','.')
        
        well      = mat[mat['Naam'].str.match ('w')]    # Dit zijn de onttrekkingsputten   
        ret       = mat[mat['Naam'].str.match ('r')]    # Dit zijn de retourputten    
        pxy       = mat[mat['Naam'].str.match ('p')]    # Het hart van het project   
        damwand   = mat[mat['Naam'].str.match ('dw')]   # De damwand om de voorbouwlocatie
        damwand2  = mat[mat['Naam'].str.match ('dx')]   # De damwand om de noordelijke moot
        Vloernoord= mat[mat['Naam'].str.match ('Vn')]   # Injectievloer aan de noordzijde
        Vloerzuid = mat[mat['Naam'].str.match ('Vz')]   # Injectievloer aan de zuidzijde
        Drainnoord= mat[mat['Naam'].str.match ('Dn')]   # drain boven injectievloer in kuip
        Drainzuid = mat[mat['Naam'].str.match ('Dz')]   # drain boven injectievloer in kuip
        Lijn1     = mat[mat['Naam'].str.match ('L1')]   # Verlagingslijn 1
        Lijn2     = mat[mat['Naam'].str.match ('L2')]   # Verlagingslijn 2
        
        x1=pxy.iloc[0,1]
        y1=pxy.iloc[0,2]
       
        """
        Hieronder de verschillende elementen in de TIMml berekening
        In de uitvoering werd 8800 m3/dag opgepomt en 75% daarvan geretourneerd
        """
        
        ##Onttrekkingsfilter
        tf =  0 # Top en bottom van het filter
        bf = -5
        nwells= len(well)
        nret  = len(ret)
        print('Aantal bronnen ', nwells)
        print('Aantal retourbronnen ', nret)
        
        tfind = df['depth'].sub(tf).abs().values.argmin()  #Bepaling laagnr voor het plaatsen van het filter op de juiste diepten vin he model
        bfind = df['depth'].sub(bf).abs().values.argmin() 
        well1 = np.arange(tfind, bfind, 1)
        for rnum in well.index:
            punt = well['Naam'][rnum]
            x    = well['x'][rnum] 
            y    = well['y'][rnum]
            punt = Well(ml, x,y,8800/19, rw=0.125, layers=well1)

        ## Retourfilters
        vstf = -15
        vsbf = -16
        vstfind = df['depth'].sub(vstf).abs().values.argmin()
        vsbfind = df['depth'].sub(vsbf).abs().values.argmin() 
        well2 = np.arange(vstfind, vsbfind, 1)
        for rnum in ret.index:
            puntr = ret['Naam'][rnum]
            x     = ret['x'][rnum] 
            y     = ret['y'][rnum]
            puntr = Well(ml, x,y, -6600/45 , rw=0.125, layers=well2)
 
        ########## Drains
        niv = 1.5
        xd = []
        yd = []
        for rnum in Drainnoord.index:
            xdn   = Drainnoord['x'][rnum] 
            ydn   = Drainnoord['y'][rnum]
            xd.append(xdn) 
            yd.append(ydn)
        drn= HeadLineSinkString(ml, xy=list(zip(xd, yd)), hls=niv, res=0, order=3, wh=1, layers=8)            

        xd1 = []
        yd1 = []
        for rnum in Drainzuid.index:
            xdz   = Drainzuid['x'][rnum] 
            ydz   = Drainzuid['y'][rnum]
            xd1.append(xdz) 
            yd1.append(ydz)
        drz= HeadLineSinkString(ml, xy=list(zip(xd1, yd1)), hls=niv, res=0, order=3, wh=1, layers=8)            
    
        ########Damwanden
        td = self.z
        od = -7
        tdf = df['depth'].sub(td).abs().values.argmin()
        odf = df['depth'].sub(od).abs().values.argmin() 
        wandlaag=np.arange(tdf, odf,1)  # wanden moeten per laag conform de aangebrachte diepten

        xp = []
        yp = []
        for rnum in damwand.index:
            xdw = damwand['x'][rnum]
            ydw = damwand['y'][rnum]
            xp.append(xdw) 
            yp.append(ydw)
        ld = LeakyLineDoubletString(ml, xy=list(zip(xp, yp)), res=1000, layers=wandlaag, order=3)        

        xp2 = []
        yp2 = []
        for rnum in damwand2.index:
            xdw = damwand2['x'][rnum]
            ydw = damwand2['y'][rnum]
            xp2.append(xdw)
            yp2.append(ydw)
        ld = LeakyLineDoubletString(ml, xy=list(zip(xp2, yp2)), res=1000, layers=wandlaag, order= 3)        

        """
        In beide bouwkuipen is een injectielaag aangebracht, deze wordt hier ingebracht
        met de inhom module 
        """
        
        botinj = -6.5
        topinj = -5.5
        df.loc[(df['depth'] >= botinj) & (df['depth'] <= topinj), 'k']    = 1e-2
        df.loc[(df['depth'] >= botinj) & (df['depth'] <= topinj), 'npor'] = 1e-2
        df.loc[(df['depth'] >= botinj) & (df['depth'] <= topinj), 'kzkh'] = 1e-2
        """"
        Hieronder het inhom element met de injectielaag.
        Buiten de injectievlakken is dus de 'originele' df nog geldig, binnen de vlakken is deze
        vigerend
        
        11  -1.32  0.35 12.86  0.15
        12  -2.32  0.35  4.33  0.08
        13  -3.32  0.35  2.03  0.07
        14  -4.32  0.33  4.97  0.08
        15  -5.32  0.31  7.16  0.10
        16  -6.32  0.01  0.01  0.01  De injectielaag met een weerstand van 100 dagen
        17  -7.32  0.35  8.55  0.11
        18  -8.32  0.28  9.62  0.12
        19  -9.32  0.30 18.46  0.21
        20 -10.32  0.35 14.85  0.17
        21 -11.32  0.35 10.91  0.13
        
        """
        
        xvn = []
        yvn = []
        for rnum in Vloernoord.index:
            xi = Vloernoord['x'][rnum]
            yi = Vloernoord['y'][rnum]
            xvn.append(xi) 
            yvn.append(yi)
        VlN = PolygonInhom3D(ml, xy=list(zip(xvn, yvn)), kaq=df['k'], z= dfn['depth'],
                                kzoverkh = df['kzkh'], npor=dfnp['npor'], 
                                topboundary='semi', hstar = hstar,  ndeg=10)        

        xvz = []
        yvz = []
        for rnum in Vloerzuid.index:
            xi = Vloerzuid['x'][rnum]
            yi = Vloerzuid['y'][rnum]
            xvz.append(xi) 
            yvz.append(yi)
        VlZ = PolygonInhom3D(ml, xy=list(zip(xvz, yvz)), kaq=df['k'], z= dfn['depth'],
                                kzoverkh = df['kzkh'], npor=dfnp['npor'], 
                                topboundary='semi', hstar = hstar,  ndeg=10)        

# ####################################
        ml.solve()
# ###################################

        """
        Hieronder worden de debieten berekend
        """
    
        qtot=int(punt.discharge().max())*nwells
        qret=int(puntr.discharge().max())*nret
        qdrz=int(drz.discharge().max())
        qdrn=int(drn.discharge().max())

        print('Bemalingsdebiet_1 dag ', qtot, ' [m3]')
        print('Retourdebiet_1 dag ', qret*1,  ' [m3]')
        print('Draindebiet_1 dag ', qdrn+qdrz, ' [m3]')
        print('Uurdebiet ',round((qtot+qdrn+qdrz)/24, 0), ' [m3/uur]')
        print('KD waarde WVP ', round((int(df['k'].sum()+(abs(dzend)-abs(df.iloc[-1,0]))*df.iloc[-1,2])),0), '  [m2/dag]')
        
        """
        Vervolgens worden de contouren gemaakt.
        Hier kan de  laag waarin de druk moet worden getoond (hier laag 8, het schuifniveau)
        en de contourniveau's. Deze kunnen onbeperkt worden aangepast.
        Verder kan de grootte van het door te rekenen oppervlakte worden aangegeven 
        (hier 400* 400 meter) en de pixelgrootte (200 pixels = 2*2 meter)
        
        """
        ml.contour(win=[x1-200,x1+200,y1-200,y1+200], ngr=200,
                    layers=[8], levels=np.arange(-8, 8, 0.5), 
                    layout=True, labels=True, decimals=1, 
                    color='blue', 
                    newfig=False, figsize=(10,10))
        
        ### Overzichtstekening met contouren, export als plaatje

        plt.fill(Vloernoord['x'], Vloernoord['y'], 'orange')
        plt.plot(damwand2['x'], damwand2['y'], 'r-', linewidth=2 )
        plt.plot(damwand['x'], damwand['y'],  'r-', linewidth=2, label='_nolegend_')
        plt.plot(well['x'], well['y'], 'rx')
        plt.plot(ret['x'], ret['y'], 'g*' )
        plt.plot(Drainnoord['x'], Drainnoord['y'], 'k--', linewidth=2)
        
        plt.plot(Lijn1['x'], Lijn1['y'], 'g-')
        plt.plot(Lijn2['x'], Lijn2['y'], 'g-', label='_nolegend_')
       
        
        plt.plot(Drainzuid['x'], Drainzuid['y'], 'k--', linewidth=2)
        plt.fill(Vloerzuid['x'], Vloerzuid['y'], 'orange')
        plt.legend(['Injectievloeren', 'Damwanden', 'Onttrekkingsbronnen', 
                    'Retourbronnen','Drains', 'Verlagingslijnen'])
        plt.xlim(-200,200)
        plt.ylim(-200,200)
        plt.grid()
        plt.savefig('TVP_TIM'+'.png', bbox_inches = "tight")
 
######   Verlagingslijnen      
        """
        Hieronder wordt van een lijn over de x en de y richting de verlagingslijnen
        op verschillende diepten gegenereerd. OOk deze wordt als apart plaatje 
        geexporteerd
        """
        
        laag = np.arange(7,29,2)
        ontwpeil = 2.0  # Schuifniveau
 
#Doorsnede over Lijn 1
        fig, ax = plt.subplots(nrows=1, ncols=2, gridspec_kw = {'width_ratios':[1,1]})

        z=int(self.z+1)
        lijn1= ml.headalongline(np.linspace(Lijn1.iloc[0,1],Lijn1.iloc[1,1],200),np.linspace(Lijn1.iloc[0,2],Lijn1.iloc[1,2],200))  #dit is uiteindelijk een x,y dwars op de bemaling
        color=iter(cm.jet(np.linspace(0,1,len(laag))))
        for i in laag:
            c=next(color)
            ax[0].plot(np.linspace(0,400,200), lijn1[i], c=c, lw=2)
        
        ax[0].grid(axis='both')        
        ax[0].set_xlabel('Lijn_1')
        ax[0].set_ylim(z-12, z+1)
        ax[0].set_yticks(np.arange(z-12, z+1,0.5))
        ax[0].set_title('TVP_TIM', pad=10)

#Doorsnede over y
        lijn2= ml.headalongline(np.linspace(Lijn2.iloc[0,1],Lijn2.iloc[1,1],200),np.linspace(Lijn2.iloc[0,2],Lijn2.iloc[1,2],200))  #dit is uiteindelijk een x,y dwars op de bemaling
        color=iter(cm.jet(np.linspace(0,1,len(laag))))
        for i in laag:
            c=next(color)
            ax[1].plot(np.linspace(0,400,200), lijn2[i], c=c, lw=2)
        leg=ax[1].legend(df.loc[df.index[laag], 'depth'], loc=(1.1,0), fontsize = 8, title='Filter op [m NAP]')
        for legobj in leg.legendHandles:
            legobj.set_linewidth(2.0)
        
        ax[1].grid(axis='both')
        ax[1].set_xlabel('Lijn_2')
        ax[1].set_ylim(z-12,z+1)
        ax[1].set_yticks(np.arange(z-12,z+1,0.5))
  
        ax[0].hlines(ontwpeil, xmin=0, xmax=400, 
                      linestyle='solid', lw=3, color='grey')
        ax[1].hlines(ontwpeil, xmin=0, xmax=400, 
                      linestyle='solid', lw=3, color='grey')
        ax[0].hlines(self.z, xmin=0, xmax=400, 
                      linestyle='solid', lw=3, color='green')
        ax[1].hlines(self.z, xmin=0, xmax=400, 
                      linestyle='solid', lw=3, color='green')
        
        totdeb = round((qtot+qdrn+qdrz))
        quur   = round(totdeb/24,1)
        
        Text=str(totdeb)+' [m\u00b3/dag]' 
        ob = offsetbox.AnchoredText(Text, loc=3, prop=dict(color='black', size=14))
        ob.patch.set(color='lightblue', alpha=0.3)
        ax[0].add_artist(ob)       
        Text=str(quur)+' [m\u00b3/uur]' 
        ob = offsetbox.AnchoredText(Text, loc=3, prop=dict(color='black', size=14))
        ob.patch.set(color='lightblue', alpha=0.3)
        ax[1].add_artist(ob)       

        fig.savefig('Verlagingen'+'.png', bbox_inches = "tight")
        plt.show()
        
 
        
for filename in os.listdir(os.getcwd()):
    if filename.endswith ('.GEF') or filename.endswith ('.gef'):
        if __name__=="__main__":
            g=GEF()
            g.readFile(filename)
            g.plt(filename)
        plt.close()
