# -*- coding: utf-8 -*-
"""
ElkForPlotting
Created on Fri Feb 22 11:47:03 2019
@author: epogue1

This code allows you to plot Elk band structures using Python

Example code is at the bottom of this file describing how to use it. 
    As-is, you can run it as a script (in the working directory containing your Elk data) and modify it as needed.
        YOU SHOULD NOT NEED TO ADJUST THE FUNCTIONS!!!
        Requires elk.in, BANDLINES.OUT, BAND.OUT at minimum
    If you put this code in the same folder as these files, it will work.
    If you write in the file path as path (ending with //) when you call these functions, it will work.
    If you want to call these functions from anywhere, you need to save it in a special location. 
        To find this location:
            type in IPython console: import sys
            press Enter
            type in IPython console: print(sys.path)
        You should get a list of file locations. You will be able to import your code 
        if you put it in any of the locations listed.
        
        With this variant, you can import the module : import ElkForPlotting as efp
            To call functions in the module, use efp.readElkInfo() or efp.readBandStructureInfoPartials()
        You can also import functions: from ElkForPlotting import readElkInfo()
            This does not change your syntax for calling the functions but you might need to import more functions.
            
"""

import pandas as pd
import scipy.constants as spc
import matplotlib.pyplot as plt
import numpy as np
import inspect
from matplotlib import rc_params
import os, sys
from matplotlib.collections import LineCollection
import matplotlib as mpl
from itertools import cycle

class ElkDataSet:
    """
    
    This object contains all the information from a given dataset related to E v. k.
    disp is the dispersion relationship, a Pandas DataFrame where
        the first column is the x-values and 
        the later columns are each energy bands in eV
    spaceGroup is the integer space group of the dataset (obtained by reading elk.in)
    S is a list of the element/species numbers for each atom in the compound
    A is a list of the atom identifiers so S[n], A[n] uniquely defines an atom
    elements is a list of the actual elements that S corresponds to
        The length of S, A, and elements equal the number of unique atoms in the compound
    description is an easier to read version of S, A, and Elements
        It is a list of strings where a given string shows the element name concatenated to A
    tot, s, p, d, and f are lists of pandas DataFrames that give the raw projections of the electron density (wavefunction) of the state at k,E onto an atomic basis
        Therefore, it doesn't capture electron density between atoms
        Each item in the list corresponds to each atom in "description" (or S,A, if you prefer)
    totNorm, sNorm, pNorm, dNorm, and fNorm are are tot, s, p, d, and f divided by the sum of all contributions at that wave vector and energy
        This is useful for seeing the percent contributions of atomic orbitals to states, 
        particularly when density has leaked into interstitial regions
    spinPolarize denotes whether these are spin-polarized calculations. If this is the case, The energy states are doubled, with one for each polarization
        Future code will plot the different polarizations differently
    hasDecomp is a boolean that is True if this object was made to include partials and False if it wasn't
    SpecPoints is a list of k-values at which special points occur   
    
    """
    def __init__(self, disp, spaceGroup=int(1), name='', S=[], A=[], elements=[], tot=[pd.DataFrame(columns=['tot'])], s=[pd.DataFrame(columns=['s'])], p=[pd.DataFrame(columns=['p'])], d=[pd.DataFrame(columns=['d'])], f=[pd.DataFrame(columns=['f'])], hasDecomp=False, spinPolarize=False, SpecPoints=[]):
        #all of these are PANDAS data frames
        self.disp=disp
        self.tot=tot
        self.s=s
        self.p=p
        self.d=d
        self.f=f
        self.S=S
        self.A=A
        self.name=name
        self.spaceGroup=spaceGroup
        self.elements=elements
        self.spinPolarize=spinPolarize
        #This sets whether 
        self.hasDecomp=hasDecomp
        self.SpecPoints=SpecPoints
        if len(elements)>0:
            self.description=[m+str(n) for m,n in zip(elements,A)]
            L=len(tot)
            g=tot[0]
            h=s[0]
            j=p[0]
            m=d[0]
            n=f[0]
            blankFrame=pd.DataFrame(0,index=g.index, columns=g.columns)
            
            totNorm=[]
            sNorm=[]
            pNorm=[]
            dNorm=[]
            fNorm=[]
            sSumNorm=pd.DataFrame(0,index=h.index, columns=h.columns)
            pSumNorm=pd.DataFrame(0,index=j.index, columns=j.columns)
            dSumNorm=pd.DataFrame(0,index=m.index, columns=m.columns)
            fSumNorm=pd.DataFrame(0,index=n.index, columns=n.columns)
            for k in range(L):
                blankFrame=blankFrame.add(tot[k])
            blankFrames=blankFrame.copy()
            blankFramep=blankFrame.copy()
            blankFramed=blankFrame.copy()
            blankFramef=blankFrame.copy()
            
            blankFrames.columns=h.columns
            blankFramep.columns=j.columns
            blankFramed.columns=m.columns
            blankFramef.columns=n.columns
            for k in range(L):
                totNorm.append(tot[k].divide(blankFrame))
                sNorm.append(s[k].divide(blankFrames))
                pNorm.append(p[k].divide(blankFramep))
                dNorm.append(d[k].divide(blankFramed))
                fNorm.append(f[k].divide(blankFramef))
                
                sSumNorm=sSumNorm.add(sNorm[k])
                pSumNorm=pSumNorm.add(pNorm[k])
                dSumNorm=dSumNorm.add(dNorm[k])
                fSumNorm=fSumNorm.add(fNorm[k])
                #build up sum of s, p, d, f contributions
                #sTotNorm=sTotNorm.add(sNorm[k])
                #pTotNorm=pTotNorm.add(pNorm[k])
                #dTotNorm=dTotNorm.add(dNorm[k])
                #fTotNorm=fTotNorm.add(fNorm[k])
            self.totSum=blankFrame
            self.totNorm=totNorm
            self.sNorm=sNorm
            self.pNorm=pNorm
            self.dNorm=dNorm
            self.fNorm=fNorm
            #sum of contributions
            self.sSumNorm=sSumNorm
            self.pSumNorm=pSumNorm
            self.dSumNorm=dSumNorm
            self.fSumNorm=fSumNorm
            
def get_var_name(var):
    """
    
    Returns the name of the variable entered here
    
    """
    callers_local_vars=inspect.currentframe().f_back.f_locals.items()
    d=[k for k, v in callers_local_vars if v is var]
    return d[0]
# read in Elk info to return list of possible elements, S, and A
#readElkInfo() is complete, returns [elements, S, A] for given file
def readElkInfo(path='.'):
    """
    
    Reads in important data from Elk file used to parse other files. (for example number of atoms)
    
    """
    
    elkin=pd.read_csv(path+'elk.in', header=None, names='Lines')
    line=elkin['L'][0]
    nLinesElkin=len(elkin)
    cnt=0
    elements=[]
    S=[]
    A=[]
    while line[0:5]!='atoms':
        cnt=cnt+1
        line=elkin['L'][cnt]
    line=elkin['L'][cnt+1]
    cnt=cnt+2
    nAt=int(line[1:5])
    atCnt=0
    specCnt=0
    aCnt=0
    while cnt<nLinesElkin:
        if elkin['L'][cnt][0]=='\'':
            line=elkin['L'][cnt]
            n=1
            while line[n]!='.':
                n=n+1
            #pull out relevant data for the first point for a given element
            aCnt=0
            numSpec=int(elkin['L'][cnt+1][1:20])
            for k in range(numSpec):
                S.append(specCnt+1)
                A.append(aCnt+1)
                elements.append(elkin['L'][cnt][1:n])
                aCnt=aCnt+1
                atCnt=atCnt+1
            #create info for atoms of a given element from line below element name
            cnt=cnt+1+numSpec+1
            specCnt=specCnt+1
        else:
                cnt=cnt+1
    cnt=0
    nLinesElkin=len(elkin)
    
    while cnt<nLinesElkin and str(elkin['L'][cnt][0:7])!='!  Herm':
        
        cnt=cnt+1
    cnt=cnt+3
    
    #cast back to int
    L=elkin['L'][cnt]
    spaceGroup=int(L[-4:])
    return [elements, S, A, spaceGroup]
def readBandStructureInfoPartials(path='.'):
    """
    imports Elk Data and initialize all values that will be stored in eventual band structure
    This should only be used if you decomposed the data into atomic and orbital contributions(task 21)
    
    """
    #path='.'
    elkin=readElkInfo(path)
    elements=elkin[0]
    S=elkin[1]
    A=elkin[2]
    hasDecomp=True
    spaceGroup=elkin[3]
    tot=[None]*len(A)
    s=[None]*len(A)
    p=[None]*len(A)
    d=[None]*len(A)
    f=[None]*len(A)
    dataChopped=pd.DataFrame()
    
    
    files=os.listdir(path)
    
    
    first=True
    for m in range(len(files)):
        fname=files[m]
        if fname[0:4]=='BAND':
           if len(fname)==18: 
               dataset=pd.read_csv(path+fname, delim_whitespace=True, header=None, names=['x', 'E-hartrees', 'tot', 's', 'p', 'd', 'f'])
               dataset['E-eV']=dataset['E-hartrees']/spc.physical_constants['electron volt-hartree relationship'][0]
               #Prepare to chop into E-eV, total, s, p, d, and f datasets
               k=True
               cnt1=0
               while k==True:
                   cnt1=cnt1+1
                   if dataset['x'][cnt1]<dataset['x'][cnt1-1]:
                       k=False
               numColTot=len(dataset.index)
               numColNew=int(numColTot/cnt1)
                
            
               #determine which element it is and where to put it---it goes at idxmin
               sFile=int(fname[6:8])
               aFile=int(fname[10:14])
               q=0
               while S[q]<sFile:
                   q=q+1
               idxmin=q
               while A[idxmin]!=aFile:
                   idxmin=idxmin+1
                    
               #start chopping the data
               dataChopped['x']=dataset['x'][0:cnt1].values
               tot[idxmin]=pd.DataFrame()
               s[idxmin]=pd.DataFrame(index=range(cnt1))
               p[idxmin]=pd.DataFrame()
               d[idxmin]=pd.DataFrame()
               f[idxmin]=pd.DataFrame()
               for n in range(numColNew):
                   #print('n=' +str(n))
                   dataChopped['E-eV' + str(n)]=dataset['E-eV'][n*cnt1:n*cnt1+cnt1].values
                   
                   tot[idxmin]['tot'+str(n)]=dataset['tot'][n*cnt1:n*cnt1+cnt1].values
                   s[idxmin]['s'+str(n)]=dataset['s'][n*cnt1:n*cnt1+cnt1].values
                   p[idxmin]['p'+str(n)]=dataset['p'][n*cnt1:n*cnt1+cnt1].values
                   d[idxmin]['d'+str(n)]=dataset['d'][n*cnt1:n*cnt1+cnt1].values
                   f[idxmin]['f'+str(n)]=dataset['f'][n*cnt1:n*cnt1+cnt1].values
                
    specPointsData=pd.read_csv(path+'BANDLINES.OUT', delim_whitespace=True, header=None, names='Coords') 
    numSpecPoints=int(len(specPointsData.index)/2)
    specPlstX=np.zeros((numSpecPoints, 1))
    for n in range(numSpecPoints):
        specPlstX[n]=specPointsData['C'][n*2]
                
    #Make ElkDataSet object for storing             
    d=ElkDataSet(dataChopped, spaceGroup=spaceGroup, S=S, A=A, elements=elements, tot=tot, s=s, p=p, d=d, f=f, hasDecomp=hasDecomp, SpecPoints=specPlstX)
    return d
#Works, returns ElkDataSet object with only dispersion disp relation defined (all else defaults)
def readBandStructureInfoNoPartials(path='.'):
    """
    
    This reads in the Band structure in the folder given by path. 
    Specifically, it looks for the BAND.OUT file and reads that in.
    It returns an ElkDataSet object containing this data as described above.
    This should only be used if you did not decompose the data into atomic and orbital contributions
    (task 20 with no task 21)
    
    """
 

    dataset=pd.read_csv(path+'BAND.OUT', delim_whitespace=True, header=None, names=['x', 'E-hartrees'])
    dataset['E-eV']=dataset['E-hartrees']/spc.physical_constants['electron volt-hartree relationship'][0]
    
     #chop dataset to include each band separately
    k=True
    cnt=0
    while k==True:
        cnt=cnt+1
        if dataset['x'][cnt]<dataset['x'][cnt-1]:
            k=False
    numColTot=dataset.size/3
    numColNew=int(numColTot/cnt)
    dataChopped=pd.DataFrame(data=dataset['x'][0:cnt], columns=['x'])
    for n in range(numColNew):
            f=dataset['E-eV'][n*cnt:n*cnt+cnt]
            dataChopped['E-eV' + str(n)]=f.values
    elkin=readElkInfo()
    spaceGroup=elkin[3]
    specPointsData=pd.read_csv(path+'BANDLINES.OUT', delim_whitespace=True, header=None, names='Coords') 
    numSpecPoints=int(len(specPointsData.index)/2)
    specPlstX=np.zeros((numSpecPoints, 1))
    for n in range(numSpecPoints):
        specPlstX[n]=specPointsData['C'][n*2]
    d=ElkDataSet(dataChopped, spaceGroup=spaceGroup, SpecPoints=specPlstX)
    return d
def labelXAxis(axs, spaceGroup, specPointsX):
    """
    
    This labels the X-axis with the default label set for your space group. Check the McQueen website for the defaults
    axs is the axis handle, spaceGroup is the integer-type spaceGroup number, and specPointsX is the set of points from Elk
    (x axis values) at which the special points are found.
    
    """
    #triclinic P
    if spaceGroup<3:
        xLbl=['$B$', '$\Gamma$', '$F$', '$\Gamma$','$G$']
    elif spaceGroup<16:
        monoCset=[5, 8, 9, 12, 15]
        if spaceGroup in monoCset:
            #mono C
            xLbl=['$V$', '$Z$','$\Gamma$', '$W$', '$L$', '$K$', '$G$']
        else:
            #mono P
            xLbl=['$\Gamma$','$Z$', '$B$', '$D$', '$Y$', '$C$', '$A$', '$E$']
    elif spaceGroup<75:
        orthoCAset=[20, 21, 35, 36, 37, 38, 39, 40, 41, 63, 64, 65, 66, 67, 68]
        orthoFset=[22, 42, 43, 69, 70]
        if spaceGroup in orthoCAset:
            #ortho C
            xLbl=['$\Gamma$', '$Z$', '$T$', '$Y$', '$\Gamma$', '$S$', '$R$']
        elif spaceGroup in orthoFset:
            #ortho F
            xLbl=['$Z$', '$\Gamma$', '$X$', '$L$', '$\Gamma$', '$Y$']
        else:
            #ortho P, I
            xLbl=['$\Gamma$','$Z$', '$T$', '$Y$', '$\Gamma$', '$X$', '$S$', '$R$', '$U$']
    elif spaceGroup<143:
        #tet P
        xLbl=['$\Gamma$','$X$', '$M$', '$\Gamma$', '$Z$', '$R$', '$A$', '$M$']
    elif spaceGroup<168:
        #trigonal or rhombohedral
        #select only rhombohedral ones
        rhomboSet=[148, 155, 160, 161, 166, 167]
        if spaceGroup in rhomboSet:
            #rhombo
            xLbl=['$L$', '$Z$', '$\Gamma$','$F$']
        else:
            #trigonal P
            xLbl=['$\Gamma$', '$K$',  '$M$', '\Gamma$', '$A$', '$L$', '$H$','$A$']
    elif spaceGroup<194:
        #hex
            xLbl=['$\Gamma$','$K$', '$M$', '$\Gamma$', '$A$', '$L$', '$H$', '$A$']
    else:
        #all cubic groups
        cubIset=[197, 199, 204, 206, 211, 214, 217, 220, 229, 230]
        cubPset=[195, 198, 200, 201, 205, 207, 208, 212, 213, 215, 218, 219, 221, 222, 223, 224]
 
        if spaceGroup in cubIset:
            #cubic I,
            xLbl=['$Gamma$', '$N$', '$P$', '$\Gamma$', '$H$', '$N$']
        elif spaceGroup in cubPset:
            #cub P
            xLbl=['$\Gamma$','$M$', '$X$', '$\Gamma$', '$R$']
        else:
            #cub F
            xLbl=['$L$', '$\Gamma$','$X$', '$W$', '$L$', '$K$', '$\Gamma$']
    #put the labels on the plot corresponding to xLbl
    [plt.axvline(_x, linestyle=':', color='k', alpha=0.5) for _x in specPointsX]
    axs.xaxis.set_major_locator(plt.FixedLocator(specPointsX))
    
    axs.set_xticklabels(xLbl)
    return xLbl
                

                
#need to get this working
def endPlot(dataForSpaceGroup,full=False, Emin=-3, Emax=3, fontSize=10):
    """
    
    This ends the plot by putting x and y lables on the plot. 
    It also labels the X-axis using the default labels for the spacegroup given in the sample
    ElkDataSet object named dataForSpaceGroup. This is also where you determine the range in energies you want to plot.
    All energies are relative to the Fermi energy
    
    
    """
    
    axs=plt.gca()
    if full==False:
        #ax=plt.gca()
        axs.set_ylim([Emin, Emax])
    spaceGroup=dataForSpaceGroup.spaceGroup
    plt.xlabel('k', fontsize=fontSize)
    plt.ylabel('E (eV)', fontsize=fontSize)
    plt.axhline(y=0, linestyle=':', color='k', alpha=0.5)
    plt.tight_layout(pad=0)
    plt.margins(x=0)
    labelXAxis(axs, spaceGroup, dataForSpaceGroup.SpecPoints)
    plt.show()
    
def plotScaledPartials(importedData, scaleToMap, thickScale=4, newEntry='', lgndOnSide=True, legsize=8, lgndList=[] ):
    """
    Plot orbital and/or atomic contributions. This should be used after beginPlot() and before endPlot()
    importedData is your ElkDataSet object and scaleToMap is what determines the thickness of the plotted lines
    
    newEntry describes the name you want to appear for this entry in the legend.
    Be sure to include lgndList if there were previous lgnd entries to ensure that newEntry is added to the old legend 
    and the old legend is not deleted
    
    lgndOnSide will place the legend on the right side of the plot.
    
    This function returns a list in the same format for lgndList that can be used subsequently the same way as in this function.
    
    """
    
    axs=plt.gca()
    if len(lgndList)==2:
        lcs=lgndList[0]
        lgnd=lgndList[1]
    else:
        lcs=[]
        lgnd=[]
    
    #leg=[c for c in axs.get_children() if isinstance(c, mpl.legend.Legend)]
    #print(len(leg))
    #if len(leg)>1:
        #print(leg[1])
        #print(leg.get_handles())
    #lcs, lgnd=axs.get_legend_handles_labels()
    #print(leg)
    #color_codes=map('C{}'.format, cycle(range(10)))
    #color_codes=next(color_codes)
    
    
    d=importedData.disp
    x=d['x'].values
    colNames=scaleToMap.columns
    L=len(colNames)
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    idx=len(lcs)
    for n in range(L):
        y=d['E-eV'+str(n)].values
        lwidths=scaleToMap[colNames[n]].values
        lwidths=scaleToMap[colNames[n]].values*thickScale
        points=np.array([x,y]).T.reshape(-1,1,2)
        segments=np.concatenate([points[:-1], points[1:]], axis=1)
        lc=LineCollection(segments, alpha=0.7, linewidths=lwidths, color=cycle[idx])
        
        axs.add_collection(lc)

        if n==0:
            lcs.append(lc)
            lgnd.append(newEntry)
            
    #print(len(lcs))
    if lgndOnSide==True:               
        leg=plt.legend(lcs, lgnd, bbox_to_anchor=(1.02, 1), loc=2, fontsize=legsize)
    else:
        leg=plt.legend(lcs, lgnd,loc=2, fontsize=legSize)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
        
    #print(lcs)
    #print(lgnd)
    #print([c for c in axs.get_children() if isinstance(c, mpl.legend.Legend)])
    return [lcs, lgnd]
        #label=label.append()
        #return

            #y=dataChopped['E-eV' + str(n)].values
            #lwidths=dataChopped[whatToPlot+str(n)].values
            #points=np.array([x,y]).T.reshape(-1,1,2)
            #segments=np.concatenate([points[:-1], points[1:]], axis=1)
            #if makeTotBlack==True:
            #    lc=LineCollection(segments, colors='k', alpha=ap, linewidths=1)  
            #if n==0:
            #    lcs.append(lc)
            #    lgnd.append('total')

    #returns a list of arrays containing the scale
    
def beginPlot(title='', width=3.25, height=2.25, titleSize=12):
    """used to start plotting--the canvas on which plots are made. Set the title here by title="This is my title", similarly with height, width, etc.  """
    fig=plt.figure(figsize=(width, height))
    ax=plt.gca()
    plt.title(title, fontsize=12)
    return fig
    


#fully working
def plotDispersion(dataObj):
    """
    #dataObj is the ElkDataSet object created by import functions
    #Use this if you did not decompose into band contributions or just want to plot energies
    """
    ax=plt.gca()
    d=dataObj.disp
    d.plot(x='x', color='C0', legend=None, ax=ax)    
    return ax
       

def addPartials(ElkDataSetObject, descriptionsToSum):
    """
    Takes ElkDataSet object "dataset" and adds up contributions of descriptions to sum
    descriptionsToSum is a list containing [ElkDataSetObject.[tot, s, p, d, f] ElkDataSetObject.Description]
    for example, [[Ba2NaIO6.tot 'O1'],
                  [Ba2NaIO6.tot 'O2']]
    
    """
    L=len(descriptionsToSum)
    d=descriptionsToSum[0][0][0]
    [r, c]=d.shape
    descriptions=ElkDataSetObject.description
    #set up dataFrame
    blankFrame=pd.DataFrame(0,index=d.index, columns=d.columns)
    #loop through all things to add together
    for k in range(L):
        #find where ElkDataSetObject.Description is
        toLook=descriptionsToSum[k]
        data=toLook[0]
        lookFor=toLook[1]
        
        idx=descriptions.index(lookFor)
        toAdd=data[idx]
        #add corresponding ElkDataSetObject to blankFrame
        blankFrame=blankFrame.add(toAdd)
        
    return blankFrame


    """Takes ElkDataSet object "dataset" and normalizes all values of tot, s, p, d, and f by the total sum of all tot atomic contributions
    This is useful for seeing the percent contributions of atomic orbitals to states, particularly when density has leaked into interstitial
    regions
    This data is stored in the ElkDataSetObject as items totNorm, sNorm, pNorm, dNorm, and fNorm
    
    """
#Ba2NaIO6=readBandStructureInfoPartials(path='.')
#import functions
Ba2AgIO6=readBandStructureInfoPartials(path='C:\\Users\\epogue1\\Documents\\postdoc\\Calculations\\Elk\\Ba2AgIO6\\')
Ba2AgIO6.name=get_var_name(Ba2AgIO6)


#Plot Dispersion Relationship
#beginPlot("$Ba_2NaIO_6$")
#plotDispersion(Ba2NaIO6)
#endPlot(Ba2NaIO6)

#Plot Partials
beginPlot("$Ba_2AgIO_6$")
#I'm plotting the total density of states for S=2 and A=1 (aka the Na atom)
dLeg=[]
dLeg=plotScaledPartials(Ba2AgIO6, Ba2AgIO6.totNorm[2], newEntry=Ba2AgIO6.description[2])# Ba2NaIO6.description[3])
#dLeg=plotScaledPartials(Ba2AgIO6, Ba2AgIO6.tot[6], newEntry=Ba2AgIO6.description[6], lgndList=dLeg)#
#dLeg=[]
#dLeg=plotScaledPartials(Ba2AgIO6, Ba2AgIO6.totNorm[3], newEntry=Ba2AgIO6.description[3], lgndList=dLeg)#
a=[[Ba2AgIO6.totNorm, 'O1'],
   [Ba2AgIO6.totNorm, 'O2'],
   [Ba2AgIO6.totNorm, 'O3'],
   [Ba2AgIO6.totNorm, 'O4'],
   [Ba2AgIO6.totNorm, 'O5'],
   [Ba2AgIO6.totNorm, 'O6']]
Ba2AgIO6AllO=addPartials(Ba2AgIO6, a)
b=[[Ba2AgIO6.tot, 'Ba1'],
   [Ba2AgIO6.tot, 'Ba2']]
#Ba2AgIO6AllBa=addPartials(Ba2AgIO6, b)
dLeg=plotScaledPartials(Ba2AgIO6, Ba2AgIO6AllO, newEntry="all O", lgndList=dLeg)
# =============================================================================
# totSum=[[Ba2AgIO6.tot, 'Ba1'],
#    [Ba2AgIO6.tot, 'Ba2'],
#    [Ba2AgIO6.tot, 'Ag1'],
#    [Ba2AgIO6.tot, 'I1'],
#    [Ba2AgIO6.tot, 'O1'],
#    [Ba2AgIO6.tot, 'O2'],
#    [Ba2AgIO6.tot, 'O3'],
#    [Ba2AgIO6.tot, 'O4'],
#    [Ba2AgIO6.tot, 'O5'],
#    [Ba2AgIO6.tot, 'O6']]
# totAddedTogether=addPartials(Ba2AgIO6, totSum)
# dLeg=plotScaledPartials(Ba2AgIO6, totAddedTogether, newEntry="all", lgndList=dLeg)
# =============================================================================
#dLeg=plotScaledPartials(Ba2AgIO6, Ba2AgIO6.tot[3], newEntry=Ba2AgIO6.description[3], lgndList=dLeg)
endPlot(Ba2AgIO6)

