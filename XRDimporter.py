# -*- coding: utf-8 -*-
"""
XRDimporter

Created on Wed Feb 13 16:19:37 2019

@author: epogue1
"""
import numpy as np
import pandas as pd

def importXRDData(fname):
    d=pd.read_csv(fname, delim_whitespace=True, names=['2Theta', 'I'])
    d['I-norm']=d['I']/max(d['I'])
    return d
    
    
def importICSDpowder(fname):
    d=pd.read_csv(fname, delim_whitespace=True, header=0)
    d['I-norm']=d['INTENSITY']/max(d['INTENSITY'])
    return d

def importVESTApowder(fname):
    d=pd.read_csv(fname, delim_whitespace=True, names=['2Theta', 'I', 'blank'])
    d['I-norm']=d['I']/max(d['I'])
    return d

def scaleXY(d, y, off):
    d['I-norm']=d['I-norm']*y-off
    return d

def importICSDxy(fname):
    d=pd.read_csv(fname, header=None, names=['2Theta', 'I'])
    d['I-norm']=d['I']/max(d['I'])
    return d
