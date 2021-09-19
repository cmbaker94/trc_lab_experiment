# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 12:51:17 2021

@author: cmbaker9

set up code to calculate the directional spectra
"""
import numpy as np
import sys
from calc_dir_sprd import calc_dir_sprd
sys.path.append("C:\\Users\\cmbaker9\\Documents\\MATLAB\\MTOOLS\\pyDIWASP")
from dirspec import dirspec

rho = 1000
g = 9.81

def calc_SM(Tinfo,inst,df,XYZ):
    SM = dict()
    EP = dict()
    ID = dict()
    
    # define instrument data type
    if inst == 'wg' or inst == 'cam' or inst == 'lid': meas = 'elev'
    elif inst == 'press': meas = 'pres'; df = df/(rho*g)
    ID['datatypes'] = [meas]
    for i in range(0,XYZ.shape[1]-1): ID['datatypes'].append(meas)
    
    if inst == 'wg' or inst == 'cam' or inst == 'lid':
        ID['depth'] = float(XYZ[2].mean()) # mean water depth of instrument
    else:
        ID['depth'] = float(Tinfo['h']-XYZ[2].mean()+0.05)
    ID['data'] = df.values # copy data
    ID['layout'] = XYZ
    ID['layout'][2,:] = Tinfo['h']-ID['layout'][2,:]
    if inst == 'cam': 
        ID['fs'] = 8
    elif inst == 'lid': 
        ID['fs'] = 10
    else: 
        ID['fs'] = 100
    
    # Spectral matrix
    nf = 301
    nd = 360
    SM['dirs'] = np.linspace(0,359,nd)
    SM['freqs'] = np.linspace(0,3,nf) 
    SM['S'] = np.zeros([nf,nd])
    
    # Resolution of estimation
    if inst == 'cam' or inst == 'lid': EP['nfft'] = 8*(2**4)
    else: EP['nfft'] = 100*(2**4)
    EP['dres'] = 360
    EP['method'] = 'BDM'#'IMLM'
    
    # Options_ = 'PLOTTYPE'=0
    Options_ = []
    [SMout,EPout] = dirspec(ID,SM,EP,Options_)
    SMout = calc_dir_sprd(SMout)
    return SMout, EPout