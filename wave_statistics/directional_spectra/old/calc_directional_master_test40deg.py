#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 16:39:08 2020

@author: cmbaker9

This code estimates the wave field mean direction and directional spread using
the sea-surface elevation measured by an array of sensors wire resistance gages 
located under the laboratory bridge in the 'inner shelf'. The code calls the 
function dirspec from the pyDIWASP (DIrectional WAve SPectra Toolbox) Toolbox.

Function input:
    - sea-surface elevation from 2 rows of 5 wire resistance gages measuring 'elevation'
    - Estimation method: IMLM: Iterated maximum liklihood method (Pawka 1983)
        Refinement of the EMLM (fast method that performs well with narrow 
        unidirectional spectra) that iteratively improves the original EMLM
        estimate. Highly dependent on the quality of the original solution so will
        tend to perform poorly in the same situations as the EMLM. Will tend to
        reduce anomalies such as negative energy in the EMLM solution.
        Computation time directly dependent on number of refining iterations but
        provides good accuracy for reasonable computing time. Can overestimate
        peaks in the directional spectra by overcorrecting the original estimate
    - Instrument data structure:
        - 'data' is 2d matrix where each column includes data from a sensor
        - 'layout' provides the x,y,z in each column
        - 'datatypes' list of sensor types (ie 'elev' for wave gages)
        - 'depth' is the mean overall depth of the measurment area (m), which 
           is close to the tank water elevation
        - 'fs' is the samplying frequency of the instrument (Hz)
    - Spectral matrix structure:
        - 'freqs' is a vector of length nf defining bin centers of the spectral
          matrix frequency axis
        - 'dirs' is a vector of length nd definding bin centers of the spectral 
           matrix direction axis
        - 'S' is a mtrix of size [nf,nd] containing the spectral density
        - 'xaxisdir' is the compass direction of the x axis from which angels are measured
    - Resolution of the estimation (EP):
        - 'nfft' and 'dres' define the maximum resolution that of spectral output
            - 'nfft' is the number of DFTs carried out in the calculation of the 
                cross-power spectra. Higher numbers = greater freq resolution.
                Bounded by SM.freqs and defaulted to 'sensible' value.
            - 'dres' is the number of directions used in the estimation calc.
                computation is carried out for a complete circle of directions.
                default of 180 with 2 deg resolution. 
        - 'smooth' is an on/off switch to determine if smoothing is applied
           to the final spectra. Default is on.
        - 'method' is the estimation method (see above)
        
Function output:
    - SMout: a spectral matrix structure containing the results
    - EPout: the estimation parameters structure with the values actually used 
       for the computation including any default settings

For more details see the manual: 
    http://152.2.92.46/dataproc/browser/DPWP/trunk/DPWP/diwasp_1_1GD/user_manual.pdf?rev=495&format=raw
"""

from IPython import get_ipython
get_ipython().magic('reset -sf')

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import numpy as np
import os
from scipy import signal
import sys
# from scipy import pi
import scipy as sp
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
# from matplotlib import rc
# font = {'family' : 'normal',
#         'weight' : 'normal',
#         'size'   : 22}
# rc('font', **font)
# rc('text', usetex=True)

# %% Add functions

def read_data(Tinfo,inst,datapath):
    nums = Tinfo[inst]
    for i in range(len(nums)):
        dname = inst + str(nums[i])
        fname = datapath + dname + '.txt'
        f = open(fname, "r")
        lines = f.readlines()
        f.close()
        x = float(lines[65].split()[2])
        y = float(lines[66].split()[2])
        if inst == 'wg': z = Tinfo['h']
        else: z = float(lines[67].split()[2])
        if i == 0: 
            temp=np.loadtxt(fname, delimiter=',', comments='%', unpack=True)
            df = pd.DataFrame(temp,columns=[dname])
            dxyz = pd.DataFrame([x, y, z],columns=[dname])
        else:
            df[dname] = pd.read_csv(fname, delimiter=',', comment='%', header=None)
            dxyz[dname] = [x, y, z]
    dXYZ = dxyz.values
    # if inst == 'press': dxyz.to_csv('innershelf_pressure_gages_xyz.csv', header=True)
    return df, dXYZ

def calc_dir_sprd(Smout):
    # Compute theta and sigma
    Sftheta       = Smout['S'].real
    Dir           = Smout['dirs']
    freq          = Smout['freqs']
    dtheta        = abs(Dir[1]-Dir[0])
    dfreq         = abs(freq[1]-freq[0])
    Sf            = np.nansum(Sftheta, axis=1)*dtheta
    Sd            = np.nansum(Sftheta, axis=0)*dfreq
    a1F           = np.matmul(Sftheta,(np.cos(np.deg2rad(Dir))*dtheta))/(np.nansum(Sftheta,axis=1)*dtheta)
    b1F           = np.matmul(Sftheta,(np.sin(np.deg2rad(Dir))*dtheta))/(np.nansum(Sftheta,axis=1)*dtheta)
    a2F           = np.matmul(Sftheta,(np.cos(np.deg2rad(2*Dir))*dtheta))/(np.nansum(Sftheta,axis=1)*dtheta)
    b2F           = np.matmul(Sftheta,(np.sin(np.deg2rad(2*Dir))*dtheta))/(np.nansum(Sftheta,axis=1)*dtheta)
    
    th_1          = np.arctan2(b1F,a1F)
    th_2          = 0.5*np.arctan2(b2F,a2F)
    sig_1         = np.sqrt(abs(2*(1-(a1F*np.cos(th_1)+b1F*np.sin(th_1)))))
    sig_2         = np.sqrt(abs(0.5*(1-(a2F*np.cos(2*th_2)+b2F*np.sin(2*th_2)))))
    
    th_1 = np.rad2deg(th_1)
    th_2 = np.rad2deg(th_2)
    sig_1 = np.rad2deg(sig_1)
    sig_2 = np.rad2deg(sig_2)
    
    Smout['Sf'] = Sf
    Smout['Sd'] = Sd
    Smout['th_2'] = th_2
    Smout['sig_2'] = sig_2
    return Smout

def calc_SM(Tinfo,inst,df,XYZ):
    SM = dict()
    EP = dict()
    ID = dict()
    
    # define instrument data type
    if inst == 'wg' or inst == 'cam': meas = 'elev'
    elif inst == 'press': meas = 'pres'; df = df/(rho*g)
    ID['datatypes'] = [meas]
    for i in range(0,XYZ.shape[1]-1): ID['datatypes'].append(meas)
    
    if inst == 'wg' or inst == 'cam':
        ID['depth'] = float(XYZ[2].mean()) # mean water depth of instrument
    else:
        ID['depth'] = float(Tinfo['h']-XYZ[2].mean()+0.05)
    ID['data'] = df.values # copy data
    ID['layout'] = XYZ
    ID['layout'][2,:] = Tinfo['h']-ID['layout'][2,:]
    if inst == 'cam': ID['fs'] = 10
    else: ID['fs'] = 100
    
    # Spectral matrix
    if inst == 'cam': nf = 3001
    else: nf = 3001
    nd = 360
    SM['dirs'] = np.linspace(0,359,nd)
    SM['freqs'] = np.linspace(0,3,nf) 
    SM['S'] = np.zeros([nf,nd])
    
    # Resolution of estimation
    if inst == 'cam': EP['nfft'] = 256
    else: EP['nfft'] = 1024
    EP['dres'] = 360
    EP['method'] = 'IMLM'
    
    Options_ = []
    [SMout,EPout] = dirspec(ID,SM,EP,Options_)
    SMout = calc_dir_sprd(SMout)
    return SMout, EPout

# %% Define file locations and path

path = '/Users/cmbaker9/Documents/Research/Lab_Experiments/'
sys.path.append(path + 'codes/pyDIWASP/')
from dirspec import dirspec
# from dirspec import dirspec as dirspec
now = datetime.now()
date = now.strftime("%m-%d-%Y")
figpath = path + 'figures/wave_statistic/directional/' + date + '/'
locpath = path + 'data/processed/insitu/'
# os.mkdir(figpath)

rho = 1000
g = 9.81

# %% Insert Trial information

Tinfo = {
    "Hs": 0.3,
    "Tp": 2,
    "h": 1.07,
    "spread": 40,
    "day": '5',
    "trial": '06',
    "camera": "08-30-2018-2119UTC_Scene1",
    "frames": "07200-11999",
    "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
    "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
    "parray":'sz'
}
datapath = path + 'data/PRJ-1873/inter/random' + Tinfo['day'] + '/Trial' + Tinfo['trial'] + '/'

instruments = ['wg','press']
for i in range(len(instruments)):
    inst = instruments[i]
    [df,XYZ] = read_data(Tinfo,inst,datapath) # read text files
    [SMout,EPout] = calc_SM(Tinfo,inst,df,XYZ)
    exec('SM' + inst + '_' + Tinfo['parray'] + '= SMout')
    del SMout

# %%  surf zone camera timeseries
campath = '/Users/cmbaker9/Documents/Research/Lab_Experiments/data/processed/DEM/TRM-' + Tinfo['camera'] + '/frames_' + Tinfo['frames'] + '/'
fpath = campath + 'pressure_gage_xy.csv'
dpxyz = pd.read_csv(fpath, header=None, names=['x', 'y'])
dpxyz['z'] = np.zeros((12,1))+Tinfo['h']
XYZ = np.transpose(dpxyz.values)

fpath = campath + 'pressure_gage_timeseries.csv'
df = pd.read_csv(fpath, header=None, names=['z1', 'z2', 'z3', 'z4', 'z5', 'z6', 'z7', 'z8', 'z9', 'z10', 'z11', 'z12'])

inst = 'cam'
[SMout,EPout] = calc_SM(Tinfo,inst,df,XYZ)
exec('SM' + inst + '_' + Tinfo['parray'] + '= SMout')

# %% inner shelf trial

Tinfo = {
    "Hs": 0.3,
    "Tp": 2,
    "h": 1.07,
    "spread": 40,
    "day": '9',
    "trial": '04',
    "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
    "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
    "parray":'is'
}
datapath = path + 'data/PRJ-1873/inter/random' + Tinfo['day'] + '/Trial' + Tinfo['trial'] + '/'

instruments = ['wg','press']
for i in range(len(instruments)):
    inst = instruments[i]
    [df,XYZ] = read_data(Tinfo,inst,datapath) # read text files
    [SMout,EPout] = calc_SM(Tinfo,inst,df,XYZ)
    exec('SM' + inst + '_' + Tinfo['parray'] + '= SMout')
    del SMout

# %% Plots

fsmall        = np.arange(1,3,0.01)
fsmall4       = (10**-2.5)*(fsmall**-4)

fig, axs = plt.subplots(2,figsize=(8,7))
axs[0].plot(SMwg_sz['freqs'],SMwg_sz['Sf'], c='k', lw=1, label='wg sz')
axs[0].plot(SMwg_is['freqs'],SMwg_is['Sf'], c='k', lw=1, linestyle='-.', label='wg is')
axs[0].plot(SMpress_sz['freqs'],SMpress_sz['Sf'], c='r', lw=1, label='press sz')
axs[0].plot(SMpress_is['freqs'],SMpress_is['Sf'], c='r', lw=1, linestyle='-.', label='press is')
# axs[0].plot(SMcam_sz['freqs'],SMcam_sz['Sf'], c='b', lw=1, label='cam sz')
axs[0].plot(fsmall,fsmall4, c='g', lw=1, label='$f^{-4}$')
axs[0].set_xlim(SMwg_sz['freqs'].min(),SMwg_sz['freqs'].max())
axs[0].set_yscale('log')
axs[0].set_ylim((10**-4.5),10**-1.5)
axs[0].set_xlabel(r'$f$ (Hz)')
axs[0].set_ylabel(r'$S_{f}$ (m$^{2}$/Hz)')
axs[0].legend(prop={'size': 10})
axs[0].grid(True, alpha = 0.2)

axs[1].plot(SMwg_sz['dirs'],SMwg_sz['Sd'], c='k', lw=1)
axs[1].plot(SMwg_is['dirs'],SMwg_is['Sd'], c='k', lw=1, linestyle='-.')
axs[1].plot(SMpress_sz['dirs'],SMpress_sz['Sd'], c='r', lw=1)
axs[1].plot(SMpress_is['dirs'],SMpress_is['Sd'], c='r', lw=1, linestyle='-.')
# axs[1].plot(SMcam_sz['dirs'],SMcam_sz['Sd'], c='b', lw=1)
axs[1].set_xlim(SMwg_sz['dirs'].min(),SMwg_sz['dirs'].max())
axs[1].set_yscale('log')
axs[1].set_ylim((10**-6.5),10**-2)
axs[1].set_xlabel('Deg. $(^{\circ})$')
axs[1].set_ylabel('$S_{d}$ (m$^{2}$/deg)')
axs[1].grid(True, alpha = 0.2)
plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sf_Sd_nocam.png', bbox_inches='tight')
plt.show()

# %% Plots

fig, axs = plt.subplots(2,figsize=(8,7))
axs[0].plot(SMwg_sz['freqs'],SMwg_sz['th_2'], c='k', lw=1, label='wg sz')
axs[0].plot(SMwg_is['freqs'],SMwg_is['th_2'], c='k', lw=1, linestyle='-.', label='wg is')
axs[0].plot(SMpress_sz['freqs'],SMpress_sz['th_2'], c='r', lw=1, label='press sz')
axs[0].plot(SMpress_is['freqs'],SMpress_is['th_2'], c='r', lw=1, linestyle='-.', label='press is')
# axs[0].plot(SMcam_sz['freqs'],SMcam_sz['th_2'], c='b', lw=1, label='cam sz')
axs[0].set_xlim(SMwg_sz['freqs'].min(),SMwg_sz['freqs'].max())
axs[0].set_ylim(-100,100)
axs[0].set_xlabel('$f$ (Hz)')
axs[0].set_ylabel('$\Theta (^{\circ})$')
axs[0].legend(prop={'size': 10})
axs[0].grid(True, alpha = 0.2)

axs[1].plot(SMwg_sz['freqs'],SMwg_sz['sig_2'], c='k', lw=1) #, label='Data')
axs[1].plot(SMwg_is['freqs'],SMwg_is['sig_2'], c='k', lw=1, linestyle='-.')
axs[1].plot(SMpress_sz['freqs'],SMpress_sz['sig_2'], c='r', lw=1)
axs[1].plot(SMpress_is['freqs'],SMpress_is['sig_2'], c='r', lw=1, linestyle='-.')
# axs[1].plot(SMcam_sz['freqs'],SMcam_sz['sig_2'], c='b', lw=1)
axs[1].set_xlim(SMwg_sz['freqs'].min(),SMwg_sz['freqs'].max())
axs[1].set_ylim(0,45)
axs[1].set_xlabel('$f$ (Hz)')
axs[1].set_ylabel('$\sigma_{\Theta} (^{\circ}$)')
axs[1].grid(True, alpha = 0.2)
plt.rcParams["font.family"] = "serif"
plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_theta_sigma_nocam.png', bbox_inches='tight')
plt.show()

# plt.close('all')