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
        - 'layout' provides the x,yz in each column
        - 'datatyles' list of sensor types (ie 'elev' for wave gages)
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

path = '/Users/cmbaker9/Documents/Research/Lab_Experiments/'
sys.path.append(path + 'codes/pyDIWASP/')
from dirspec import dirspec
# from dirspec import dirspec as dirspec
now = datetime.now()
date = now.strftime("%m-%d-%Y")
figpath = path + 'figures/wave_statistic/directional/' + date + '/'
# os.mkdir(figpath)

# %% Insert Trial information

Tinfo = {
    "Hs": 0.3,
    "Tp": 2,
    "tide": 1.07,
    "spread": 40,
    "day": '7',
    "trial": '14',
    "insitu": "09-01-2018-2213UTC",
    "wgnos": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15]
}
datapath = path + 'data/raw/insitu/Trials/random' + Tinfo['day'] + '/Trial' + Tinfo['trial'] + '/'
wgno = Tinfo['wgnos']


# %% Read text files

for i in range(len(wgno)):
    wgname = 'wg' + str(wgno[i])
    wgfile = datapath + wgname + '.txt'
    
    f = open(wgfile, "r")
    lines = f.readlines()
    f.close()
    x = float(lines[65].split()[2])
    y = float(lines[66].split()[2])
    
    if i == 0: 
        temp=np.loadtxt(wgfile, delimiter=',', comments='%', unpack=True)
        df = pd.DataFrame(temp,columns=['wg1'])
        dxy = pd.DataFrame([x, y],columns=['wg1'])
    else:
        df[wgname] = pd.read_csv(wgfile, delimiter=',', comment='%', header=None)
        dxy[wgname] = [x, y]
        
dXY = dxy.values

# %% Set variables and run dirspec function

SM = dict()
EP = dict()
ID = dict()

# Instrument data
ID['datatypes'] = ['elev']
for i in range(0,13):
    ID['datatypes'].append('elev')
ID['depth'] = 1.07
ID['data'] = df.values
N = dxy.shape[1]
ID['layout'] = np.concatenate((np.array([np.transpose(dXY[0]), np.transpose(dXY[1])]), np.zeros([1,N])), axis=0)
ID['fs'] = 100

# Spectral matrix
nf = 3001
nd = 360
SM['dirs'] = np.linspace(0,359,nd)
SM['freqs'] = np.linspace(0,3,nf) 
SM['S'] = np.zeros([nf,nd])

# Resolution of estimation
EP['nfft'] = 1024
EP['dres'] = 360
EP['method'] = 'IMLM'

Options_ = []
[Smout,EPout] = dirspec(ID,SM,EP,Options_)

# save as csv

# %% Compute theta and sigma

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

# %% Plots

fsmall        = np.arange(1,3,0.01)
fsmall4       = (10**-2.5)*(fsmall**-4)

fig, axs = plt.subplots(2,figsize=(8,7))
axs[0].plot(Smout['freqs'],Sf, c='k', lw=1) #, label='Data')
axs[0].plot(fsmall,fsmall4, c='b', lw=1)
axs[0].set_xlim(Smout['freqs'].min(),Smout['freqs'].max())
axs[0].set_yscale('log')
axs[0].set_ylim((10**-4.5),10**-1.5)
axs[0].set_xlabel(r'$f$ (Hz)')
axs[0].set_ylabel(r'$S_{f}$ (m$^{2}$/Hz)')
axs[0].legend(prop={'size': 12})
axs[0].grid(True, alpha = 0.2)

axs[1].plot(Smout['dirs'],Sd, c='k', lw=1) #, label='Data')
axs[1].set_xlim(Smout['dirs'].min(),Smout['dirs'].max())
axs[1].set_yscale('log')
axs[1].set_ylim((10**-6.5),10**-2)
axs[1].set_xlabel('Deg. $(^{\circ})$')
axs[1].set_ylabel('$S_{d}$ (m$^{2}$/deg)')
axs[1].legend(prop={'size': 12})
axs[1].grid(True, alpha = 0.2)
plt.savefig(figpath + Tinfo['insitu'] + '_Sf_Sd.png', bbox_inches='tight')
plt.show()

fig, axs = plt.subplots(2,figsize=(8,7))
axs[0].plot(freq,th_2, c='k', lw=1) #, label='Data')
axs[0].set_xlim(Smout['freqs'].min(),Smout['freqs'].max())
axs[0].set_ylim(-100,100)
axs[0].set_xlabel('$f$ (Hz)')
axs[0].set_ylabel('$\Theta (^{\circ})$')
axs[0].legend(prop={'size': 12})
axs[0].grid(True, alpha = 0.2)

axs[1].plot(freq,sig_2, c='k', lw=1) #, label='Data')
axs[1].set_xlim(Smout['freqs'].min(),Smout['freqs'].max())
axs[1].set_ylim(0,45)
axs[1].set_xlabel('$f$ (Hz)')
axs[1].set_ylabel('$\sigma_{\Theta} (^{\circ}$)')
axs[1].legend(prop={'size': 12})
axs[1].grid(True, alpha = 0.2)
plt.rcParams["font.family"] = "serif"
plt.savefig(figpath + Tinfo['insitu'] + '_theta_sigma.png', bbox_inches='tight')
plt.show()

plt.close('all')