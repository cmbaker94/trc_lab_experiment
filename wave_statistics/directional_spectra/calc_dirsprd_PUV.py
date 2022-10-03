
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
plt.close('all')
import matplotlib as mpl
import pandas as pd
from datetime import datetime
from itertools import combinations
import operator as op
from functools import reduce
import random
import copy
# import seaborn as sns
# sns.set()
from matplotlib import rc
# font = {'family' : 'normal',
#         'weight' : 'normal',
#         'size'   : 18}

# font = {'weight' : 'normal',
#         'size'   : 18}
font = {'size':18}
rc('font', **font)


# rc('text', usetex=True)
# mpl.rcParams['text.usetex'] = True

sys.path.append("E:\\code\\trc_lab_experiment\\toolbox")
from trial_files import trial_files
from calc_SM import calc_SM
from calc_dir_sprd import calc_dir_sprd
from calc_energy_weighted_mean import calc_energy_weighted_mean
from read_insitu_data import read_insitu_data

# %% Add functions

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom  # or / in Python 2

# %% Define file locations and path

path = 'E:/'
now = datetime.now()
date = now.strftime("%m-%d-%Y")
figpath = path + 'figures/wave_statistics/directional/' + date + '/'
locpath = path + 'data/processed/insitu/'
if not os.path.exists(figpath):
    os.mkdir(figpath)

# %% Insert Trial information
spread = 0
Hs = 0.3
h = 1.07
Tp = 2

Tinfo = trial_files(Hs,Tp,spread,h)

conditions =  'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_tide' + str(int(Tinfo['h']*100)) + '_spread' + str(Tinfo['spread'])
datapath = path + 'data/PRJ-1873/inter/random' + Tinfo['day'] + '/Trial' + Tinfo['trial'] + '/'

# %%  surf zone camera timeseries

puvpath = 'E:/data/processed/insitu/' + Tinfo['insitu'] + '/'
dname = 'PUV_p11.csv'
PUV11 = pd.read_csv(puvpath+dname,header=None,names=['freq','SSE','dir2','sig2']) 
dname = 'PUV_p06.csv'
PUV06 = pd.read_csv(puvpath+dname,header=None,names=['freq','SSE','dir2','sig2']) 

PUV11['SSE']=PUV11['SSE'].replace(0, np.nan)
PUV06['SSE']=PUV06['SSE'].replace(0, np.nan)

# %% Create scatter plot
freqrange = [0.1, 1.0]

dirmean = dict()
sprdmean = dict()

# dirmean['PUV06'] = calc_energy_weighted_mean(PUV06['dir2'], PUV06['freq'], freqrange)
# dirmean['PUV11'] = calc_energy_weighted_mean(PUV11['dir2'], PUV11['freq'], freqrange)
# dirmean['wg'] = calc_energy_weighted_mean(SMwg_sz['th_2'], SMwg_sz['freqs'], freqrange)
dirmean['press_sz'] = calc_energy_weighted_mean(SMpress_sz['th_2'], SMpress_sz['freqs'], freqrange)
# dirmean['press_is'] = calc_energy_weighted_mean(SMpress_is['th_2'], SMpress_is['freqs'], freqrange)
# dirmean['cam_sz'] = calc_energy_weighted_mean(SMcam_sz['th_2'], SMcam_sz['freqs'], freqrange)
# dirmean['lid_sz'] = calc_energy_weighted_mean(SMlid_sz['th_2'], SMlid_sz['freqs'], freqrange)


# sprdmean['PUV06']  = calc_energy_weighted_mean(PUV06['sig2']*180/np.pi, PUV06['freq'], freqrange)
# sprdmean['PUV11']  = calc_energy_weighted_mean(PUV11['sig2']*180/np.pi, PUV06['freq'], freqrange)
# sprdmean['wg'] = calc_energy_weighted_mean(SMwg_sz['sig_2'], SMwg_sz['freqs'], freqrange)
sprdmean['press_sz'] = calc_energy_weighted_mean(SMpress_sz['sig_2'], SMpress_sz['freqs'], freqrange)
# sprdmean['press_is'] = calc_energy_weighted_mean(SMpress_is['sig_2'], SMpress_is['freqs'], freqrange)
# sprdmean['cam_sz'] = calc_energy_weighted_mean(SMcam_sz['sig_2'], SMcam_sz['freqs'], freqrange)
# sprdmean['lid_sz'] = calc_energy_weighted_mean(SMlid_sz['sig_2'], SMlid_sz['freqs'], freqrange)

ddir = pd.DataFrame()
ddir = ddir.append(dirmean, ignore_index=True)

dsprd = pd.DataFrame()
dsprd = dsprd.append(sprdmean, ignore_index=True)

save_dir =savepath + 'mean_direction_' + inst + '_' + ffaddition + '.csv'
ddir.to_csv(save_dir, index = False, header=True)

save_sprd = savepath + 'directional_spread_' + inst + '_' + ffaddition + '.csv'
dsprd.to_csv(save_sprd, index = False, header=True)
