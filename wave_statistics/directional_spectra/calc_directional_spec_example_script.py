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
import random
import copy
from matplotlib import rc
font = {'size':18}
rc('font', **font)

sys.path.append("E:\\code\\trc_lab_experiment\\toolbox")
from trial_files import trial_files
from calc_SM import calc_SM
from calc_dir_sprd import calc_dir_sprd
from calc_energy_weighted_mean import calc_energy_weighted_mean

# %% Define file locations and path

path = 'E:/'
now = datetime.now()
date = now.strftime("%m-%d-%Y")
figpath = path + 'figures/wave_statistics/directional/' + date + '/'
locpath = path + 'data/processed/insitu/'
if not os.path.exists(figpath):
    os.mkdir(figpath)

rho = 1000
g = 9.81

# %% Insert Trial information and analysis information

# trial conditions
spread    = [30];
Hs = [0.25];
Tp = [2];
h = [1.07];

# analysis information
inst = 'cam' # instrument type
pick = False # if true: code selects locations to compute spectra, if false: code randomly selects locations
if pick == False:
    itnum = 200 # number of iterations
    instnum = 14 # number of instruments selected each iteration


# %% Start computing directional sprectra

    
Tinfo = trial_files(Hs,Tp,spread,h) # get trial information based on conditions
conditions =  'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_tide' + str(int(Tinfo['h']*100)) + '_spread' + str(Tinfo['spread'])
datapath = path + 'data/PRJ-1873/inter/random' + Tinfo['day'] + '/Trial' + Tinfo['trial'] + '/'

# load stereo sea-surface elevation gridded data
campath = 'E:/data/processed/conditions/' + conditions + '/' + Tinfo['clpath'] + '/time_' + Tinfo['timesection'] + '/'
fprocinfo = 'regx28-33_regy-10-4_dx50_dy50'
# load xy position
fpath = campath + 'cam_grid_xy_' + fprocinfo + '.csv'
dpxyz = pd.read_csv(fpath, header='infer')
dpid = dpxyz.transpose()
dpid.columns = ['x','y']
# load data
fpath = campath + 'cam_grid_timeseries_' + fprocinfo + '.csv'
dfall = pd.read_csv(fpath, header='infer')
dfall = dfall.interpolate(method='linear', limit_direction='both')

if pick == True:
    locxy = [30, 0],[30,-1],[30,-3],[29.5, 0],[29.5,-1],[29.5,-3],[28, 0],[28,-1],[28,-3]
    ffaddition = 'grid_pickloc' 
    loc2calc = np.array([])
    for idxy in locxy:
        print(idxy)
        ide = (dpid == idxy).all(1)
        print(np.where(ide)[0])
        loc2calc = np.append(loc2calc,np.where(ide)[0])
    dploc = dpxyz.iloc[:,loc2calc]
    XYZ = dploc.values
    z = np.zeros((1,XYZ.shape[1]))+Tinfo['h']
    XYZ = np.vstack([XYZ, z])
    df = dfall.iloc[:,loc2calc]
    [SMout,EPout] = calc_SM(Tinfo,inst,df,XYZ)
    exec('SM' + inst + '_' + Tinfo['parray'] + '= copy.deepcopy(SMout)')
elif pick == False:
    ffaddition = fprocinfo + '_grid_rand_' + str(instnum) + 'inst_' + str(itnum) + 'iter'
    Ssum = np.empty([301,360,itnum])
    stdxt = np.empty([itnum])
    stdyt = np.empty([itnum])
    for i in range(itnum):
        while True:
            loc2calc = random.sample(range(0, dpxyz.shape[1]), instnum)
            dploc = dpxyz.iloc[:,loc2calc]
            XYZ = dploc.values
            z = np.zeros((1,XYZ.shape[1]))+Tinfo['h']
            XYZ = np.vstack([XYZ, z])
            stdx = np.std(XYZ[0,:])
            stdy = np.std(XYZ[1,:])
            stdxt[i] = np.std(XYZ[0,:])
            stdyt[i] = np.std(XYZ[1,:])
            # code will select a new set of locations if less than or more than a specified standard deviation
            if stdx > 1.3 and stdx < 2 and stdy > 3 and stdy < 7:
                break
        df = dfall.iloc[:,loc2calc]
        [SMout,EPout] = calc_SM(Tinfo,inst,df,XYZ)
        Ssum[:,:,i] = SMout['S'].real
    SMavg = dict()
    Savg = Ssum.mean(axis=2)
    SMavg['S'] = Savg
    SMavg['dirs'] = SMout['dirs']
    SMavg['freqs'] = SMout['freqs']
    SMavg = calc_dir_sprd(SMavg)
    exec('SM' + inst + '_' + Tinfo['parray'] + '= copy.deepcopy(SMavg)')

savepath_S = campath + 'spec_' + inst + '_' + ffaddition + '.csv'
savepath_f = campath + 'freq_' + inst + '_' + ffaddition + '.csv'
savepath_d = campath + 'dirs_' + inst + '_' + ffaddition + '.csv'
np.savetxt(savepath_f, SMcam_sz['freqs'], delimiter=",")
np.savetxt(savepath_d, SMcam_sz['dirs'], delimiter=",")
np.savetxt(savepath_S, SMcam_sz['S'], delimiter=",")

# %% Compute directional spread and spectra from directional spectra
freqrange = [0.1, 1.0]

dirmean = dict()
sprdmean = dict()

dirmean = calc_energy_weighted_mean(SMcam_sz['th_2'], SMcam_sz['freqs'], freqrange)
sprdmean = calc_energy_weighted_mean(SMcam_sz['sig_2'], SMcam_sz['freqs'], freqrange)

# %% Plots

fsmall        = np.arange(1,3,0.01)
fsmall4       = (10**-2.5)*(fsmall**-4)

fig, axs = plt.subplots(2,figsize=(8,7))
axs[0].plot(SMcam_sz['freqs'],SMcam_sz['Sf'], c='b', lw=1, label='cam sz')
axs[0].plot(fsmall,fsmall4, c='m', lw=1, label=r'$f^{-4}$')
axs[0].set_xlim(SMcam_sz['freqs'].min(),2.5)
axs[0].set_yscale('log')
axs[0].set_ylim((10**-4.5),10**-1.5)
axs[0].legend(prop={'size': 10})
axs[0].grid(True, alpha = 0.2)

axs[1].plot(SMcam_sz['dirs'],SMcam_sz['Sd'], c='b', lw=1)
axs[1].set_yscale('log')
axs[1].set_ylim((10**-6.5),10**-2)
axs[1].set_xlabel(r'Deg. $(^{\circ})$')
axs[1].set_ylabel(r'$S_{d}$ (m$^{2}$/deg)')
axs[1].grid(True, alpha = 0.2)
plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sf_Sd_cam_' + ffaddition +'.png', bbox_inches='tight')
plt.show()

# %% Plots

fig, axs = plt.subplots(2,figsize=(8,7))
axs[0].plot(SMcam_sz['freqs'],SMcam_sz['th_2'], c='b', lw=1, label='cam sz')
axs[0].set_xlim(SMcam_sz['freqs'].min(),2.5)
axs[0].set_ylim(-100,100)
axs[0].set_xlabel(r'$f$ (Hz)')
axs[0].set_ylabel(r'$\Theta (^{\circ})$')
axs[0].legend(prop={'size': 10})
axs[0].grid(True, alpha = 0.2)

axs[1].plot(SMcam_sz['freqs'],SMcam_sz['sig_2'], c='b', lw=1)
axs[1].set_xlim(SMcam_sz['freqs'].min(),2.5)
axs[1].set_ylim(0,45)
axs[1].set_xlabel(r'$f$ (Hz)')
axs[1].set_ylabel(r'$\sigma_{\Theta} (^{\circ}$)')
axs[1].grid(True, alpha = 0.2)
plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_theta_sigma_cam_' + ffaddition +'.png', bbox_inches='tight')
plt.show()

# %% camera polar plot log scale

from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm
cax = [10**-6, 10**-4.65]

values = np.transpose(copy.copy(Savg))
values[values<cax[0]]=cax[0]
values[values>cax[1]]=cax[1]

theta = np.array(SMcam_sz['dirs'])
azimuths = np.radians(SMcam_sz['dirs'])
zeniths = np.array(SMcam_sz['freqs'][0:141])
 
values = np.array(values[:, 0:141])
values = values.reshape(len(azimuths), len(zeniths))
 
r, theta = np.meshgrid(zeniths, azimuths)
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'),figsize=(12, 9), dpi=80)
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

cnorm= LogNorm(vmin=cax[0], vmax=cax[1])
CS = plt.contourf(theta, r, values, levels=np.logspace(-6,-4.4,31),  cmap='viridis', norm=cnorm)
cbaxes = fig.add_axes([0.86, 0.1, 0.03, 0.7])
cbar = fig.colorbar(CS, ticks=[1e-6, 1e-5], cax = cbaxes)
cbar.ax.set_yticklabels([r'$10^{-6}$', '$10^{-5}$'])
cbar.set_label(r'$S(f,\Theta)$',labelpad=-40, y=1.1, rotation=0)

rlabels = ax.get_ymajorticklabels()
for label in rlabels[:-1]:
    label.set_color('white')

plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sfd_polar_camera_sz_log_' + ffaddition +'.png', bbox_inches='tight')