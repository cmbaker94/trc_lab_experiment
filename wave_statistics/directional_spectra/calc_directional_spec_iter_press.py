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
Tspreads    = [30];
THs = [0.25];
TTp = [2];
Ttide = [1.07];

for ispreads in range(len(Tspreads)):
    # %% Insert Trial information
    spread = Tspreads[ispreads]
    Hs = THs[ispreads]
    h = Ttide[ispreads]
    Tp = TTp[ispreads]

    Tinfo = trial_files(Hs,Tp,spread,h)
    
    conditions =  'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_tide' + str(int(Tinfo['h']*100)) + '_spread' + str(Tinfo['spread'])
    datapath = path + 'data/PRJ-1873/inter/random' + Tinfo['day'] + '/Trial' + Tinfo['trial'] + '/'
    
    
    # %%  surf zone camera timeseries
    
    # # campath = '/Users/cmbaker9/Documents/Research/Lab_Experiments/data/processed/cameras/TRM-' + Tinfo['camera'] + '/frames_' + Tinfo['frames'] + '/'
    # campath = 'E:/data/processed/conditions/' + conditions + '/' + Tinfo['clpath'] + '/time_' + Tinfo['timesection'] + '/'
    # fpath = campath + 'cam_grid_xy_dx50_dy100.csv' #_5cm
    # dpxyz = pd.read_csv(fpath, header='infer')#, names=['x', 'y'])
    
    # dpid = dpxyz.transpose()
    # dpid.columns = ['x','y']
    
    instnum = 5
    itnummax = 200
    itnum = ncr(12,instnum)
    inst = 'press'
    
    # Get all combinations of [1, 2, 3]
    # and length 2
    comb = combinations([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], instnum)
    
    # Print the obtained combinations
    count = 0
    sinst = np.empty([itnum, instnum])
    for i in list(comb):
        sinst[count,:] = i 
        count = count+1
    
    np.random.shuffle(sinst)
    
    if itnum > itnummax:
        itnum = itnummax
    
    conditions =  'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_tide' + str(int(Tinfo['h']*100)) + '_spread' + str(Tinfo['spread'])
    datapath = 'G:/PRJ-1873/inter/random' + Tinfo['day'] + '/Trial' + Tinfo['trial'] + '/'
    savepath = 'E:/data/processed/conditions/' + conditions + '/' + Tinfo['clpath'] + '/time_' + Tinfo['timesection'] + '/'
    ffaddition = 'press_sz_' + str(instnum) + 'inst_' + str(itnum) + 'iter'
      
    Ssum = np.empty([301,360,itnum])
    for i in range(itnum):
        ipick = sinst[i,:].astype(int)
        [df,XYZ] = read_insitu_data(Tinfo,inst,datapath) # read text files
        df2 = df.iloc[:, ipick]
        XYZ2 = XYZ[:, ipick]
        [SMout,EPout] = calc_SM(Tinfo,inst,df2,XYZ2)
        # [SMout,EPout] = calc_SM(Tinfo,inst,df,XYZ)
        Ssum[:,:,i] = SMout['S'].real
        # SMout['S'] = Ssum/itnum
        
    SMavg = dict()
    Savg = Ssum.mean(axis=2)
    SMavg['S'] = Savg
    SMavg['dirs'] = SMout['dirs']
    SMavg['freqs'] = SMout['freqs']
    SMavg = calc_dir_sprd(SMavg)
    exec('SM' + inst + '_' + Tinfo['parray'] + '= copy.deepcopy(SMavg)')
    
    # exec('SMsave = SM' + inst + '_' + Tinfo['parray'])
    savepath_S = savepath + 'spec_' + inst + '_' + ffaddition + '.csv'
    savepath_f = savepath + 'freq_' + inst + '_' + ffaddition + '.csv'
    savepath_d = savepath + 'dirs_' + inst + '_' + ffaddition + '.csv'
    np.savetxt(savepath_f, SMpress_sz['freqs'], delimiter=",")
    np.savetxt(savepath_d, SMpress_sz['dirs'], delimiter=",")
    np.savetxt(savepath_S, SMpress_sz['S'].real, delimiter=",")
    
    # %% Plots
    
    fsmall        = np.arange(1,3,0.01)
    fsmall4       = (10**-2.5)*(fsmall**-4)
    
    fig, axs = plt.subplots(2,figsize=(8,7))
    axs[0].plot(SMpress_sz['freqs'],SMpress_sz['Sf'], c='r', lw=1, label='press sz')
    axs[0].plot(fsmall,fsmall4, c='m', lw=1, label=r'$f^{-4}$')
    axs[0].set_xlim(SMpress_sz['freqs'].min(),2.5)
    axs[0].set_yscale('log')
    axs[0].set_ylim((10**-4.5),10**-1.5)
    axs[0].legend(prop={'size': 10})
    axs[0].grid(True, alpha = 0.2)
    
    axs[1].plot(SMpress_sz['dirs'],SMpress_sz['Sd'], c='r', lw=1)
    axs[1].set_yscale('log')
    axs[1].set_ylim((10**-6.5),10**-2)
    axs[1].set_xlabel(r'Deg. $(^{\circ})$')
    axs[1].set_ylabel(r'$S_{d}$ (m$^{2}$/deg)')
    axs[1].grid(True, alpha = 0.2)
    plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sf_Sd_cam_' + ffaddition +'.png', bbox_inches='tight')
    plt.show()
    
    # %% Plots
    
    fig, axs = plt.subplots(2,figsize=(8,7))
    axs[0].plot(SMpress_sz['freqs'],SMpress_sz['th_2'], c='r', lw=1, label='press sz')
    axs[0].set_xlim(SMpress_sz['freqs'].min(),2.5)
    axs[0].set_ylim(-100,100)
    axs[0].set_xlabel(r'$f$ (Hz)')
    axs[0].set_ylabel(r'$\Theta (^{\circ})$')
    axs[0].legend(prop={'size': 10})
    axs[0].grid(True, alpha = 0.2)
    
    axs[1].plot(SMpress_sz['freqs'],SMpress_sz['sig_2'], c='r', lw=1)
    axs[1].set_xlim(SMpress_sz['freqs'].min(),2.5)
    axs[1].set_ylim(0,45)
    axs[1].set_xlabel(r'$f$ (Hz)')
    axs[1].set_ylabel(r'$\sigma_{\Theta} (^{\circ}$)')
    axs[1].grid(True, alpha = 0.2)
    # plt.rcParams["font.family"] = "serif"
    plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_theta_sigma_cam_' + ffaddition +'.png', bbox_inches='tight')
    plt.show()
    
    
    # %% log scale
    
    from matplotlib import colors, ticker, cm
    from matplotlib.colors import LogNorm
    cax = [10**-6, 10**-4.65]
    #cax = [1e-6, 1e-5]
    
    # values = np.transpose(SMcam_sz['S'].real)
    #test = SMcam_sz['S'].real
    values = np.transpose(copy.copy(Savg))
    values[values<cax[0]]=cax[0]
    values[values>cax[1]]=cax[1]
    
    theta = np.array(SMpress_sz['dirs'])
    azimuths = np.radians(SMpress_sz['dirs'])
    zeniths = np.array(SMpress_sz['freqs'][0:141])
     
    values = np.array(values[:, 0:141])
    values = values.reshape(len(azimuths), len(zeniths))
     
    r, theta = np.meshgrid(zeniths, azimuths)
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'),figsize=(12, 9), dpi=80)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    # CS = plt.contourf(theta, r, values, levels=[10**-5, 10**-4, 10**-3, 10**-2],  cmap='plasma')
    
    
    cnorm= LogNorm(vmin=cax[0], vmax=cax[1])
    
    
    # CS = plt.contourf(theta, r, values, locator=ticker.LogLocator(subs=range(1,10)),  cmap='autumn', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
    CS = plt.contourf(theta, r, values, levels=np.logspace(-6,-4.4,31),  cmap='viridis', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
    #CS = plt.contourf(theta, r, values, levels=np.logspace(-6,-5,31),  cmap='viridis', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
    cbaxes = fig.add_axes([0.86, 0.1, 0.03, 0.7])
    #cbar = fig.colorbar(CS, ticks=[1e-6, 1e-5, 1e-4], cax = cbaxes)
    cbar = fig.colorbar(CS, ticks=[1e-6, 1e-5], cax = cbaxes)
    #cbar.ax.set_yticklabels([r'$10^{-6}$', '$10^{-5}$', '$10^{-4}$'])
    cbar.ax.set_yticklabels([r'$10^{-6}$', '$10^{-5}$'])
    cbar.set_label(r'$S(f,\Theta)$',labelpad=-40, y=1.1, rotation=0)
    
    
    rlabels = ax.get_ymajorticklabels()
    for label in rlabels[:-1]:
        label.set_color('white')
    
    plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sfd_polar_log_' + ffaddition +'.png', bbox_inches='tight')
    
    # plt.close('all')
    
    # %% camera no log
    
    cax = [10**-6, 10**-4.55]
    cax = [0,0.000022]
    # values = np.transpose(SMcam_sz['S'].real)
    values = np.transpose(copy.copy(Savg))
    values[values<cax[0]]=cax[0]
    values[values>cax[1]]=cax[1]
    
    theta = np.array(SMpress_sz['dirs'])
    azimuths = np.radians(SMpress_sz['dirs'])
    zeniths = np.array(SMpress_sz['freqs'][0:121])
     
    values = np.array(values[:, 0:121])
    values = values.reshape(len(azimuths), len(zeniths))
     
    r, theta = np.meshgrid(zeniths, azimuths)
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'),figsize=(12, 9), dpi=80)
    plt.rcParams['font.size'] = '26'
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    # ax.set_rlabel_position(-5)
    # CS = plt.contourf(theta, r, values, levels=[10**-5, 10**-4, 10**-3, 10**-2],  cmap='plasma')
    
    ax.xaxis.set_tick_params(pad=20) # pad the theta away from the plot
    # cnorm= LogNorm(vmin=cax[0], vmax=cax[1])
    
    
    # CS = plt.contourf(theta, r, values, locator=ticker.LogLocator(subs=range(1,10)),  cmap='autumn', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
    # CS = plt.contourf(theta, r, values, levels=np.logspace(-6,-4,31),  cmap='viridis', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
    CS = plt.contourf(theta, r, values, levels=np.linspace(cax[0],cax[1],41), cmap='viridis')#, vmin=0.0001, vmax=0.01, cmap='plasma')
    
    # CS = plt.contourf(theta, r, values, levels=np.logspace(-6,-5,31),  cmap='viridis', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
    cbaxes = fig.add_axes([0.86, 0.14, 0.03, 0.6])
    # cbar = fig.colorbar(CS, ticks=[1e-6, 1e-5, 1e-4], cax = cbaxes)
    # cbar = fig.colorbar(CS, ticks=[0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.00007, 0.00008, 0.00009, 0.0001], cax = cbaxes)
    cbar = fig.colorbar(CS, cax = cbaxes)
    
    ax.set_rticks([0, 0.3, 0.6, 0.9, 1.2])
    
    # thetaticks = np.arange(0,360,45)
    # ax.set_thetagrids(thetaticks, frac=1.3)
    # cbar.ax.set_yticklabels([r'$10^{-6}$', '$10^{-5}$', '$10^{-4}$'])
    # cbar.ax.set_yticklabels([r'$10^{-6}$', '$10^{-5}$'])
    cbar.set_label(r'$S(f,\Theta)$',labelpad=-50, y=1.2, rotation=0, fontsize=28)
    # sns.set_style("whitegrid")
    
    rlabels = ax.get_ymajorticklabels()
    # rlabels.set_fontsize(28)
    # ax.set_xmajorticklabels(fontsize=22)
    for label in rlabels[:-1]:
        label.set_color('white')
        label.set_fontsize(28)
        # label.set_position(1)
        
    
    plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sfd_polar_nolog_' + ffaddition +'_clim1.png', bbox_inches='tight')
    
    
    # # plt.close('all')
    
    
    
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
