#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 14:32:01 2020

@author: cmbaker9
"""
## Sync in situ, lidar, and cameras
# This codes will use a time-lagged cross-correlation analysis to find the time offset between the in situ gages (pressure gage used here), Velodyne lidar, and camera.

# Author: C.M. Baker \
# Last updated: November 23, 2020
from IPython import get_ipython
get_ipython().magic('reset -sf')

import numpy as np
import os
from scipy import signal
import scipy as sp
from scipy import signal
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
from datetime import datetime
procpath = '/Users/cmbaker9/Documents/Research/Lab_Experiments/data/processed/conditions/'
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

rc('font', **font)
rc('text', usetex=True)
plt.close('all')

def dispersi(omega,h):
    """ Solution of the linear dispersion relationship for surface gravity waves.
    Usage: kh = dispersi(const)
    Input: omega = radian frequency, h = water depth
    Output: k = wavenumber
    Given: g = acceleration due to gravity
    The solution satisfies const=kh*tanh(kh) within 1e-06 via Newton-Raphson iteration.  
    The routine returns kh = NaN for negative values of const."""
    g = 9.81
    const = (omega**2)*h/g
    # initialize kh
    kh = np.empty(omega.shape)
    kh[:] = np.nan
    # solution for zero values of const
    m = np.where(const==0)[0]
    if not m:
        kh[m] = np.zeros(len(m))
    # initial value for positive const
    m = np.where(const>0)[0]
    kh[m] = np.sqrt(const[m])
    m = np.where(const>1)[0]
    kh[m] = const[m]
    # iterative solution for positive const
    m = np.where(const>0)[0]
    while np.absolute(const[m] - kh[m]*np.tanh(kh[m])).max() > 10**(-6):
        f = kh[m]*np.tanh(kh[m]) - const[m]
        fprime = kh[m]/np.cosh(kh[m])**2 + np.tanh(kh[m])
        kh[m] = kh[m] - f/fprime
    k = kh/h
    return k

def timeseries_Kp_correction(press,fs,maxfac,offset):
    g = 9.81
    rho = 1000
    etaKp = (press/(rho*g))
    hp = np.nanmean(etaKp)
    h = hp + offset # Total depth
    etaKp = etaKp - hp # subtract MSSE to get etaKp
    temp_fft = sp.fftpack.fft(etaKp)
    temp_psd = np.abs(temp_fft) ** 2
    fftfreq = sp.fftpack.fftfreq(len(temp_psd), 1. / fs)

    omega = np.array(2*np.pi*np.abs(fftfreq))
    k = dispersi(omega,h) # find wavenumbers
    Kp = np.square(np.cosh(k*(h-hp))/np.cosh(k*h)) # compute pressure response factor

    # add in maxfactor
    temp_fft_filt = temp_fft.copy()
    See = temp_fft_filt/Kp
    See[np.abs(fftfreq) > maxfac] = 0
    # nyq = 0.5 * fs
    # lowcut = maxfac / nyq
    # order = 2
    # sos = sp.signal.butter(order, lowcut, analog=False, btype='low', output='sos')
    # w, h = sp.signal.sosfreqz(sos, worN=len(press))
    # See = See*h.real
    # See[np.abs(fftfreq) > maxfac] = 0
    eta_inst = np.real(sp.fftpack.ifft(See))
    return eta_inst

def crosscorr(datax, datay, lag=0, wrap=False):
    """ Lag-N cross correlation. 
    Shifted data filled with NaNs 
    
    Parameters
    ----------
    lag : int, default 0
    datax, datay : pandas.Series objects of equal length
    Returns
    ----------
    crosscorr : float
    """
    if wrap:
        shiftedy = datay.shift(lag)
        shiftedy.iloc[:lag] = datay.iloc[-lag:].values
        return datax.corr(shiftedy)
    else: 
        return datax.corr(datay.shift(lag))


" --------- "
Tinfo = dict()
# Wave conditions
Tinfo['Hs'] = 0.30
Tinfo['Tp'] = 2
Tinfo['watlev'] = 1.07
Tinfo['spread'] = 40

if Tinfo['spread'] == 0:
    Tinfo['timerange'] = '2229-2238'
    Tinfo['trialstart'] = '09-01-2018-2214UTC'
elif Tinfo['spread'] == 20:
    Tinfo['timerange'] = '2237-2246'
    Tinfo['trialstart'] = '08-30-2018-2222UTC'
elif Tinfo['spread'] == 40:
    Tinfo['timerange'] = '2144-2153'
    Tinfo['trialstart'] = '08-30-2018-2129UTC'

# Create compiled name
temp1 = str(Tinfo['Hs']*100).rstrip('0').rstrip('.')
temp2 = str(Tinfo['watlev']*100).rstrip('0').rstrip('.')
Tinfo['comp'] = 'Hs' + temp1 + '_Tp' + str(Tinfo['Tp']) +'_tide'+ temp2 + '_spread' + str(Tinfo['spread'])

# Generate string where data is stored
datafolder = procpath + Tinfo['comp'] + '/' + Tinfo['trialstart'] + '/time_' + Tinfo['timerange'] + '/'

# os.makedirs('/Users/cmbaker9/Documents/Research/Lab_Experiments/figures/conditions/' + Tinfo['comp'] + '/' + Tinfo['trialstart']+ '/time_' + Tinfo['timerange'])
now = datetime.now()
date = now.strftime("%m-%d-%Y")
figfolder = '/Users/cmbaker9/Documents/Research/Lab_Experiments/figures/conditions/' + Tinfo['comp'] + '/' + Tinfo['trialstart']+ '/time_' + Tinfo['timerange'] + '/' + date + '/'
if not os.path.exists(figfolder):
    os.makedirs(figfolder)

# Read csv files
xyz = pd.read_csv(datafolder + 'pressure_array_xyz.csv')
dp = pd.read_csv(datafolder + 'pressure_array_timeseries_full.csv')
dc = pd.read_csv(datafolder + 'cam_array_timeseries.csv')
dl = pd.read_csv(datafolder + 'lidar_array_timeseries_full.csv')
# !head -n 15 $

# In situ data calculations
fs = 100
maxfac = 1.2
offset = 0.05

# Depth-attenuation correction
di = pd.DataFrame()
for (columnName, columnData) in dp.iteritems():
    print('Colunm Name : ', columnName)
    presstemp = columnData.values
    di[columnName] = timeseries_Kp_correction(presstemp,fs,maxfac,offset)
    
# Hydrostatic Approximation
g = 9.81
rho = 1000
dpe = (dp/(rho*g))+offset

# Demean dataframes
dpe = dpe - dpe.mean()
dc = dc-dc.mean()
dl = dl-dl.mean()

# Define time and Interopolate to 8 Hz:
# pressure gages at 100 Hz, cameras at 8 Hz, lidar at 10 Hz
# dpe['date'] = pd.date_range('9/1/2018 22:14:00', periods=len(dpe), freq='10ms')
# dpes = dpe.set_index('date').resample('125ms').mean()
# di['date'] = pd.date_range('9/1/2018 22:14:00', periods=len(dpe), freq='10ms')
# dis = di.set_index('date').resample('125ms').mean()

# dc['date'] = pd.date_range('9/1/2018 22:14:00', periods=len(dc), freq='125ms')
# dcs = dc.set_index('date').resample('125ms').mean()
# dl['date'] = pd.date_range('9/1/2018 22:14:00', periods=len(dl), freq='100ms')
# dls = dl.set_index('date').resample('125ms').mean()

dpe['date'] = pd.date_range('00:00', periods=len(dpe), freq='10ms')
dpes = dpe.set_index('date').resample('125ms').mean()
di['date'] = pd.date_range('00:00', periods=len(dpe), freq='10ms')
dis = di.set_index('date').resample('125ms').mean()

dc['date'] = pd.date_range('00:00', periods=len(dc), freq='125ms')
dcs = dc.set_index('date').resample('125ms').mean()
dl['date'] = pd.date_range('00:00', periods=len(dl), freq='100ms')
dls = dl.set_index('date').resample('125ms').mean()

# plot raw data
fig, axs = plt.subplots(3, figsize=(15,10))
axs[0].plot(dpe['date'], dpe['p11'], c='k', lw=0.3, label='press')
# axs[0].plot(dpes['date'], dpes['p11'], c='r', lw=0.2, label='press resamp')
axs[1].plot(dc['date'], dc['c11'], c='k', lw=0.3, label='camera')
axs[2].plot(dl['date'], dl['l11'], c='k', lw=0.3, label='lidar')
axs[2].plot(dc['date'], dls['l11'], c='r', lw=0.2, label='press resamp')
# axs[0].set_ylim((-.3,.3))
axs[0].set_title('Pressure Timeseries')
axs[1].set_title('Stereo Timeseries')
axs[2].set_title('Lidar Timeseries')
for ax in fig.get_axes():
    ax.set_xlim((dpe['date'].iloc[0],dpe['date'].iloc[-1]))
#     ax.legend()
    ax.grid(True, alpha = 0.3)
    ax.set_xlabel('t (s)')
    ax.set_ylabel('$\eta$ (m)')
# plt.savefig(figfolder + 'timeseries.png', bbox_inches='tight')
plt.show()

# %% Find time lag
tpes = dpes['p11'].reset_index(drop = True) # same result for pressure and depth-attenuation corrected
tls = dls['l11'].reset_index(drop = True) # referenced to lidar since that is started first
tc = dc['c11'].reset_index(drop = True)

seconds = 60*1
fps = 8

rs_pc = [crosscorr(tpes,tc, lag) for lag in range(0,len(tpes))]
rs_pl = [crosscorr(tpes,tls, lag) for lag in range(0,len(tpes))]
rs_cp = [crosscorr(tc,tpes, lag) for lag in range(0,len(tpes))]
rs_cl = [crosscorr(tc,tls, lag) for lag in range(0,len(tls))]
rs_lp = [crosscorr(tls,tpes, lag) for lag in range(0,len(tls))]
rs_lc = [crosscorr(tls,tc, lag) for lag in range(0,len(tls))]

rs = [(np.nan, np.argmax(rs_pc), np.argmax(rs_pl)) , 
      (np.argmax(rs_cp), np.nan, np.argmax(rs_cl)), 
      (np.argmax(rs_lp), np.argmax(rs_lc), np.nan)]

drs = pd.DataFrame(rs, columns = ['press', 'camera', 'lidar'], index=['press', 'camera', 'lidar'])

# %%plot examples
f,ax=plt.subplots(figsize=(14,3))
ax.plot(rs_lp)
ax.axvline(np.argmax(rs_lp),color='r',linestyle='--',label='Peak synchrony')
plt.legend()
plt.show()

# %%

# Define the index that offsetting timeseries by
if Tinfo['spread'] == 0:
    last = 'l'
elif Tinfo['spread'] == 20:
    # lidar before pressure gage
    # camera??
    last = 'l'
elif Tinfo['spread'] == 40:
    last = 'l'
    # lidar before pressure gage
    # camera??

# New timerseries
if last == 'p':
    cstart = int(drs.loc['camera']['press']) # plus 5 to adjust for the middle
    lstart = int(drs.loc['lidar']['press'])
    dc = dc.truncate(before=cstart)
    dl = dl.truncate(before=lstart)
    dc['date'] = pd.date_range('00:00', periods=len(dc), freq='125ms')
    dl['date'] = pd.date_range('00:00', periods=len(dl), freq='100ms')
elif last == 'c':
    pstart = int((drs.loc['press']['camera']/8)*100)+5 # plus 5 to adjust for the middle
    lstart = int(drs.loc['lidar']['camera'])
    di = di.truncate(before=pstart)
    dpe = dpe.truncate(before=pstart)
    dl = dl.truncate(before=lstart)
    dpe['date'] = pd.date_range('00:00', periods=len(dpe), freq='10ms')
    di['date'] = pd.date_range('00:00', periods=len(dpe), freq='10ms')
    dl['date'] = pd.date_range('00:00', periods=len(dl), freq='100ms')
elif last == 'l':
    pstart = int((drs.loc['press']['lidar']/8)*100)+5 # plus 5 to adjust for the middle
    cstart = int(drs.loc['camera']['lidar'])
    di = di.truncate(before=pstart)
    dpe = dpe.truncate(before=pstart)
    dc = dc.truncate(before=cstart)
    dpe['date'] = pd.date_range('00:00', periods=len(dpe), freq='10ms')
    di['date'] = pd.date_range('00:00', periods=len(dpe), freq='10ms')
    dc['date'] = pd.date_range('00:00', periods=len(dc), freq='125ms')
elif last == 'u':
    pstart = int((drs.loc['press']['lidar']/8)*100)+5 # plus 5 to adjust for the middle
    di = di.truncate(before=pstart)
    dpe = dpe.truncate(before=pstart)
    dpe['date'] = pd.date_range('00:00', periods=len(dpe), freq='10ms')

    # print('pressure gages are '+ str(p2l) + ' timesteps ahead of lidar')
    # print(dpes.index[int(drs.loc['press']['lidar'])])
    # print(dpe['date'].iloc[p2l])
    # print('cameras are '+ str(c2l) + ' timesteps ahead of lidar')
    # print(dcs.index[int(drs.loc['camera']['lidar'])])
# dpe['date'] = pd.date_range('00:00', periods=len(dpe), freq='10ms')
# di['date'] = pd.date_range('00:00', periods=len(dpe), freq='10ms')

# dc['date'] = pd.date_range('00:00', periods=len(dc), freq='125ms')
# dl['date'] = pd.date_range('00:00', periods=len(dl), freq='100ms')

# %% Plot timeseries with subplots

press_gauges = ['06','11']

for inst in press_gauges:
    # fig, axs = plt.subplots(3, figsize=(20,12), facecolor='white')
    # axs[0].plot(dpe['date'], dpe['p' + inst], c='k', lw=0.1, label='press resamp')
    # axs[0].scatter(dpe['date'], dpe['p' + inst], s = 0.3, marker='.', c='k', label='press resamp')
    # axs[0].plot(di['date'], di['p' + inst], c='grey', lw=0.1, label='press resamp')
    # axs[0].scatter(di['date'], di['p' + inst], s = 0.3, marker='.', c='grey', label='press resamp')
    # axs[1].plot(dc['date'], dc['c' + inst], c='k', lw=0.1, label='camera')
    # axs[1].scatter(dc['date'], dc['c' + inst], marker='.',s = 3,  c='k', label='camera')
    # axs[2].plot(dl['date'], dl['l' + inst], c='k', lw=0.2, label='press resamp')
    # axs[2].scatter(dl['date'], dl['l' + inst], marker='.',s = 3,  c='k', label='press resamp')
    # axs[0].set_title('(a) Pressure Gage',loc='left')
    # axs[1].set_title('(b) Stereo Reconstruction',loc='left')
    # axs[2].set_title('(c) Lidar',loc='left')
    # axs[0].get_xaxis().set_ticklabels([])
    # axs[1].get_xaxis().set_ticklabels([])
    # axs[2].set_xlabel('time')
    # for ax in fig.get_axes():
    #     ax.set_xlim((dpe['date'].iloc[0],dpe['date'].iloc[4000]))
    #     ax.set_ylim(-0.35,0.35)
    #     ax.grid(True, alpha = 0.1)
    #     ax.set_ylabel('$\eta$ (m)')
    # plt.savefig(figfolder + 'timeseries_p' + inst + '.png', bbox_inches='tight')
    # plt.show()
    
    
    for tp in range(0,9):
        trange = [tp, tp+1]
        if inst == '11':
            fig, axs = plt.subplots(1, figsize=(80,3), facecolor='white') # was 20x4
        elif inst == '06':
            fig, axs = plt.subplots(1, figsize=(80,4), facecolor='white') # was 20x6
        axs.plot(dpe['date'], dpe['p' + inst], c='y', lw=0.4, label='hyd press')
        axs.scatter(dpe['date'], dpe['p' + inst], s = 0.3, marker='.', c='y')
        axs.plot(di['date'], di['p' + inst], c='k', lw=0.4, label='dac press')
        axs.scatter(di['date'], di['p' + inst], s = 0.3, marker='.', c='k')
        axs.plot(dc['date'], dc['c' + inst], c='r', lw=0.7, label='camera')
        axs.scatter(dc['date'], dc['c' + inst], marker='.',s = 4,  c='r')
        # axs[2].plot(dl['date'], dl['l12'], c?='k', lw=0.3, label='lidar')
        axs.plot(dl['date'], dl['l' + inst], c='b', lw=0.7, label='lidar')
        axs.scatter(dl['date'], dl['l' + inst], marker='.',s = 4,  c='b')
        axs.hlines(0, dl['date'].min(), dl['date'].max(), colors='grey', linestyles='-',alpha = 0.25)
        # axs[0].set_ylim((-.3,.3))
        axs.set_title('$H_s$ = ' + str(Tinfo['Hs']) + ' m, $T_p$ = ' + str(Tinfo['Tp']) +r' s, $\langle \eta \rangle $ = '
                      + str(Tinfo['watlev']) + r' m, $\sigma_{\theta}$ = ' + str(Tinfo['spread']) + '$^{\circ}$, Instrument: ' 
                      + 'p' + inst + ', $x$ = ' + str(round(xyz['p' + inst].iloc[0],1)) + ' m',loc='left')
        # axs.get_xaxis().set_ticklabels([])
        axs.set_xlabel('time (hh:mm:ss)')
        axs.set_xlim((dpe['date'].iloc[6000*trange[0]],dpe['date'].iloc[6000*trange[1]])) # was to 4000
        if inst == '11':
            axs.set_ylim(-0.2,0.2)
        elif inst == '06':
            axs.set_ylim(-0.2,0.55)
        axs.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
        axs.grid(True, alpha = 0.2)
        axs.set_ylabel('$\eta$ (m)')
        plt.savefig(figfolder + Tinfo['comp'] + '_timeseries_p' + inst + str(trange[0]) + '.png', bbox_inches='tight')
        plt.show()
        
# plt.close('all')

# press_gauges = ['06','11']

# for inst in press_gauges:
#     fig, axs = plt.subplots(3, figsize=(20,12), facecolor='white')
#     axs[0].plot(dpe['date'][:-p2l], dpe['p' + inst][p2l:], c='k', lw=0.1, label='press resamp')
#     axs[0].scatter(dpe['date'][:-p2l], dpe['p' + inst][p2l:], s = 0.3, marker='.', c='k', label='press resamp')
#     axs[0].plot(di['date'][:-p2l], di['p' + inst][p2l:], c='grey', lw=0.1, label='press resamp')
#     axs[0].scatter(di['date'][:-p2l], di['p' + inst][p2l:], s = 0.3, marker='.', c='grey', label='press resamp')
#     axs[1].plot(dc['date'][:-c2l], dc['c' + inst][c2l:], c='k', lw=0.1, label='camera')
#     axs[1].scatter(dc['date'][:-c2l], dc['c' + inst][c2l:], marker='.',s = 3,  c='k', label='camera')
#     axs[2].plot(dl['date'], dl['l' + inst], c='k', lw=0.2, label='press resamp')
#     axs[2].scatter(dl['date'], dl['l' + inst], marker='.',s = 3,  c='k', label='press resamp')
#     axs[0].set_title('(a) Pressure Gage',loc='left')
#     axs[1].set_title('(b) Stereo Reconstruction',loc='left')
#     axs[2].set_title('(c) Lidar',loc='left')
#     axs[0].get_xaxis().set_ticklabels([])
#     axs[1].get_xaxis().set_ticklabels([])
#     axs[2].set_xlabel('time')
#     for ax in fig.get_axes():
#         ax.set_xlim((dpe['date'].iloc[0],dpe['date'].iloc[4000]))
#         ax.set_ylim(-0.35,0.35)
#         ax.grid(True, alpha = 0.1)
#         ax.set_ylabel('$\eta$ (m)')
#     plt.savefig(figfolder + 'timeseries_p' + inst + '.png', bbox_inches='tight')
#     plt.show()
    
#     if inst == '11':
#         fig, axs = plt.subplots(1, figsize=(20,4), facecolor='white')
#     elif inst == '06':
#         fig, axs = plt.subplots(1, figsize=(20,6), facecolor='white')
#     # axs[0].plot(dpe['date'], dpe['p12'], c='k', lw=0.3, label='press')
#     axs.plot(dpe['date'][:-p2l], dpe['p' + inst][p2l:], c='y', lw=0.4, label='hyd press')
#     axs.scatter(dpe['date'][:-p2l], dpe['p' + inst][p2l:], s = 0.3, marker='.', c='y')
#     axs.plot(di['date'][:-p2l], di['p' + inst][p2l:], c='k', lw=0.4, label='dac press')
#     axs.scatter(di['date'][:-p2l], di['p' + inst][p2l:], s = 0.3, marker='.', c='k')
#     axs.plot(dc['date'][:-c2l], dc['c' + inst][c2l:], c='r', lw=0.7, label='camera')
#     axs.scatter(dc['date'][:-c2l], dc['c' + inst][c2l:], marker='.',s = 4,  c='r')
#     # axs[2].plot(dl['date'], dl['l12'], c?='k', lw=0.3, label='lidar')
#     axs.plot(dl['date'], dl['l' + inst], c='b', lw=0.7, label='lidar')
#     axs.scatter(dl['date'], dl['l' + inst], marker='.',s = 4,  c='b')
#     axs.hlines(0, dl['date'][0], dl['date'][-1:], colors='grey', linestyles='-',alpha = 0.25)
#     # axs[0].set_ylim((-.3,.3))
#     axs.set_title('$H_s$ = ' + str(Tinfo['Hs']) + ' m, $T_p$ = ' + str(Tinfo['Tp']) +r' s, $\langle \eta \rangle $ = '
#                   + str(Tinfo['watlev']) + r' m, $\sigma_{\theta}$ = ' + str(Tinfo['spread']) + '$^{\circ}$, Instrument: ' 
#                   + 'p' + inst + ', $x$ = ' + str(round(xyz['p' + inst].iloc[0],1)) + ' m',loc='left')
#     # axs.get_xaxis().set_ticklabels([])
#     axs.set_xlabel('time (hh:mm:ss)')
#     axs.set_xlim((dpe['date'].iloc[0],dpe['date'].iloc[4000]))
#     if inst == '11':
#         axs.set_ylim(-0.2,0.2)
#     elif inst == '06':
#         axs.set_ylim(-0.2,0.55)
#     axs.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
#     axs.grid(True, alpha = 0.2)
#     axs.set_ylabel('$\eta$ (m)')
#     plt.savefig(figfolder + Tinfo['comp'] + '_timeseries_p' + inst + '.png', bbox_inches='tight')
#     plt.show()
    
    
# %%
# fig, axs = plt.subplots(2, figsize=(20,14), facecolor='white', wratios=(3, 2))
# # gs = plt.gridspec.GridSpec(2, 1, width_ratios=[2, 3]) 
# axs[0].plot(dpe['date'][:-p2l], dpe['p' + inst][p2l:], c='y', lw=0.4, label='hyd press')
# axs[0].scatter(dpe['date'][:-p2l], dpe['p' + inst][p2l:], s = 0.3, marker='.', c='y')
# axs[0].plot(di['date'][:-p2l], di['p' + inst][p2l:], c='k', lw=0.4, label='dac press')
# axs[0].scatter(di['date'][:-p2l], di['p' + inst][p2l:], s = 0.3, marker='.', c='k')
# axs[0].plot(dc['date'][:-c2l], dc['c' + inst][c2l:], c='r', lw=0.7, label='camera')
# axs[0].scatter(dc['date'][:-c2l], dc['c' + inst][c2l:], marker='.',s = 4,  c='r')
# # axs[2].plot(dl['date'], dl['l12'], c?='k', lw=0.3, label='lidar')
# axs[0].plot(dl['date'], dl['l' + inst], c='b', lw=0.7, label='lidar')
# axs[0].scatter(dl['date'], dl['l' + inst], marker='.',s = 4,  c='b')
# axs[0].set_title('$H_s$ = ' + str(Tinfo['Hs']) + ' m, $T_p$ = ' + str(Tinfo['Tp']) +r' s, $\langle \eta \rangle $ = '
#               + str(Tinfo['watlev']) + r' m, $\sigma_{\theta}$ = ' + str(Tinfo['spread']) + '$^{\circ}$, Instrument: ' 
#               + 'p' + inst + ', $x$ = ' + str(round(xyz['p' + inst].iloc[0],1)) + ' m',loc='left')


# axs.hlines(0, dl['date'][0], dl['date'][-1:], colors='grey', linestyles='-',alpha = 0.25)
# # axs[0].set_ylim((-.3,.3))
# # axs.get_xaxis().set_ticklabels([])


# axs.set_xlabel('time (hh:mm:ss)')
# axs.set_xlim((dpe['date'].iloc[0],dpe['date'].iloc[4000]))
# # axs.set_ylim(-0.2,0.2)
# axs.set_ylim(-0.2,0.55)
# axs.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
# axs.grid(True, alpha = 0.2)
# axs.set_ylabel('$\eta$ (m)')
# plt.savefig(figfolder + Tinfo['comp'] + '_timeseries_p' + inst + '.png', bbox_inches='tight')
# plt.show()






plt.close('all')