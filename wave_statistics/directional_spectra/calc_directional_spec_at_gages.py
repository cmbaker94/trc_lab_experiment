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
plt.close('all')
import matplotlib as mpl
import pandas as pd
from datetime import datetime
from matplotlib import rc
# font = {'family' : 'normal',
#         'weight' : 'normal',
#         'size'   : 18}
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
    
    Options_ = []
    [SMout,EPout] = dirspec(ID,SM,EP,Options_)
    SMout = calc_dir_sprd(SMout)
    return SMout, EPout

def calc_energy_weighted_mean(spec,freq,freqrange):
    idmin = min(range(len(freq)), key=lambda i: abs(freq[i]-freqrange[0]))
    idmax = min(range(len(freq)), key=lambda i: abs(freq[i]-freqrange[1]))
    specint = np.trapz(spec[idmin:idmax],x=freq[idmin:idmax])
    emean = specint/(freq[idmax]-freq[idmin])
    return emean


# %% Define file locations and path

path = 'E:/'
sys.path.append("C:\\Users\\cmbaker9\\Documents\\MATLAB\\MTOOLS\\pyDIWASP")
from dirspec import dirspec
# from dirspec import dirspec as dirspec
now = datetime.now()
date = now.strftime("%m-%d-%Y")
figpath = path + 'figures/wave_statistic/directional/' + date + '/'
locpath = path + 'data/processed/insitu/'
ffaddition = 'allinst'
if not os.path.exists(figpath):
    os.mkdir(figpath)

rho = 1000
g = 9.81

# %% Insert Trial information
spread = 40
Hs = 0.3
h = 1.07
Tp = 2

# "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
# "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],

if spread == 0:
    if Hs == 0.3:
        Tinfo = {
            "Hs": 0.3,
            "Tp": 2,
            "h": 1.07,
            "spread": 0,
            "day": '7',
            "trial": '14',
            "insitu": "09-01-2018-2213UTC",
            "clpath": "09-01-2018-2214UTC",
            "timesection": "2229-2238",
            "camera": "09-01-2018-2155UTC_Scene1",
            "frames": "07200-11999",
            "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
            "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            "parray":'sz'
            }
    elif Hs == 0.25:
        Tinfo = {
            "Hs": 0.25,
            "Tp": 2,
            "h": 1.07,
            "spread": 0,
            "day": '7',
            "trial": '12',
            "insitu": "09-01-2018-2001UTC",
            "clpath": "09-01-2018-2002UTC",
            "timesection": "2017-2026",
            "camera": "09-01-2018-1950UTC_Scene1",
            "frames": "07200-11999",
            "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
            "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            "parray":'sz'
            }
elif spread == 20:
    if Tp == 2:
         Tinfo = {
             "Hs": 0.3,
             "Tp": 2,
             "h": 1.07,
             "spread": 20,
             "day": '5',
             "trial": '07',
             "insitu": "08-30-2018-2222UTC",
             "clpath": "08-30-2018-2222UTC",
             "timesection": "2237-2246",
             "camera": "08-30-2018-2222UTC_Scene1",
             "frames": "07200-11999",
             "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
             "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
             "parray":'sz'
             }
    elif Tp == 3:
        Tinfo = {
            "Hs": 0.3,
            "Tp": 3,
            "h": 1.07,
            "spread": 20,
            "day": '7',
            "trial": '09',
            "insitu": "09-01-2018-1626UTC",
            "clpath": "09-01-2018-1627UTC",
            "timesection": "1642-1651",
            "camera": "09-01-2018-1530UTC_Scene1",
            "frames": "07200-11999",
            "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
            "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            "parray":'sz'
            }
elif spread == 30:
    Tinfo = {
    "Hs": 0.3,
    "Tp": 2,
    "h": 1.07,
    "spread": 30,
    "day": '4',
    "trial": '06',
    "insitu": "08-29-2018-2255UTC",
    "clpath": "08-29-2018-2255UTC",
    "timesection": "2310-2319",
    "camera": "08-29-2018-2236UTC_Scene1",
    "frames": "07200-11999",
    "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
    "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
    "parray":'sz'
    }
elif spread == 40:
    if Hs == 0.3:
        if h == 1.07:
            Tinfo = {
                "Hs": 0.3,
                "Tp": 2,
                "h": 1.07,
                "spread": 40,
                "day": '5',
                "trial": '06',
                "insitu": "08-30-2018-2129UTC",
                "clpath": "08-30-2018-2129UTC",
                "timesection": "2144-2153",
                "camera": "08-30-2018-2119UTC_Scene1",
                "frames": "07200-11999",
                "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
                "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                "parray":'sz'
                }
        elif h == 1.00:
            Tinfo = {
                "Hs": 0.3,
                "Tp": 2,
                "h": 1.00,
                "spread": 40,
                "day": '6',
                "trial": '07',
                "insitu": "08-31-2018-2232UTC",
                "clpath": "08-31-2018-2232UTC",
                "timesection": "2247-2256",
                "camera": "08-31-2018-2225UTC_Scene1",
                "frames": "07200-11999",
                "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
                "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                "parray":'sz'
                }
    elif Hs == 0.25:
        Tinfo = {
        "Hs": 0.25,
        "Tp": 2,
        "h": 1.07,
        "spread": 40,
        "day": '5',
        "trial": '02',
        "insitu": "08-30-2018-1655UTC",
        "clpath": "08-30-2018-1655UTC",
        "timesection": "1710-1719",
        "camera": "08-30-2018-1634UTC_Scene1",
        "frames": "07200-11999",
        "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
        "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        "parray":'sz'
        }
    elif Hs == 0.2:
        Tinfo = {
        "Hs": 0.2,
        "Tp": 2,
        "h": 1.07,
        "spread": 40,
        "day": '5',
        "trial": '01',
        "insitu": "08-30-2018-1534UTC",
        "clpath": "08-30-2018-1534UTC",
        "timesection": "1549-1558",
        "camera": "08-30-2018-1518UTC_Scene1",
        "frames": "07200-11999",
        "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
        "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        "parray":'sz'
        }


conditions =  'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_tide' + str(int(Tinfo['h']*100)) + '_spread' + str(Tinfo['spread'])
datapath = 'G:/PRJ-1873/inter/random' + Tinfo['day'] + '/Trial' + Tinfo['trial'] + '/'
savepath = 'E:/data/processed/conditions/' + conditions + '/' + Tinfo['clpath'] + '/time_' + Tinfo['timesection'] + '/'

instruments = ['wg','press']
for i in range(len(instruments)):
    inst = instruments[i]
    [df,XYZ] = read_data(Tinfo,inst,datapath) # read text files
    # df2 = df.iloc[:, 3:11]
    # XYZ2 = XYZ[:,2:10]
    # [SMout,EPout] = calc_SM(Tinfo,inst,df2,XYZ2)
    [SMout,EPout] = calc_SM(Tinfo,inst,df,XYZ)
    exec('SM' + inst + '_' + Tinfo['parray'] + '= SMout')
    savepath_S = savepath + 'spec_' + inst + '_' + ffaddition + '.csv'
    savepath_f = savepath + 'freq_' + inst + '_' + ffaddition + '.csv'
    savepath_d = savepath + 'dirs_' + inst + '_' + ffaddition + '.csv'
    np.savetxt(savepath_f, SMout['freqs'], delimiter=",")
    np.savetxt(savepath_d, SMout['dirs'], delimiter=",")
    np.savetxt(savepath_S, SMout['S'].real, delimiter=",")
    del SMout
    
#%% Read surfzone array in situ PUV

puvpath = 'E:/data/processed/insitu/' + Tinfo['insitu'] + '/'
dname = 'PUV_p11.csv'
PUV11 = pd.read_csv(puvpath+dname,header=None,names=['freq','SSE','dir2','sig2']) 
dname = 'PUV_p06.csv'
PUV06 = pd.read_csv(puvpath+dname,header=None,names=['freq','SSE','dir2','sig2']) 

PUV11['SSE']=PUV11['SSE'].replace(0, np.nan)
PUV06['SSE']=PUV06['SSE'].replace(0, np.nan)

# %%  surf zone camera timeseries
# campath = '/Users/cmbaker9/Documents/Research/Lab_Experiments/data/processed/DEM/TRM-' + Tinfo['camera'] + '/frames_' + Tinfo['frames'] + '/'
# stepx = [0, 1, -1, 0, 0]
# stepy = [0, 0, 0, 1, -1]

# for i in range(len(stepx)):
#     fpath = campath + 'pressure_gage_xy_stepx' + str(stepx[i]) + '_stepy' + str(stepy[i]) + '.csv'
#     dpxyz = pd.read_csv(fpath, header=None, names=['x', 'y'])
#     dpxyz['z'] = np.zeros((12,1))+Tinfo['h']
#     XYZ = np.transpose(dpxyz.values)
    
#     fpath = campath + 'pressure_gage_timeseries_stepx' + str(stepx[i]) + '_stepy' + str(stepy[i]) + '.csv'
#     df = pd.read_csv(fpath, header=None, names=['z1', 'z2', 'z3', 'z4', 'z5', 'z6', 'z7', 'z8', 'z9', 'z10', 'z11', 'z12'])
    
#     inst = 'cam'
#     [SMout,EPout] = calc_SM(Tinfo,inst,df,XYZ)
#     exec('SM' + inst + '_' + Tinfo['parray'] + str(i) + '= SMout')
    
# # calc average
# SMcam_sz = dict()
# SMcam_sz['freqs'] = SMcam_sz0['freqs']
# SMcam_sz['dirs'] = SMcam_sz0['dirs']
# SMcam_sz['Sf'] = np.nanmean(np.vstack((SMcam_sz0['Sf'], SMcam_sz1['Sf'], SMcam_sz2['Sf'], SMcam_sz3['Sf'], SMcam_sz4['Sf'])),axis=0)
# SMcam_sz['Sd'] = np.nanmean(np.vstack((SMcam_sz0['Sd'], SMcam_sz1['Sd'], SMcam_sz2['Sd'], SMcam_sz3['Sd'], SMcam_sz4['Sd'])),axis=0)
# SMcam_sz['th_2'] = np.nanmean(np.vstack((SMcam_sz0['th_2'], SMcam_sz1['th_2'], SMcam_sz2['th_2'], SMcam_sz3['th_2'], SMcam_sz4['th_2'])),axis=0)
# SMcam_sz['sig_2'] = np.nanmean(np.vstack((SMcam_sz0['sig_2'], SMcam_sz1['sig_2'], SMcam_sz2['sig_2'], SMcam_sz3['sig_2'], SMcam_sz4['sig_2'])),axis=0)

# campath = '/Users/cmbaker9/Documents/Research/Lab_Experiments/data/processed/cameras/TRM-' + Tinfo['camera'] + '/frames_' + Tinfo['frames'] + '/'
campath = 'E:/data/processed/conditions/' + conditions + '/' + Tinfo['clpath'] + '/time_' + Tinfo['timesection'] + '/'
fpath = campath + 'cam_array_xy_5cm.csv' #_5cm
dpxyz = pd.read_csv(fpath, header='infer')#, names=['x', 'y'])
# dpxyz.append(np.zeros((1,12))+Tinfo['h'],ignore_index=True)
XYZ = dpxyz.values
z = np.zeros((1,12))+Tinfo['h']
XYZ = np.vstack([XYZ, z])


fpath = campath + 'cam_array_timeseries_5cm.csv'
# df = pd.read_csv(fpath, header=None, names=['z1', 'z2', 'z3', 'z4', 'z5', 'z6', 'z7', 'z8', 'z9', 'z10', 'z11', 'z12'])
df = pd.read_csv(fpath, header='infer')#, names=['z1', 'z10', 'z11', 'z12', 'z2', 'z3','z4', 'z5', 'z6', 'z7', 'z8', 'z9'])
df = df.interpolate(method='linear', limit_direction='both')#, axis='columns')

inst = 'cam'
[SMout,EPout] = calc_SM(Tinfo,inst,df,XYZ)
exec('SM' + inst + '_' + Tinfo['parray'] + '= SMout')

# %% lidar

lidarpath = 'E:/data/processed/conditions/' + conditions + '/' + Tinfo['clpath'] + '/time_' + Tinfo['timesection'] + '/'
fpath = lidarpath + 'lidar_array_xy.csv'
dpxyz = pd.read_csv(fpath, header='infer')#, names=['x', 'y'])
# dpxyz.append(np.zeros((1,12))+Tinfo['h'],ignore_index=True)
XYZ = dpxyz.values
z = np.zeros((1,12))+Tinfo['h']
XYZ = np.vstack([XYZ, z])


fpath = campath + 'lidar_array_timeseries.csv'
# df = pd.read_csv(fpath, header=None, names=['z1', 'z2', 'z3', 'z4', 'z5', 'z6', 'z7', 'z8', 'z9', 'z10', 'z11', 'z12'])
df = pd.read_csv(fpath, header='infer')#, names=['z1', 'z10', 'z11', 'z12', 'z2', 'z3','z4', 'z5', 'z6', 'z7', 'z8', 'z9'])
df = df.interpolate(method='linear', limit_direction='both')#, axis='columns')
dfin = df

inst = 'lid'
[SMout,EPout] = calc_SM(Tinfo,inst,df,XYZ)
exec('SM' + inst + '_' + Tinfo['parray'] + '= SMout')

# %% inner shelf trial
if spread == 0:
    if Hs == 0.3:
        Tinfo = {
            "Hs": 0.3,
            "Tp": 2,
            "h": 1.07,
            "spread": 0,
            "day": '9',
            "trial": '01',
            "insitu": "09-06-2018-1535UTC",
            "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
            "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            "parray":'is'
        }
    elif Hs == 0.25:
        Tinfo = {
            "Hs": 0.25,
            "Tp": 2,
            "h": 1.07,
            "spread": 0,
            "day": '9',
            "trial": '05',
            "insitu": "09-06-2018-1952UTC",
            "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
            "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            "parray":'is'
        }
elif spread == 20:
    Tinfo = {
    "Hs": 0.3,
    "Tp": 2,
    "h": 1.07,
    "spread": 20,
    "day": '9',
    "trial": '02',
    "insitu": "09-06-2018-1655UTC",
    "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
    "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
    "parray":'is'
    }
elif spread == 30:
    Tinfo = {
    "Hs": 0.3,
    "Tp": 2,
    "h": 1.07,
    "spread": 30,
    "day": '9',
    "trial": '03',
    "insitu": "09-06-2018-1748UTC",
    "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
    "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
    "parray":'is'
    }
elif spread == 40:
    if Hs == 0.3:
        Tinfo = {
        "Hs": 0.3,
        "Tp": 2,
        "h": 1.07,
        "spread": 40,
        "day": '9',
        "trial": '04',
        "insitu": '09-06-2018-1841UTC',
        "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
        "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        "parray":'is'
        }
    elif Hs == 0.25:
        Tinfo = {
        "Hs": 0.25,
        "Tp": 2,
        "h": 1.07,
        "spread": 40,
        "day": '9',
        "trial": '07',
        "insitu": '09-06-2018-2150UTC',
        "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
        "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        "parray":'is'
        }
    elif Hs == 0.2:
        Tinfo = {
        "Hs": 0.2,
        "Tp": 2,
        "h": 1.07,
        "spread": 40,
        "day": '9',
        "trial": '08',
        "insitu": '09-06-2018-2248UTC',
        "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
        "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        "parray":'is'
        }

datapath = 'G:/PRJ-1873/inter/random' + Tinfo['day'] + '/Trial' + Tinfo['trial'] + '/'

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
# axs[0].plot(SMwg_is['freqs'],SMwg_is['Sf'], c='k', lw=1, linestyle='-.', label='wg is')
axs[0].plot(SMpress_sz['freqs'],SMpress_sz['Sf'], c='r', lw=1, label='press sz')
axs[0].plot(SMpress_is['freqs'],SMpress_is['Sf'], c='r', lw=1, linestyle='-.', label='press is')
axs[0].plot(SMcam_sz['freqs'],SMcam_sz['Sf'], c='b', lw=1, label='cam sz')
axs[0].plot(SMlid_sz['freqs'],SMlid_sz['Sf'], c='g', lw=1, label='lidar sz')
axs[0].plot(fsmall,fsmall4, c='m', lw=1, label='$f^{-4}$')
axs[0].set_xlim(SMwg_sz['freqs'].min(),2.5)
axs[0].set_yscale('log')
axs[0].set_ylim((10**-4.5),10**-1.5)
axs[0].set_xlabel(r'$f$ (Hz)')
axs[0].set_ylabel(r'$S_{f}$ (m$^{2}$/Hz)')
axs[0].legend(prop={'size': 10})
axs[0].grid(True, alpha = 0.2)

axs[1].plot(SMwg_sz['dirs'],SMwg_sz['Sd'], c='k', lw=1)
# axs[1].plot(SMwg_is['dirs'],SMwg_is['Sd'], c='k', lw=1, linestyle='-.')
axs[1].plot(SMpress_sz['dirs'],SMpress_sz['Sd'], c='r', lw=1)
axs[1].plot(SMpress_is['dirs'],SMpress_is['Sd'], c='r', lw=1, linestyle='-.')
axs[1].plot(SMcam_sz['dirs'],SMcam_sz['Sd'], c='b', lw=1)
axs[1].plot(SMlid_sz['dirs'],SMlid_sz['Sd'], c='g', lw=1)
axs[1].set_xlim(SMwg_sz['dirs'].min(),SMwg_sz['dirs'].max())
axs[1].set_yscale('log')
axs[1].set_ylim((10**-6.5),10**-2)
axs[1].set_xlabel('Deg. $(^{\circ})$')
axs[1].set_ylabel('$S_{d}$ (m$^{2}$/deg)')
axs[1].grid(True, alpha = 0.2)
plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sf_Sd_cam_' + ffaddition +'.png', bbox_inches='tight')
plt.show()

# %% Plots

fig, axs = plt.subplots(2,figsize=(8,7))
axs[0].plot(SMwg_sz['freqs'],SMwg_sz['th_2'], c='k', lw=1, label='wg sz')
# axs[0].plot(SMwg_is['freqs'],SMwg_is['th_2'], c='k', lw=1, linestyle='-.', label='wg is')
axs[0].plot(SMpress_sz['freqs'],SMpress_sz['th_2'], c='r', lw=1, label='press sz')
axs[0].plot(SMpress_is['freqs'],SMpress_is['th_2'], c='r', lw=1, linestyle='-.', label='press is')
axs[0].plot(SMcam_sz['freqs'],SMcam_sz['th_2'], c='b', lw=1, label='cam sz')
axs[0].plot(SMlid_sz['freqs'],SMlid_sz['th_2'], c='g', lw=1, label='lid sz')
axs[0].set_xlim(SMwg_sz['freqs'].min(),2.5)
axs[0].set_ylim(-100,100)
axs[0].set_xlabel('$f$ (Hz)')
axs[0].set_ylabel('$\Theta (^{\circ})$')
axs[0].legend(prop={'size': 10})
axs[0].grid(True, alpha = 0.2)

axs[1].plot(SMwg_sz['freqs'],SMwg_sz['sig_2'], c='k', lw=1) #, label='Data')
# axs[1].plot(SMwg_is['freqs'],SMwg_is['sig_2'], c='k', lw=1, linestyle='-.')
axs[1].plot(SMpress_sz['freqs'],SMpress_sz['sig_2'], c='r', lw=1)
axs[1].plot(SMpress_is['freqs'],SMpress_is['sig_2'], c='r', lw=1, linestyle='-.')
axs[1].plot(SMcam_sz['freqs'],SMcam_sz['sig_2'], c='b', lw=1)
axs[1].plot(SMlid_sz['freqs'],SMlid_sz['sig_2'], c='g', lw=1)
axs[1].set_xlim(SMwg_sz['freqs'].min(),2.5)
axs[1].set_ylim(0,45)
axs[1].set_xlabel('$f$ (Hz)')
axs[1].set_ylabel('$\sigma_{\Theta} (^{\circ}$)')
axs[1].grid(True, alpha = 0.2)
plt.rcParams["font.family"] = "serif"
plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_theta_sigma_cam_' + ffaddition +'.png', bbox_inches='tight')
plt.show()

# %% Plot 2d spectra plot
# cax = [10**-6, 10**-4]

# Specwg = SMwg_sz['S'].real
# Specwg[Specwg<cax[0]] = cax[0] 

# fig, axs = plt.subplots(1,figsize=(8,7))
# c = axs.pcolor(SMwg_sz['freqs'], SMwg_sz['dirs'], np.transpose(Specwg), cmap='plasma', vmin=cax[0], vmax=cax[1]) # norm=mpl.colors.LogNorm(),
# axs.set_title('Directional Spread ' + str(Tinfo['spread']) + '$^{\circ}$' + ', Wave Gauges')
# plt.rcParams["font.family"] = "serif"
# cbar = fig.colorbar(c, ax=axs)
# cbar.set_label('$S(f,\Theta)$', rotation=270)
# axs.set_xlim(0,2.5)
# axs.set_xlabel('$f$ (Hz)')
# axs.set_ylabel('$\Theta (^{\circ}$)')
# plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sfd_wg_' + ffaddition +'.png', bbox_inches='tight')

# Specpress = SMpress_sz['S'].real
# Specpress[Specpress<cax[0]] = cax[0] 

# fig, axs = plt.subplots(1,figsize=(8,7))
# c = axs.pcolor(SMpress_sz['freqs'], SMpress_sz['dirs'], np.transpose(Specpress), cmap='plasma', vmin=cax[0], vmax=cax[1]) # norm=mpl.colors.LogNorm(),
# axs.set_title('Directional Spread ' + str(Tinfo['spread']) + '$^{\circ}$' + ', Surfzone Pressure Gauges')
# plt.rcParams["font.family"] = "serif"
# cbar = fig.colorbar(c, ax=axs)
# cbar.set_label('$S(f,\Theta)$', rotation=270)
# axs.set_xlim(0,2.5)
# axs.set_xlabel('$f$ (Hz)')
# axs.set_ylabel('$\Theta (^{\circ}$)')
# plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sfd_press_' + ffaddition +'.png', bbox_inches='tight')

# Speccam = SMcam_sz['S'].real
# Speccam[Speccam<cax[0]] = cax[0] 

# fig, axs = plt.subplots(1,figsize=(8,7))
# c = axs.pcolor(SMpress_sz['freqs'], SMpress_sz['dirs'], np.transpose(Speccam), cmap='plasma', vmin=cax[0], vmax=cax[1]) # norm=mpl.colors.LogNorm(),
# axs.set_title('Directional Spread ' + str(Tinfo['spread']) + '$^{\circ}$' + ', Surfzone Pressure Gauges')
# plt.rcParams["font.family"] = "serif"
# cbar = fig.colorbar(c, ax=axs)
# cbar.set_label('$S(f,\Theta)$', rotation=270)
# axs.set_xlim(0,2.5)
# axs.set_xlabel('$f$ (Hz)')
# axs.set_ylabel('$\Theta (^{\circ}$)')
# plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sfd_cam_' + ffaddition +'.png', bbox_inches='tight')


# %%
from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm
cax = [1e-6, 1e-3]


values = np.transpose(SMwg_sz['S'].real)
values[values<cax[0]]=cax[0]
values[values>cax[1]]=cax[1]

theta = np.array(SMwg_sz['dirs'])
azimuths = np.radians(SMwg_sz['dirs'])
zeniths = np.array(SMwg_sz['freqs'][0:2001])
 
values = np.array(values[:, 0:2001])
values = values.reshape(len(azimuths), len(zeniths))
 
r, theta = np.meshgrid(zeniths, azimuths)
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
# CS = plt.contourf(theta, r, values, levels=[10**-5, 10**-4, 10**-3, 10**-2],  cmap='plasma')
cnorm= LogNorm(vmin=cax[0], vmax=cax[1])
# CS = plt.contourf(theta, r, values, locator=ticker.LogLocator(subs=range(1,10)),  cmap='autumn', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
CS = plt.contourf(theta, r, values, levels=np.logspace(-6,-3,31),  cmap='viridis', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
cbaxes = fig.add_axes([0.86, 0.1, 0.03, 0.7])
cbar = fig.colorbar(CS, ticks=[1e-6, 1e-5, 1e-4, 1e-3], cax = cbaxes)
cbar.ax.set_yticklabels(['$10^{-6}$','$10^{-5}$', '$10^{-4}$', '$10^{-3}$'])
cbar.set_label(r'$S(f,\Theta)$',labelpad=-40, y=1.1, rotation=0)


rlabels = ax.get_ymajorticklabels()
for label in rlabels[:-1]:
    label.set_color('white')

plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sfd_polar_wg_10-6_' + ffaddition +'.png', bbox_inches='tight')

#%%

from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm
cax = [1e-5, 1e-3]


values = np.transpose(SMwg_sz['S'].real)
values[values<cax[0]]=cax[0]
values[values>cax[1]]=cax[1]

theta = np.array(SMwg_sz['dirs'])
azimuths = np.radians(SMwg_sz['dirs'])
zeniths = np.array(SMwg_sz['freqs'][0:2001])
 
values = np.array(values[:, 0:2001])
values = values.reshape(len(azimuths), len(zeniths))
 
r, theta = np.meshgrid(zeniths, azimuths)
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
# CS = plt.contourf(theta, r, values, levels=[10**-5, 10**-4, 10**-3, 10**-2],  cmap='plasma')
cnorm= LogNorm(vmin=cax[0], vmax=cax[1])
# CS = plt.contourf(theta, r, values, locator=ticker.LogLocator(subs=range(1,10)),  cmap='autumn', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
CS = plt.contourf(theta, r, values, levels=np.logspace(-5,-3,31),  cmap='viridis', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
cbaxes = fig.add_axes([0.86, 0.1, 0.03, 0.7])
cbar = fig.colorbar(CS, ticks=[1e-5, 1e-4, 1e-3], cax = cbaxes)
cbar.ax.set_yticklabels(['$10^{-5}$', '$10^{-4}$', '$10^{-3}$'])
cbar.set_label(r'$S(f,\Theta)$',labelpad=-40, y=1.1, rotation=0)


rlabels = ax.get_ymajorticklabels()
for label in rlabels[:-1]:
    label.set_color('white')

plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sfd_polar_wg_' + ffaddition +'.png', bbox_inches='tight')


# %%ax = [1e-5, 1e-3]
from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm
cax = [1e-6, 1e-4]

values = np.transpose(SMpress_sz['S'].real)
values[values<cax[0]]=cax[0]
values[values>cax[1]]=cax[1]

theta = np.array(SMpress_sz['dirs'])
azimuths = np.radians(SMpress_sz['dirs'])
zeniths = np.array(SMpress_sz['freqs'][0:2001])
 
values = np.array(values[:, 0:2001])
values = values.reshape(len(azimuths), len(zeniths))
 
r, theta = np.meshgrid(zeniths, azimuths)
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
# CS = plt.contourf(theta, r, values, levels=[10**-5, 10**-4, 10**-3, 10**-2],  cmap='plasma')
cnorm= LogNorm(vmin=cax[0], vmax=cax[1])
# CS = plt.contourf(theta, r, values, locator=ticker.LogLocator(subs=range(1,10)),  cmap='autumn', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
CS = plt.contourf(theta, r, values, levels=np.logspace(-6,-4,31),  cmap='viridis', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
cbaxes = fig.add_axes([0.86, 0.1, 0.03, 0.7])
cbar = fig.colorbar(CS, ticks=[1e-6, 1e-5, 1e-4], cax = cbaxes)
cbar.ax.set_yticklabels(['$10^{-6}$', '$10^{-5}$', '$10^{-4}$'])
cbar.set_label(r'$S(f,\Theta)$',labelpad=-40, y=1.1, rotation=0)


rlabels = ax.get_ymajorticklabels()
for label in rlabels[:-1]:
    label.set_color('white')

plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sfd_polar_press_sz_' + ffaddition +'.png', bbox_inches='tight')

## is
values = np.transpose(SMpress_is['S'].real)
values[values<cax[0]]=cax[0]
values[values>cax[1]]=cax[1]

theta = np.array(SMpress_is['dirs'])
azimuths = np.radians(SMpress_is['dirs'])
zeniths = np.array(SMpress_is['freqs'][0:2001])
 
values = np.array(values[:, 0:2001])
values = values.reshape(len(azimuths), len(zeniths))
 
r, theta = np.meshgrid(zeniths, azimuths)
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
# CS = plt.contourf(theta, r, values, levels=[10**-5, 10**-4, 10**-3, 10**-2],  cmap='plasma')
cnorm= LogNorm(vmin=cax[0], vmax=cax[1])
# CS = plt.contourf(theta, r, values, locator=ticker.LogLocator(subs=range(1,10)),  cmap='autumn', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
CS = plt.contourf(theta, r, values, levels=np.logspace(-6,-4,31),  cmap='viridis', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
cbaxes = fig.add_axes([0.86, 0.1, 0.03, 0.7])
cbar = fig.colorbar(CS, ticks=[1e-6, 1e-5, 1e-4], cax = cbaxes)
cbar.ax.set_yticklabels(['$10^{-6}$', '$10^{-5}$', '$10^{-4}$'])
cbar.set_label(r'$S(f,\Theta)$',labelpad=-40, y=1.1, rotation=0)


rlabels = ax.get_ymajorticklabels()
for label in rlabels[:-1]:
    label.set_color('white')

plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sfd_polar_press_is_' + ffaddition +'.png', bbox_inches='tight')


# %% camera

from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm
cax = [1e-6, 1e-4]

values = np.transpose(SMcam_sz['S'].real)
values[values<cax[0]]=cax[0]
values[values>cax[1]]=cax[1]

theta = np.array(SMcam_sz['dirs'])
azimuths = np.radians(SMcam_sz['dirs'])
zeniths = np.array(SMcam_sz['freqs'][0:2001])
 
values = np.array(values[:, 0:2001])
values = values.reshape(len(azimuths), len(zeniths))
 
r, theta = np.meshgrid(zeniths, azimuths)
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
# CS = plt.contourf(theta, r, values, levels=[10**-5, 10**-4, 10**-3, 10**-2],  cmap='plasma')
cnorm= LogNorm(vmin=cax[0], vmax=cax[1])
# CS = plt.contourf(theta, r, values, locator=ticker.LogLocator(subs=range(1,10)),  cmap='autumn', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
CS = plt.contourf(theta, r, values, levels=np.logspace(-6,-4,31),  cmap='viridis', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
cbaxes = fig.add_axes([0.86, 0.1, 0.03, 0.7])
cbar = fig.colorbar(CS, ticks=[1e-6, 1e-5, 1e-4], cax = cbaxes)
cbar.ax.set_yticklabels(['$10^{-6}$', '$10^{-5}$', '$10^{-4}$'])
cbar.set_label(r'$S(f,\Theta)$',labelpad=-40, y=1.1, rotation=0)


rlabels = ax.get_ymajorticklabels()
for label in rlabels[:-1]:
    label.set_color('white')

plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sfd_polar_camera_sz_' + ffaddition +'.png', bbox_inches='tight')

# plt.close('all')

# %% lidar

from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm
cax = [1e-6, 1e-4]

values = np.transpose(SMlid_sz['S'].real)
values[values<cax[0]]=cax[0]
values[values>cax[1]]=cax[1]

theta = np.array(SMlid_sz['dirs'])
azimuths = np.radians(SMlid_sz['dirs'])
zeniths = np.array(SMlid_sz['freqs'][0:2001])
 
values = np.array(values[:, 0:2001])
values = values.reshape(len(azimuths), len(zeniths))
 
r, theta = np.meshgrid(zeniths, azimuths)
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
# CS = plt.contourf(theta, r, values, levels=[10**-5, 10**-4, 10**-3, 10**-2],  cmap='plasma')
cnorm= LogNorm(vmin=cax[0], vmax=cax[1])
# CS = plt.contourf(theta, r, values, locator=ticker.LogLocator(subs=range(1,10)),  cmap='autumn', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
CS = plt.contourf(theta, r, values, levels=np.logspace(-6,-4,31),  cmap='viridis', norm=cnorm)#, vmin=0.0001, vmax=0.01, cmap='plasma')
cbaxes = fig.add_axes([0.86, 0.1, 0.03, 0.7])
cbar = fig.colorbar(CS, ticks=[1e-6, 1e-5, 1e-4], cax = cbaxes)
cbar.ax.set_yticklabels(['$10^{-6}$', '$10^{-5}$', '$10^{-4}$'])
cbar.set_label(r'$S(f,\Theta)$',labelpad=-40, y=1.1, rotation=0)


rlabels = ax.get_ymajorticklabels()
for label in rlabels[:-1]:
    label.set_color('white')

plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_Sfd_polar_lidar_sz_' + ffaddition +'.png', bbox_inches='tight')

# plt.close('all')

#%% plot puv camera pressure gage


fig, axs = plt.subplots(3,figsize=(8,10))
# axs[0].plot(SMwg_sz['freqs'],SMwg_sz['Sf'], c='k', lw=1, label='wg sz')
# axs[0].plot(SMwg_is['freqs'],SMwg_is['Sf'], c='k', lw=1, linestyle='-.', label='wg is')
axs[0].plot(PUV06['freq'],PUV06['SSE']/(10**4), c='k', lw=1.5, label='PUV p06')
axs[0].plot(PUV11['freq'],PUV11['SSE']/(10**4), c='k',linestyle='-.', lw=1.5, label='PUV p11')
axs[0].plot(SMpress_sz['freqs'],SMpress_sz['Sf'], c='r', lw=1.5, label='Press Array')
# axs[0].plot(SMpress_is['freqs'],SMpress_is['Sf'], c='r', lw=1, linestyle='-.', label='press is')
axs[0].plot(SMcam_sz['freqs'],SMcam_sz['Sf'], c='b', lw=1.5, label='Stereo Array')
axs[0].plot(SMlid_sz['freqs'],SMlid_sz['Sf'], c='g', lw=1.5, label='Lidar Array')
# axs[0].plot(fsmall,fsmall4, c='g', lw=1, label='$f^{-4}$')
axs[0].set_xlim(SMwg_sz['freqs'].min(),1.5)
axs[0].set_yscale('log')
axs[0].set_ylim((10**-4.5),10**-1.5)
axs[0].set_ylabel(r'$S_{f}$ (m$^{2}$/Hz)')
axs[0].xaxis.set_ticklabels([])
axs[0].legend(prop={'size': 13},loc=(0.48,1.04),ncol=2)#bbox_to_anchor=(0,1.02,1,0.2),loc='upper right')
# axs[0].legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", borderaxespad=0, ncol=4)
axs[0].grid(True, alpha = 0.2)

# axs[1].plot(SMwg_sz['freqs'],SMwg_sz['th_2'], c='k', lw=1, label='wg sz')
# axs[0].plot(SMwg_is['freqs'],SMwg_is['th_2'], c='k', lw=1, linestyle='-.', label='wg is')
axs[1].plot(PUV06['freq'],PUV06['dir2'], c='k', lw=1.5, label='PUV')
axs[1].plot(PUV11['freq'],PUV11['dir2'], c='k', lw=1.5,linestyle='-.', label='PUV')
axs[1].plot(SMpress_sz['freqs'],SMpress_sz['th_2'], c='r', lw=1.5, label='press sz')
# axs[1].plot(SMpress_is['freqs'],SMpress_is['th_2'], c='r', lw=1, linestyle='-.', label='press is')
axs[1].plot(SMcam_sz['freqs'],SMcam_sz['th_2'], c='b', lw=1.5, label='cam sz')
axs[1].plot(SMlid_sz['freqs'],SMlid_sz['th_2'], c='g', lw=1.5, label='lid sz')
axs[1].set_xlim(SMwg_sz['freqs'].min(),1.5)
axs[1].set_ylim(-100,100)
axs[1].set_ylabel('$\Theta~(^{\circ})$')
axs[1].xaxis.set_ticklabels([])
# axs[1].legend(prop={'size': 10})
axs[1].grid(True, alpha = 0.2)

# axs[2].plot(SMwg_sz['freqs'],SMwg_sz['sig_2'], c='k', lw=1) #, label='Data')
# axs[1].plot(SMwg_is['freqs'],SMwg_is['sig_2'], c='k', lw=1, linestyle='-.'
axs[2].plot(PUV06['freq'],PUV06['sig2']*180/np.pi, c='k', lw=1.5, label='p06')
axs[2].plot(PUV11['freq'],PUV11['sig2']*180/np.pi, c='k', lw=1.5, linestyle='-.', label='p11')
axs[2].plot(SMpress_sz['freqs'],SMpress_sz['sig_2'], c='r', lw=1.5)
# axs[2].plot(SMpress_is['freqs'],SMpress_is['sig_2'], c='r', lw=1, linestyle='-.')
axs[2].plot(SMcam_sz['freqs'],SMcam_sz['sig_2'], c='b', lw=1.5)
axs[2].plot(SMlid_sz['freqs'],SMlid_sz['sig_2'], c='g', lw=1.5)
axs[2].set_xlim(SMwg_sz['freqs'].min(),1.5)
axs[2].set_ylim(0,45)
axs[2].set_xlabel('$f$ (Hz)')
axs[2].set_ylabel('$\sigma_{\Theta}~(^{\circ}$)')
axs[2].grid(True, alpha = 0.2)
plt.rcParams["font.family"] = "serif"
plt.savefig(figpath + 'Hs' + str(int(Tinfo['Hs']*100)) + '_Tp' + str(Tinfo['Tp']) + '_sprd' + str(Tinfo['spread']) + '_h' + str(int(Tinfo['h']*100)) + '_SEE_dir_sig_sz_' + ffaddition +'.png', bbox_inches='tight')
plt.show()


#%% Create scatter plot
freqrange = [0.1, 1.0]

dirmean = dict()
sprdmean = dict()

dirmean['PUV06'] = calc_energy_weighted_mean(PUV06['dir2'], PUV06['freq'], freqrange)
dirmean['PUV11'] = calc_energy_weighted_mean(PUV11['dir2'], PUV11['freq'], freqrange)
dirmean['wg'] = calc_energy_weighted_mean(SMwg_sz['th_2'], SMwg_sz['freqs'], freqrange)
dirmean['press_sz'] = calc_energy_weighted_mean(SMpress_sz['th_2'], SMpress_sz['freqs'], freqrange)
dirmean['press_is'] = calc_energy_weighted_mean(SMpress_is['th_2'], SMpress_is['freqs'], freqrange)
dirmean['cam_sz'] = calc_energy_weighted_mean(SMcam_sz['th_2'], SMcam_sz['freqs'], freqrange)
dirmean['lid_sz'] = calc_energy_weighted_mean(SMlid_sz['th_2'], SMlid_sz['freqs'], freqrange)


sprdmean['PUV06']  = calc_energy_weighted_mean(PUV06['sig2']*180/np.pi, PUV06['freq'], freqrange)
sprdmean['PUV11']  = calc_energy_weighted_mean(PUV11['sig2']*180/np.pi, PUV06['freq'], freqrange)
sprdmean['wg'] = calc_energy_weighted_mean(SMwg_sz['sig_2'], SMwg_sz['freqs'], freqrange)
sprdmean['press_sz'] = calc_energy_weighted_mean(SMpress_sz['sig_2'], SMpress_sz['freqs'], freqrange)
sprdmean['press_is'] = calc_energy_weighted_mean(SMpress_is['sig_2'], SMpress_is['freqs'], freqrange)
sprdmean['cam_sz'] = calc_energy_weighted_mean(SMcam_sz['sig_2'], SMcam_sz['freqs'], freqrange)
sprdmean['lid_sz'] = calc_energy_weighted_mean(SMlid_sz['sig_2'], SMlid_sz['freqs'], freqrange)

ddir = pd.DataFrame()
ddir = ddir.append(dirmean, ignore_index=True)

dsprd = pd.DataFrame()
dsprd = dsprd.append(sprdmean, ignore_index=True)

savepath = campath + 'mean_direction.csv'
ddir.to_csv(savepath, index = False, header=True)

savepath = campath + 'directional_spread.csv'
dsprd.to_csv(savepath, index = False, header=True)
