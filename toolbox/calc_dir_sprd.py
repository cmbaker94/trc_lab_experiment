# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 12:50:21 2021

@author: cmbaker9

Function to compute directional spectra
"""
import numpy as np

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