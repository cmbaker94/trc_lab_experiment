# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 12:51:59 2021

@author: cmbaker9
calc energy weighted mean
"""

import numpy as np

def calc_energy_weighted_mean(spec,freq,freqrange):
    idmin = min(range(len(freq)), key=lambda i: abs(freq[i]-freqrange[0]))
    idmax = min(range(len(freq)), key=lambda i: abs(freq[i]-freqrange[1]))
    specint = np.trapz(spec[idmin:idmax],x=freq[idmin:idmax])
    emean = specint/(freq[idmax]-freq[idmin])
    return emean