# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 13:18:03 2021

@author: cmbaker9

read in situ data
"""

import numpy as np
import pandas as pd

def read_insitu_data(Tinfo,inst,datapath):
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