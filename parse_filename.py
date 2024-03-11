# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 13:31:02 2022

@author: Glove
"""

from collections import namedtuple
import re
import os

Info = namedtuple('Info', 'mode_polariton_wl pump_mode pump_wl shape spatial_filter axis polarization power filter')


def parse_filename(f):
    
    f =  os.path.basename(f)
    
    info = f.split()[0]
    info = info.split('_') 
 
    mode_polariton_wl = info[0]
    m2 = re.match(r'(\D*)([0-9.]+)nm' ,info[1])
    pump_mode = m2.group(1)
    pump_wl = float(m2.group(2))
    shape = {"f0": "2-3μm", "f1": "10-20μm", "f2": "ring pump"}  [info[2]]
    spatial_filter = {"s0": True, "s1": False} [info[3]]
    axis = info[4]
    polarization = info[5] 
    def m3(group) :
        if info[6].endswith('uW'):
          return re.match(r'([0-9.]+)uW', info[6]),1e-3
        else:
          return re.match(r'([0-9.]+)mW', info[6]),1
    m3,factor = m3(1)     
    power = float(m3.group(1))*factor
    if len(info) > 7 :
        filter = info[7]
    else:
        filter = ''
        
    
    
    parameters = [mode_polariton_wl, pump_mode, pump_wl, shape, spatial_filter, axis, polarization, power, filter]
    
    
    return Info (*parameters)
    





print(parse_filename('2021-03-17 g1\\g817.3nm_p770nm_f2_s1_ky_TE_336uW 2021-03-19 21_23_31 8093.speProcessed.npz'))
