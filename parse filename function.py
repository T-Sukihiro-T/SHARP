# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 13:31:02 2022

@author: Jarrett & Christian
"""

from collections import namedtuple
import re

Info = namedtuple('Info', 'mode polariton_wl pump_mode pump_wl shape spatial_filter axis polarization power filter')


def parse_filename(f):
    info = f.split()[0]
    info = info.split('_') 
 
    
    m = re.match(r'(\D*)([0-9.]+)nm' ,info[0])
    mode = m.group(1)
    polariton_wl = float(m.group(2))
    m2 = re.match(r'(\D*)([0-9.]+)nm' ,info[1])
    pump_mode = m2.group(1)
    pump_wl = float(m2.group(2))
    shape = {"f0": "2-3μm", "f1": "10-20μm", "f2": "ring pump"}  [info[2]]
    spatial_filter = {"s0": True, "s1": False} [info[3]]
    axis = info[4]
    polarization = info[5]
    m3 = re.match(r'([0-9.]+)mW', info[6])
    power = float(m3.group(1))
    if len(info) > 7 :
        filter = info[7]
    else:
        filter = ''
        
    
    
    parameters = [mode, polariton_wl, pump_mode, pump_wl, shape, spatial_filter, axis, polarization, power, filter]
    
    
    return Info (*parameters)
    





parse_filename('g817.3nm_p770nm_f2_s1_ky_TE_8.39mW 2021-03-19 21_21_25 8087.spe')