# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 09:32:52 2022

@author: Jarrett & Christian 
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA




def upper_polariton(k, k0, Ex, Ec0, Ac, G):
    
    k = k - k0
    a = 1
    b = -Ex-Ac*k**2-Ec0
    c = (Ac*k**2+Ec0)*Ex-G**2
    
    return  (-b+np.sqrt(b**2-4*a*c))/(2*a)

def lower_polariton(k, k0, Ex, Ec0, Ac, G):
    
    k = k - k0
    a = 1
    b = -Ex-Ac*k**2-Ec0
    c = (Ac*k**2+Ec0)*Ex-G**2
    
    return  (-b-np.sqrt(b**2-4*a*c))/(2*a)

def lower_polariton_with_fixed_arg(**kwargs):
    args = ['k', 'k0', 'Ex', 'Ec0', 'Ac', 'G']
    params = []
    for key in kwargs:
        if not key in args:
            raise KeyError(f'{key} not an allowable argument: {args}')
        params.append((key,kwargs[key]))
        args.remove(key)
    def fixed(*remaining):
        assert(len(remaining) == len(args))
        call_with = dict(params)
        for k,v in zip(args,remaining):
            call_with[k] = v
        return lower_polariton(**call_with)
        
    return fixed

def constants(k, Ex):
    
    return np.full_like(k, Ex, dtype=float)

def quadratic(k, k0, Ac,  Ec0):
    
    return (Ac*(k-k0)**2+Ec0)


def plot_polariton(k, k0, Ex, Ec0, Ac, G, color = 'black', label = None):
    plt.plot(k,lower_polariton(k, k0, Ex, Ec0, Ac, G),'--', color = color, label = label)
    plt.plot(k,upper_polariton(k, k0, Ex, Ec0, Ac, G),'--' ,color = color)
    plt.plot(k,constants(k, Ex), ':', color = color)
    plt.plot(k,quadratic(k, k0, Ac,  Ec0),':', color = color)
    
def N_Polaritons(k, k0, Ec0, Ac, Ex, G, *ExG, plot=False): 
    if len(ExG) % 2 == 1:
        raise TypeError('N_Polaritons requires a matching number of Exciton Energy values and Coupling strength values')
    C = Ec0 + Ac * (k - k0)**2
    Ex = [Ex] + [Ex for Ex in ExG[0::2]]
    G = [G] + [G for G in ExG[1::2]]
    # print(C,Ex,k)
    H = np.array([np.diag([Ck] + Ex) for Ck in C])
    # print(H)
    H[:,0,1:] = G
    H[:,1:,0] = G
    # print(H)
    Polariton_Energies = LA.eigvalsh(H)
    
    if plot:
        Cavity = C
        Exciton_lines = Ex
        plt.plot(k, Polariton_Energies,'k-')
        plt.plot(k, Cavity, 'k--')
        plt.hlines(Exciton_lines, xmin = k[0], xmax = k[-1], linestyle = '--', color = 'k')
        
    return Polariton_Energies

if __name__ == '__main__':
    N_Polaritons(np.linspace(-1,1,50),0,1,1,1,0.1,1.2,0.13,plot=1)
    
    fig, ax = plt.subplots()
    k = np.linspace(-2,2,1000)
    k0, Ex, Ec0, Ac, G =  0, 1.5, 1, 1, 0.3
    e = lower_polariton(k,k0, Ex, Ec0, Ac, G)
    d = upper_polariton(k,k0, Ex, Ec0, Ac, G)
    j = constants (k, Ex)
    g = quadratic(k, k0, Ac, Ec0)
    
    plt.plot(k,lower_polariton_with_fixed_arg(Ac=1)(k,k0,Ex,Ec0,G),'r-*')

   
    plt.plot(k,e)
    plt.plot(k,d)
    plt.plot(k,j,':')
    plt.plot(k,g,'--')
    plt.show()

