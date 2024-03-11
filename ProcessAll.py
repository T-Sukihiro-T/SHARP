# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 13:52:37 2022

@author: Jarrett Ames & Christian Glover

(P.S.) We thank you all for the kind emails and pictures given to us by the members of 
Deng-Group, we are truley grateful for all of the support and hospitality when preforming our duties 
for the Data Analysis and our research. 

We hope for another great year in the office and to see you all in the comming year as we apply to the University Of Michigan.

And most of all, thank you to Professor Deng and Mr. Lydick for making all of this possible, and also making this summer one of the best that we've had in a long time.

 -Christian and Jarrett 
 
(P.P.S, Samantha will live on in our hearts forever)

"""

from ProcessSpectrum import ProcessSpectrum

from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from parse_filename import parse_filename
from adjustText import adjust_text
from getfilter import getfilter
from Polaritons import lower_polariton, plot_polariton, lower_polariton_with_fixed_arg, upper_polariton
from Polaritons import N_Polaritons
from scipy.optimize import curve_fit
from fitting import multi_curve_fit
from showimg import showimg

plot3texts = []

plot4texts = []

powertexts = []

points = []

MI = []

directory = '2021-03-17 g1'

'''
The for-loop below takes the given directory and determines if the filename is 
is equal to f, in which it then calls the process spectrum function to prepare
the individual spectrum graphs for plotting
'''

if 0:
    for i,f in enumerate(sorted(glob(directory+'/*.spe'))[::]):
        # if i<= 1:
        if 1 or f == directory:
            "Process spectrum uses the filename and its gif to turn the files into images and place them onto individual figures"
            x,y,w,background,model,parameters = ProcessSpectrum(f, gif = True, wavelengthmin=805,wavelengthmax=825)
            #break
        else:
            print(repr(f))
# raise Exception()

'''
The for-loop below takes the parameters from the Polariton functions
for each filename at different graph positions and plots them with a line of best fit.
'''

for idx,f in enumerate(sorted(glob(directory+'/*.npz'))[1:]):
    with np.load(f) as data:
        params = data['fit_parameters']
        exposuretime = data['exposuretime']
        wavelength = data['wavelength']
        img = data['img'][0,0]
    info = parse_filename(f)    
    spefilename = f[:-len('Processed.npz')]
    
    #Here, the graphs with multiple polaritons are plotted by using fixed parameters that make the lines of best fit plot more acurately.
        
    background = params[:,0]
    plt.figure(1)
    plt.plot(background)
    plt.figure(2)
    intensity = params[:,np.arange(1,params.shape[1],3,dtype=int)]
    width = np.abs(params[:,np.arange(2,params.shape[1],3,dtype=int)])
    position = params[:,np.arange(3,params.shape[1],3,dtype=int)]
    
    plt.plot(position,'o')
    plt.gca().invert_yaxis()
    
    #Here below are the fixed values that correct the lines of best fit seen in the multi-polariton graphs
    
    k0 = 128
    Ec0 = 1.517
    Ac = 0.023
    Ex = 1.525
    G = -2.380
    Ex2 = 1.531
    G2 = -2.250
    
    p0 = (128, 1.517, 0.023/128**2)
        
    #Variables for corrected intensity, max indicies of intensity, the position of the graph at max intensity, the max intensity, and a variable that can find nan values in lists are defined below.
    
    correctedI = np.pad(intensity,((0,0),(0,1)))
    maxIdcs = np.nanargmax(correctedI,axis=1)
    
    position_at_max = np.pad(position,((0,0),(0,1)))[[*range(len(maxIdcs))],maxIdcs]
    
    maxintensity = correctedI[[*range(len(maxIdcs))],maxIdcs]
    
    isnotnan = ~(np.isnan(maxintensity) | (maxintensity < 0)) & (position_at_max > 0)
    
    #The code below takes all of the maxintensity and position values, and clears the nan values that may be included in their indexes after being bounded between 650 and 65000
   
    xpts = np.arange(len(maxintensity))[30:230]
    allx = np.stack([xpts]*position.shape[1],axis=-1)
    use = (intensity[30:230] >= 650) &  (intensity[30:230] <= 65000)
    x = allx[use]
    y = 1240/position[30:230][use]
    # nn=np.all(np.isnan(position),axis=1)
    notnan = ~np.isnan(y)
    x = x[notnan]
    y = y[notnan]    
    "N_Polaritons is a function that takes its parameters from a filename and uses them to plot a polaritons line of best fit"
    #Lambda function used for replacing the parameters in the N_Polaritons function with fixed values
    fixed_coupling_strength = lambda k, k0, Ec0, Ac, plot = False: N_Polaritons(k, k0, Ec0, Ac, Ex, 10**-2.380, Ex2, 10**-2.250, plot = plot)
    
    # N_Polaritons(x,*p0,plot=True)
    
    if len(x) < len(p0):
        print(f"Too few points remaining in {spefilename}. Not fitting.")
        continue
    try:     
    
        "Multi-Curve fit utilizes fixed values for the N_Polaritons function, and creates a line of best fit for the plotted Polaritons"
        parameters, _ = multi_curve_fit(fixed_coupling_strength, x, y, p0=p0, which = np.digitize(y, bins = [Ex, Ex2]))#, bounds = ([120,1, 1, 0, 0.001], [140,2, 2, np.inf, .2] ))
        #lower_parameters, _ = multi_curve_fit(lower_polariton, intensity, position, p0, which = position <= 815, bounds = ([120,1, 1, 0, 0.001], [140,2, 2, np.inf, .2] ))
        #upper_parameters, _ = multi_curve_fit(upper_polariton, intensity, position, p0, which = position <= 815, bounds = ([120,1, 1, 0, 0.001], [140,2, 2, np.inf, .2] ))
        # lower_polariton.add(min = 0.001)
        # if lower_polariton(Ec0) > 5:
        #     ex_energy = partial(lower_polariton, 1, 0, 1.5, 1, 0.3)
        # showimg(np.arange(len(img)),1240/wavelength,np.log10(img.T)[::-1,:],cmap='Blues',nonlinear=False)
        
    except Exception as e: 
        print(e)
        # raise e
        continue
    
    #Below, the 20 energy graphs are plotted first in the program, being labeled with their rows and energies. 
    #These values are also bound between the x values of 30 and 230
    
    plt.figure()
    showimg(np.arange(len(img)),1240/wavelength,np.log10(img.T),cmap='Blues',nonlinear=True)
    plt.gca().set_prop_cycle(None)
    fixed_coupling_strength(xpts,*parameters, plot = True)
    plt.scatter(allx,1240/position[30:230],c=np.log10(intensity[30:230]), s=1, marker='x',cmap='summer')
    plt.scatter(x,y,c=np.log10(intensity[30:230][use]), s=2, cmap='spring')
    plt.title("Graph Position: " + info.mode_polariton_wl)
    plt.xlabel("Row")
    plt.ylabel("Energy (eV)")
    
    plt.ylim(1.51,1.541)
    
    
plt.show()
'''
The for-loop below takes the file name from the directory and plots the peaks in the intensity graphs
'''
for idx,f in enumerate(sorted(glob(directory+'/*.npz'))[::]):
    with np.load(f) as data:
        params = data['fit_parameters']
        exposuretime = data['exposuretime']
    info = parse_filename(f)    
    spefilename = f[:-len('Processed.npz')]
    
    #Here the figures that plot the peak graphs are created, along with the variables of intensity, width, and position being defined for the x and y axis.
        
    background = params[:,0]
    plt.figure(1)
    plt.plot(background)
    plt.figure(2)
    intensity = params[:,np.arange(1,params.shape[1],3,dtype=int)]
    width = np.abs(params[:,np.arange(2,params.shape[1],3,dtype=int)])
    position = params[:,np.arange(3,params.shape[1],3,dtype=int)]
    
    #Here the intensity is crected and its maximums are defined inside of maxintensity
    
    correctedI = np.pad(intensity,((0,0),(0,1)))
    maxIdcs = np.nanargmax(correctedI,axis=1)
    #intensity = imgs.intensity
    maxintensity = correctedI[[*range(len(maxIdcs))],maxIdcs]
    x = maxintensity
    peaks, properties = find_peaks(x, prominence = (10000,90000), width =(0,5))
    plt.plot(peaks, x[peaks], "x", color="k")
    
    plt.xlabel("Row")
    plt.ylabel("Intensity")
    
    #This loop takes the left and right properties of the ips values and fills empty ones with NaN's 
    for i in range(len(peaks)):
        leftips = properties["left_ips"][i]
        rightips = properties["right_ips"][i]
        left = int(np.ceil(leftips))
        right = int(np.ceil(rightips))
        maxintensity[left:right]=np.nan
        
    #paramlist[:,np.arange(3,paramlist.shape[1],3,dtype=int)],'o',ms=2, color = 'white', linewidth=1
    
    #The code below plots the polaritons with the maxintensity in each row on the graph, and uses different markers to show the difference in each type of pump mode (also utilizes a log scale to present more efficient results).
    if any(maxintensity>75e3):
        plt.figure()
        plt.plot(maxintensity, label = f)
        plt.xlabel("Row")
        plt.ylabel("Energy (eV)")
        print(f)
        plt.legend()
    # plt.plot(intensity)
    
    polarization = info.polarization
    
    MI.append(maxintensity)
    
    maxintensity[maxintensity>75e3]=np.nan
    
    
    plt.plot(maxintensity)
    print(info.power)
    x = info.power
    y = np.nanmax(maxintensity)*getfilter(name = info.filter)/exposuretime
    width_at_max = np.pad(width,((0,0),(0,1)))[[*range(len(maxIdcs))],maxIdcs]
    position_at_max = np.pad(position,((0,0),(0,1)))[[*range(len(maxIdcs))],maxIdcs]
    z = width_at_max[np.nanargmax(maxintensity)]
    print(x,y,z)
    
    marker = 'o'
    if info.polarization == 'TE':
        marker = '+'
    elif info.polarization == 'TM':
        marker == 'x'
    else:
        marker == "o"
    
    plt.figure(3)
    plt.ylabel("MaxIntensity")
    plt.xlabel("Power (mV)")
    plt.plot(x,y,marker,ms=4,color='k',linewidth=1)
    maxintensity_label = plt.annotate(f"{info.pump_mode} {info.shape}", [x,y], fontsize = 8)
    plot3texts.append(maxintensity_label)
    
    plt.xscale("log")
    plt.yscale("log")
    
    plt.figure(4)
    
    plt.plot(x,z, marker,label = 'Line Width')   
    plt.ylabel("Width at MaxIntensity")
    plt.xlabel("Power(mV)")
    width_label = plt.annotate(f"{info.pump_mode} {info.shape}", [x,z], fontsize = 8)
    plot4texts.append(width_label)
    
    points.append((info,x,y,z))
    
    plt.xscale("log")
    plt.yscale("log")
    
    #The code below filters the maxintensity list of nans values, and only plots maxintensity values that have a power that is less than or equal to 60
    
    isnotnan = ~(np.isnan(maxintensity) | (maxintensity < 0)) & (position_at_max > 0)
    isnotnan[:40] = False
    isnotnan[-40:] = False
    
    x = np.arange(len(maxintensity))[isnotnan]
    y = 1240/position_at_max[isnotnan]
    w = np.sqrt(maxintensity)[isnotnan]
    
    if info.power <= 60 and isnotnan.any():
        model = np.poly1d(np.polyfit(x, y, 2, w=w))
        plt.figure()
        plt.scatter(x,y, c = maxintensity[isnotnan], s = 3)
        yl = plt.ylim()
        plt.xlabel("Row")
        plt.ylabel("Energy (eV)")
        for i,p0 in enumerate([(128, 1.524, model[0], model[2], 0.08),
                   (128, 1.524, model[0], model[2], 0.15),
                   (128, 1.524, 0.9*model[0], model[2], 0.08),
                   (128, 1.524, model[0], 0.9*model[2], 0.08)]):
            try:     
            
                parameters, _=curve_fit(lower_polariton, x, y, p0, 1/w, bounds = ([120,1, 1, 0, 0.001], [140,2, 2, np.inf, .2] ))
                # lower_polariton.add(min = 0.001)
                # if lower_polariton(Ec0) > 5:
                #     ex_energy = partial(lower_polariton, 1, 0, 1.5, 1, 0.3)
                
            except Exception as e: 
                print(e)
                continue
            
            
            
            residuals = y - lower_polariton(x, *parameters)    
            np.sum(abs(residuals)**2)
            if info.power <= 60:
                plot_polariton(np.arange(len(maxintensity)), *parameters, color=f'C{i}', label = np.sum(abs(residuals)**2))
                plt.ylim(*yl)
            plt.legend()
            #Just so you know, it was not the high-schoolers that named the next set of variables
            #It was a certain someone from the deng-group office. 
        
        #The code below sets fixed arguements for the polariton functions to plot the lines of best fit more acurately
        fixable_args = ('k0', 'Ex', 'Ec0', 'Ac', 'G')
        for idx, fixed_arg in enumerate(fixable_args):
            residuals_list = []
            fixedparams_list = []
            p0=[128, 1.524, model[0], model[2] if model[2] >= 0 else 1, 0.05]
            low = [120,1,1, 0, 0.0001]
            high = [140,2,2, np.inf, .07]
                
            L_spaceidx = np.linspace(low[idx], high[idx] if np.isfinite(high[idx]) else 2, 25)
            del p0[idx]
            del low[idx]
            del high[idx]
            fig,((p1,p2),(p3,p4))=plt.subplots(2,2,figsize=(10,10))
            for v in L_spaceidx:
                #samantha = lambda k, k0, Ec0, Ac, G : lower_polariton(k, k0, v, Ec0, Ac, G)
                # samantha = lower_polariton_with_fixed_arg(Ex=v)
                fix_argument = {fixable_args[idx] : v}
                samantha = lower_polariton_with_fixed_arg(**fix_argument)
                
                try:     
                    
                    parameters, _=curve_fit(samantha, x, y, p0, 1/w, bounds = 
                                            (low, high))
                    
                except Exception as e: 
                    print(e)
                    residuals_list.append(np.nan)
                    fixedparams_list.append(None)
                    continue
                
                residuals = y - samantha(x, *parameters)
                fixed_params = list(parameters)
                fixed_params.insert(idx, v)
                fixedparams_list.append(fixed_params)
                #The code below plots the polaritons with the fixed arguements and the included residuals.
                if v == np.min(L_spaceidx):
                    p3.scatter(x,y, c = maxintensity[isnotnan], s = 3)
                    fixed_params = list(parameters)
                    fixed_params.insert(idx, v)
                    plt.sca(p3)
                    plot_polariton(np.arange(len(maxintensity)), *fixed_params, color=f'C{i}', label = np.sum(abs(residuals)**2))
                    plt.title(f'k0={fixed_params[0]:3g}, Ex={fixed_params[1]:3g}, Ec0={fixed_params[2]:3g}, Ac={fixed_params[3]:3g}, G={fixed_params[4]:3g}', fontsize = 9)
                    plt.xlabel("Row")
                    plt.ylabel("Energy (eV)")
                    plt.ylim(*yl)
                    
                    
                if v == np.max(L_spaceidx):
                    p4.scatter(x,y, c = maxintensity[isnotnan], s = 3)
                    fixed_params = list(parameters)
                    fixed_params.insert(idx, v)
                    plt.sca(p4)
                    plot_polariton(np.arange(len(maxintensity)), *fixed_params, color=f'C{i}', label = np.sum(abs(residuals)**2))
                    plt.title(f'k0={fixed_params[0]:3g}, Ex={fixed_params[1]:3g}, Ec0={fixed_params[2]:3g}, Ac={fixed_params[3]:3g}, G={fixed_params[4]:3g}', fontsize = 9)
                    plt.xlabel("Row")
                    plt.ylabel("Energy (eV)")
                    plt.ylim(*yl)
                    
                residuals = y - samantha(x, *parameters)
                residuals_list.append(np.sum(abs(residuals)**2))
            
            plt.xlim(40, 220)
            
            p1.plot(L_spaceidx,residuals_list)
            p1.set_xlabel(fixable_args[idx])
            p1.set_ylabel('$\sum |residual|^2$')
            p1.set_yscale('log')
            
            if not np.isnan(residuals_list).all():
                plt.sca(p2)
                p2.scatter(x,y, c = maxintensity[isnotnan], s = 3)
                fixed_params = fixedparams_list[np.nanargmin(residuals_list)]
                plot_polariton(np.arange(len(maxintensity)), *fixed_params, color=f'C{i}')
                plt.title(f'k0={fixed_params[0]:3g}, Ex={fixed_params[1]:3g}, Ec0={fixed_params[2]:3g}, Ac={fixed_params[3]:3g}, G={fixed_params[4]:3g}', fontsize = 9)
                plt.xlabel("Row")
                plt.ylabel("Energy (eV)")
                plt.ylim(*yl)
        
        
#This if-statement adjusts the plotted texts in the pump mode graphs for clarity 
        
if 0:
    plt.figure(3)
    adjust_text(plot3texts, arrowprops=dict(arrowstyle="->", color="red", lw=1.5))    
    plt.figure(4)
    adjust_text(plot4texts, arrowprops=dict(arrowstyle="->", color="green", lw=1.5))

plt.figure(5)
kinds = set()    

'''
This for-loop takes the variables given by the presented 
filename, and adds them to the set in kinds
'''

for info ,x,y,z in points:
    kinds.add((info.polarization, info.pump_mode, info.shape, info.axis))

'''
The for-loops below plot the data inside of the set in kinds. 
This code also plots a multi-axis chart that plots the data for intensity, power
and linewidth; that has appropriatly colored lines to show the axis that the data
coresponds to. 
'''

for kind in kinds:
    plt.figure()
    plt.plot([1]) #so it will show the title
    plt.title(kind)
    marker = 'o'
    intensity_x = []
    intensity_y = []
    linewidth_x = []
    linewidth_z = []
    
    for info,x,y,z in points:
        if kind == (info.polarization, info.pump_mode, info.shape, info.axis):
            plt.plot(x,y,marker,ms=5,color='k',linewidth=1)
            plt.plot()
            plt.xscale("log")
            plt.yscale("log")
            intensity_x.append(x)
            intensity_y.append(y)
    plt.xlabel("Power")
    plt.ylabel("Intensity")
    
    plt.figure()
    for info,x,y,z in points:
        if kind == (info.polarization, info.pump_mode, info.shape, info.axis):
            plt.plot(x,z,marker,ms=5,color='k',linewidth=1)
            plt.xscale("log")
            linewidth_x.append(x)
            linewidth_z.append(z)
            plt.xlabel("Power")
            plt.ylabel("Linewidth")
    plt.figure() 
    
    linewidth_x, linewidth_z = zip(*sorted(zip(linewidth_x, linewidth_z), key=lambda pair: pair[0]))
    intensity_x, intensity_y = zip(*sorted(zip(intensity_x, intensity_y), key=lambda pair: pair[0]))
    plt.plot(intensity_x, intensity_y, '-*', color = 'red', lw=0.5)
    
    plt.xlabel("Power")
    plt.ylabel("Intensity", color = 'red')
    plt.gca().tick_params(axis = 'y', labelcolor = 'red')
    plt.yscale('log')
    plt.xscale('log')
    #color red
    
    plt.twinx().plot(linewidth_x, linewidth_z, '--o',lw=0.5)  
    
    # color C0
    
    plt.xlabel("Power")
    plt.ylabel("Linewidth", color = 'C0')
    plt.gca().tick_params(axis = 'y', labelcolor = 'C0')
    plt.title(kind)
    
A = [x]
B = [y]
C = [z]

data_evaluation = np.array([A,B,C])

#Prints the covarience matrix of the data set in kinds
covMatrix = np.cov(data_evaluation,bias=True)
print (covMatrix)

