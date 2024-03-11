# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 14:25:31 2022

@author: Jarrett & Christian
"""

import os

import numpy as np



import spe_loader

from scipy.optimize import curve_fit

from scipy.signal import find_peaks

from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.use('agg')
plt.ion()

from Polaritons import lower_polariton
from parse_filename import parse_filename

def ProcessSpectrum (filename,gif = False,wavelengthmin=813,wavelengthmax=825):

    if gif:
        import imageio as io
    
    paramlist = []
    
    filenames = []
    
    
    imgs=spe_loader.load_from_files([filename])
    exposuretime_s = float(imgs.footer.SpeFormat.DataHistories.DataHistory.Origin.Experiment.Devices.Cameras.Camera.ShutterTiming.ExposureTime.cdata)/1000
    
    
    
    plt.imshow(imgs.data[0][0])
    plt.figure()
    #plt.plot(imgs.wavelength)
    
    
    def lorentzian (wavelengths, amplitude, width, position):
        return amplitude/(1+((wavelengths - position)/width)**2)
    
    def lorentzian_offset(wavelengths, offset, amplitude, width, position):
        return offset +  lorentzian (wavelengths, amplitude, width, position)
    
    def Nlorentzians (wavelengths, offset, *a):
        if (len(a)%3 != 0):
            raise Exception("Needs a multiple of three arguements")
        out = np.full_like(wavelengths, offset)
        for i in range (0, len(a), 3):
            out += lorentzian(wavelengths, *a[i: i + 3])
        return out
    
    
    for row in range(0,imgs.data[0][0].shape[0]-0):
    
    
        spectrum=(imgs.data[0][0][row,:])
    
        plt.figure()
        plt.plot(imgs.wavelength,spectrum)
        plt.title(f"Photoluminescence Row {row}")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Intensity (arb. units)")
        
        i = np.argmax(spectrum)
        wavelengths = imgs.wavelength
        maxwavelength = wavelengths[i]
        amplitude = 7000
        width = 0.3
        position = maxwavelength
        
        x = spectrum
        peaks, properties = find_peaks(x, prominence=25, height = (None, None), width =(2))
        plt.plot(wavelengths[peaks], x[peaks], "x", color="orange")
        #print(wavelengths[peaks])
        
       
        
        parameters = [np.average(spectrum)]
        bounds = [(0,parameters[0]*2)]
        if len(peaks) >= 1: 
            for i in range(len(peaks)):
                peakindex = peaks[i]
                peakheight = properties["peak_heights"][i]
                leftips = properties["left_ips"][i]
                rightips = properties["right_ips"][i]
                peakposition = wavelengths[peakindex]
                left_pos = wavelengths[round(leftips)]
                right_pos = wavelengths[round(rightips)]
                width = abs(right_pos-left_pos)
                parameters.extend([max(1,peakheight-parameters[0]), width, peakposition])
                bounds.extend([(0,80000), (0,wavelengthmax-wavelengthmin), (wavelengths[0],wavelengths[-1])])
        else:
             parameters.extend([max(1,np.max(spectrum)-parameters[0]), 2, maxwavelength])
             bounds.extend([(0,80000), (0,wavelengthmax-wavelengthmin), (wavelengths[0],wavelengths[-1])])
                            
        
        
        y = Nlorentzians(wavelengths, 600, amplitude, width, position)
        #plt.plot (wavelengths, y)
         
       
        
        #amp2 = 4000
        #width2 = 0.1
        #pos2 = maxwavelength + 0.5
        
        def Precision (wavelengths, offset, amp2, width2, pos2, amplitude, width, position):
            return offset + (amp2/(1+((wavelengths - pos2)/width2)**2)) + lorentzian(wavelengths, amplitude, width, position)
        
        try:
       #  parameters, _=curve_fit(Nlorentzians, wavelengths, spectrum, (1000, amplitude, width, position, amp2, width2, pos2))
            parameters, _=curve_fit(Nlorentzians, wavelengths, spectrum, parameters,bounds=tuple(zip(*bounds)))
        except Exception as e:
            print("Row Number:",row)
            print(e)
            print(parameters)
            print(bounds)
            paramlist.append((np.nan,))
            continue
        #print(parameters)
            
        y2 = Nlorentzians(wavelengths, *parameters)
        paramlist.append(parameters)
        
        plt.plot (wavelengths, y2)
        plt.xlim(wavelengthmin,wavelengthmax)
        
        print("Row Number:",row)
        
        
        
        imagename = f'{row}.png'
        filenames.append(imagename)
        if gif:
            plt.savefig(imagename)
        #plt.show()
        plt.close('all')
    if gif:
        with io.get_writer(f'{filename}{row}.gif', mode='I') as writer:
            for f in filenames:
                image = io.imread(f)
                writer.append_data(image)
                    
            for f in set(filenames):
                os.remove(f)   
      

    
    #%%
    
    arrayparams = np.full((len(paramlist), np.max(list(map(len,paramlist)))), np.nan)
    for i, params in enumerate(paramlist):
        arrayparams[i,0:len(params)] = params
        
    paramlist = arrayparams
    np.savez_compressed(f'{filename}Processed.npz',fit_parameters=paramlist,img=imgs.data,wavelength=imgs.wavelength,roi=imgs.roi, exposuretime = exposuretime_s)
    
    
    #%%
        
    plt.figure()
    plt.xlabel("Row")
    plt.ylabel("Intensity (arb. units)")
    
    with PdfPages(f'{filename}ProcessSpectrumOutput.pdf') as pdf:       
         
                
        print(paramlist)
        paramlist = np.array(paramlist)
        plt.plot(range(len(paramlist)),paramlist[:,np.arange(1,paramlist.shape[1],3,dtype=int)], label = 'Intensity')
        plt.plot(range(len(paramlist)),paramlist[:,[0]],'k', label = 'Background')
        background = np.nanmean(paramlist[:,[0]])
        plt.legend(loc='upper left')
        plt.twinx().plot(range(len(paramlist)),np.abs(paramlist[:,np.arange(2,paramlist.shape[1],3,dtype=int)]),':',label = 'Line Width')
        plt.legend(loc='upper right')
        
        
        
        plt.ylabel("Line Width (nm)")
        
        pdf.savefig()  # saves the current figure into a pdf page
          
        
        plt.figure()
        plt.xlabel("Row")
        plt.ylabel("Wavelength (nm)")
        
        
        
        
        plt.plot(range(len(paramlist)),paramlist[:,np.arange(3,paramlist.shape[1],3,dtype=int)],"o",ms=3)
        plt.ylim(wavelengthmax,wavelengthmin)
        pdf.savefig()  # saves the current figure into a pdf page
        print('index',i)
        
        plt.figure()
        plt.imshow(imgs.data[0][0].T)
        
        pdf.savefig()
        
        plt.figure()
        X,Y = np.meshgrid(np.arange(len(imgs.data[0][0])),wavelengths)
        plt.pcolormesh(X,Y, imgs.data[0][0].T,shading='nearest')
        plt.plot(range(len(paramlist)),paramlist[:,np.arange(3,paramlist.shape[1],3,dtype=int)],'o',ms=2, color = 'white', linewidth=1)
        plt.gca().invert_yaxis()
        plt.xlabel("Row")
        plt.ylabel("Wavelength (nm)")
        plt.ylim (wavelengthmax,wavelengthmin)
        
        pdf.savefig()
        
        plt.figure()
        sel=np.arange(3,paramlist.shape[1],3,dtype=int)
        #x=np.tile(np.arange(len(paramlist)),len(sel)).flatten()
        x = np.array([[i]*len(sel) for i in range(len(paramlist))]).flatten()
        y=1240/paramlist[:,sel].flatten()
        w=paramlist[:,sel-2].flatten() #Peak amplitude
        n=np.isnan(y)|np.isnan(w)
        n=~n
        n&=w>background
        x=x[n]
        y=y[n]
        w=w[n]
        
        
        np.polynomial.polynomial.Polynomial.fit(x, y, deg=2, w=w, window=None)
        #np.poly(np.array[np.arange(len(paramlist))])
        plt.show()
        model = np.poly1d(np.polyfit(x, y, 2, w=w))
        polyline = np.linspace(0, 255, 200)
        plt.scatter(x,y,c=w)
        plt.colorbar()
        plt.plot(polyline, model(polyline),'k')
        plt.xlabel("Row")
        plt.ylabel("Energy (eV)")
        plt.ylim (1240/wavelengthmax,1240/wavelengthmin)
        pdf.savefig()
        
        try:     
            
            n = w>background
            assert(n.any())
            parameters, _=curve_fit(lower_polariton, x[n], y[n], (128, 1.524, model[0], model[2], 0.03))#1/w[n])
        except Exception as e:
            print(e)
            print(background)
            print(np.max(w))
            parameters = (128, 1.524, model[0], model[2], 0.03)
        plt.figure()
        plt.plot(x,lower_polariton(x,*parameters), 'r:')
        plt.scatter(x[n],y[n],c=w[n])
        
        
        y2 = lower_polariton(x, *parameters)
        
        
        
        
        #plt.title(f'${model[2]:g}x^2 + {model[1]}x + {model[0]}$',fontsize=8)
        plt.xlabel("Row")
        plt.ylabel("Energy (eV)")
        plt.ylim (1240/wavelengthmax,1240/wavelengthmin)
        print(model)
        pdf.savefig()
        plt.show()
        
        
        return(x,y,w,background,model,parameters)

if __name__ == '__main__':
    from glob import glob
    
    for i,f in enumerate(sorted(glob('*.spe'))[::]):
        if i>= 0:
            x,y,w,background,model,parameters = ProcessSpectrum(f)
            #break
        
