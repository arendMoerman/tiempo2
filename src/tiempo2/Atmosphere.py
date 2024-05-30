"""!
@file
This file handles all atmospheric screen generation. 

In particular, the EPL maps from ARIS are converted to Gaussian-smoothed PWV maps.
This script also contains functions that read atmospheric transmission curves, as function of frequency and PWV.
"""

import os
import numpy as np
import csv

from astropy.io import fits
from scipy.ndimage import gaussian_filter

def prepAtmospherePWV(atmosphereDict, telescopeDict, clog=None, number=1):
    """!   
    Prepare ARIS atmospheric screens for usage in TiEMPO2.
    This function works in the following way:
        1. Load an ARIS subchunk, stored under the 'path' key in atmosphereDict.
           This requires the 'filename' key to include the extension of the file.
        2. Remove ARIS metadata
        3. Convolve with 2D Gaussian
        4. Store ARIS data in subfolder (/prepd/) with same name.
    """
   
    filename = atmosphereDict.get("filename")
    path = atmosphereDict.get("path")
    pwv0 = atmosphereDict.get("PWV0")

    Rtel = telescopeDict.get("Dtel") / 2
    std = Rtel/np.sqrt( 2.*np.log(10.) )
    truncate = Rtel/std
   

    # Conversion from dEPL to dPWV from Smith-Weintraub relation.
    a = 6.3003663 

    prepd_path = os.path.join(path, "prepd")
    if not os.path.isdir(prepd_path):
        os.mkdir(prepd_path)

    #### AREND VOOR JOU MORGEN: Sla chunks op in groepen van 5, met 10 m/s heb je dan een screen van 1 uur

    test_l = []
    ny = 0
    nx = 0
    for subdir, dirs, files in os.walk(path):
        idx_l = np.array([int(x.split("-")[-1]) for x in files])
        idx_l_srt = np.argsort(idx_l)
        l_srt = np.array(files)[idx_l_srt]

        for i in range(number):
            subpath = os.path.join(path, l_srt[i])
            subchunk = np.loadtxt(subpath, delimiter=',')

            if i == 0: 
                ny = np.unique(subchunk[:,1]).size
            nx += np.unique(subchunk[:,0]).size

            dEPL = subchunk[:,2].reshape((np.unique(subchunk[:,0]).size, ny))
            PWV = pwv0 + (1./a * dEPL*1e-6)*1e+3 #in mm

            if i == 0:
                PWV_st = PWV
            else:
                PWV_st = np.concatenate((PWV_st, PWV), axis=0)
        PWV_Gauss = gaussian_filter(PWV_st, std, mode='mirror', truncate=truncate)

        return PWV_Gauss, nx, ny
        #for file in files:
        #    clog.info(f"Preparing {file}...")
        #
        #    file_split = file.split("-")

        #    file_idx = int(file_split[-1])
        #    
        #    subpath = os.path.join(path, file)
        #    subchunk = np.loadtxt(subpath, delimiter=',')

        #    n_subx = np.unique(subchunk[:,0]).size
        #    n_suby = np.unique(subchunk[:,1]).size

        #    x_file = os.path.join(prepd_path, f"{file_idx}_x.datp")
        #    np.savetxt(x_file, np.unique(subchunk[:,0]))

        #    if file_idx == 0: 
        #        y_file = os.path.join(prepd_path, f"y.datp")
        #        np.savetxt(y_file, np.unique(subchunk[:,1]))

        #    dEPL = subchunk[:,2].reshape((n_subx, n_suby))
        #    PWV = pwv0 + (1./a * dEPL*1e-6)*1e+3 #in mm
        #    
        #    PWV_Gauss = gaussian_filter(PWV, std, mode='mirror', truncate=truncate)

        #    PWV_path = os.path.join(prepd_path, f"{file_idx}_{file_split[0]}p")
        #    np.savetxt(PWV_path, PWV_Gauss)
            
            #extent = [0, 32, 0, 32]

            #import matplotlib.pyplot as pt
            #fig, ax = pt.subplots(1,2)
            #ax[0].imshow(PWV[:n_suby, :].T, extent=extent)
            #ax[0].set_xlabel(r"$x$ [m]")
            #ax[0].set_ylabel(r"$y$ [m]")
            #
            #ax[1].imshow(PWV_Gauss[:n_suby, :].T, extent=extent)
            #ax[1].set_xlabel(r"$x$ [m]")
            #ax[1].set_ylabel(r"$y$ [m]")
            #pt.show()
    
            #test_l.append(n_suby)
            #return test_l
    
    if clog is not None:
        clog.info(f"Finished preparing atmospheric screens.")
    else:
        print(f"Finished preparing atmospheric screens.")

