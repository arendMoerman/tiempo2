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

def generateAtmospherePWV(atmosphereDict, telescopeDict, clog):
    """!   
    Generation of PWV maps from FITS file.
    NOTE: these FITS files are merged ARIS files.

    @param atmosphereDict Dictionary containing atmosphere parameters.
    @param telescopeDict Dictionary containing telescope parameters
    @param clog Custom logger.

    @returns A 2D array containing Gaussian-smoothed PWV maps of the atmosphere.
    @returns Number of x-coordinates.
    @returns Number of y-coordinates.
    """
    
    filename = atmosphereDict.get("filename")
    path = atmosphereDict.get("path")
    pwv0 = atmosphereDict.get("PWV0")

    Rtel = telescopeDict.get("Dtel") / 2

    flist = os.path.join(path, filename + '.fits')

    clog.info(f"Reading atmospheric screen from: {flist}")

    if len(flist):
        Loadfits = True
    else:
        Loadfits = False

    # Conversion from dEPL to dPWV from Smith-Weintraub relation.
    a = 6.3003663 

    num_strips = 1
    # Interpolate on PWV_Gauss
    
    for i in range(num_strips):
        if Loadfits:
            hdul = fits.open(os.path.join(path, '%s.fits' %filename), memmap=True)
            try:
                d = hdul[1].data
            except:
                d = hdul[0].data
        else:
            d = np.loadtxt(os.path.join(path, filename), delimiter=',')

        if i==0:
            nx = int( max(d[:, 0]) ) + 1
            ny = int( max(d[:, 1]) ) + 1
        epl= np.zeros( [nx, ny] )

        for j in range(len(d)):
            epl[ int(d[j, 0])-int(d[0, 0]), int(d[j, 1]) ] = int( d[j, 2] )

        if i == 0:
            dEPL_matrix = epl[:, :] # 31 corresponds to y-direction: 31*pwvgrid (m)
        else:
            dEPL_matrix = np.concatenate( (dEPL_matrix, epl[:, :]), axis=0)

        if Loadfits: hdul.close()
    PWV = pwv0 + (1./a * dEPL_matrix*1e-6)*1e+3 #in mm
    
    std = Rtel/np.sqrt( 2.*np.log(10.) )
    truncate = Rtel/std
   
    PWV_Gauss = gaussian_filter(PWV, std, mode='mirror', truncate=truncate)
    
    clog.info(f"Finished reading atmospheric screen.")
    import matplotlib.pyplot as pt
    pt.plot(PWV_Gauss[:,80])
    pt.show()

    return PWV_Gauss, nx, ny

def prepAtmospherePWV(atmosphereDict, telescopeDict, clog, number=0):
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

        print(l_srt)

        for i in range(number, number+5):
            subpath = os.path.join(path, l_srt[i])
            subchunk = np.loadtxt(subpath, delimiter=',')

            if i == number: 
                ny = np.unique(subchunk[:,1]).size
            print(ny)
            print(nx)
            nx += np.unique(subchunk[:,0]).size

            dEPL = subchunk[:,2].reshape((np.unique(subchunk[:,0]).size, ny))
            PWV = pwv0 + (1./a * dEPL*1e-6)*1e+3 #in mm

            if i == number:
                PWV_st = PWV
            else:
                PWV_st = np.concatenate((PWV_st, PWV), axis=0)
        PWV_Gauss = gaussian_filter(PWV_st, std, mode='mirror', truncate=truncate)
        #import matplotlib.pyplot as pt
        #pt.plot(PWV_Gauss[:90000,80])
        #pt.show()

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
    
    clog.info(f"Finished preparing atmospheric screens.")

def readAtmTransmissionText():
    """!
    Parser for reading atmospheric transmission curves from text file.
    
    @returns A 2D array containing atmospheric transmission.
    @returns Array containing frequencies at which atmospheric transmission is defined.
    @returns Array containing PWV values at which atmospheric transmission is defined.
    """

    dat_loc = os.path.join(os.path.dirname(__file__), "resources", "trans_data.dat")

    with open(dat_loc, 'r') as file:
        eta_atm = []
        freqs = []

        for line in file:
            if line[0] == "#":
                continue
            elif line[0] == " ":
                line = list(line.strip().split(" "))
                while("" in line):
                    line.remove("")
                pwv_curve = np.array(line, dtype=np.float64)
                continue
            
            line = list(line.strip().split(" "))
            while("" in line):
                line.remove("")
            
            eta_atm.append(np.array(line[1:]))
            freqs.append(line[0])

        freqs = np.array(freqs, dtype=np.float64)
        eta_atm = np.array(eta_atm, dtype=np.float64)

    return eta_atm.T, freqs, pwv_curve

def readAtmTransmissionCSV():
    """!
    Parser for reading atmospheric transmission curves from csv file.
    
    @returns A 2D array containing atmospheric transmission.
    @returns Array containing frequencies at which atmospheric transmission is defined.
    @returns Array containing PWV values at which atmospheric transmission is defined.
    """
    csv_loc = os.path.join(os.path.dirname(__file__), "resources", "atm.csv")

    with open(csv_loc, 'r') as file:
        csvreader = list(csv.reader(file, delimiter=" "))

        data = []

        for row in csvreader:
            if row[0] == "#":
                continue
            elif row[0] == "F":
                pwv_curve = np.array(row[1:], dtype=np.float64)
                continue
            while("" in row):
                row.remove("")
            data.append(np.array(row, dtype=np.float64))

        _arr = np.array(data)
        freqs = _arr[:,0]
        eta_atm = _arr[:,1:]
    return eta_atm, freqs, pwv_curve

            
