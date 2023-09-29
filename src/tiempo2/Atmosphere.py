"""!
@file
This file handles all atmospheric screen generation. 
In particular, the EPL maps from ARIS are converted to Gaussian-smoothed PWV maps.
"""

import glob
import os
import numpy as np
import csv

from astropy.io import fits
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from dataclasses import dataclass, field


#def generateAtmospherePWV(prefix_filename, path_data, pwv0, pwvgrid, max_windspeed, obs_duration, dish_radius,
 #           max_num_strips=1, x_length_strip=0, OFFdist=233.6, Rscan=60., EL0=60., h_column=1000., load_spline=False):
def generateAtmospherePWV(atmosphereDict, telescopeDict):
    """!   
    Generation of PWV maps from ARIS data

    @param atmosphereDict Dictionary containing atmosphere parameters.
    @param telescopeDict Dictionary containing telescope parameters
    """
    
    filename = atmosphereDict.get("filename")
    path = atmosphereDict.get("path")
    pwv0 = atmosphereDict.get("PWV0")

    Rtel = telescopeDict.get("Dtel") / 2

    flist = os.path.join(path, filename + '.fits')
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

    return PWV_Gauss, nx, ny

def readAtmTransmission():
    """!
    Parser for reading atmospheric transmission curves from csv file.

    @param line Line of text 
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

            
