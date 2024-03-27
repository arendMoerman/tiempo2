"""!
@file
This file contains the possible source signals for a TiEMPO2 simulation.

Currently, can only use MockSZ.
"""

from dataclasses import dataclass, field
import shutil
import os
from tqdm import tqdm
import numpy as np
import scipy.special as scp
import scipy.ndimage as sci
from scipy.integrate import dblquad
from multiprocessing import Pool
from functools import partial
import matplotlib.pyplot as pt

def generateSZMaps(SZsourceDict, clog, convolve_beam=True, telescopeDict=None, save=False, sourceName=None, path=None, trace_src=None):
    """!
    Generate SZ maps on-sky over a range of frequencies.

    @param SZsourceDict An SZsource dictionary.
    @param convolve_beam Convolve generated with diffraction-limited beams at specified frequencies. If True, need to specifiy a telescopeDict in order to get Dtel.
    @param telescopeDict A telescope dictionary. Need not be given if convolve_beam is False.
    @param trace_src Wether or not to evaluate SZ signal on a b grid (None) or on a trace of coordinates (numpy array of size 2x3).

    @returns SZ Datacube containing SZ intensity sky maps across the frequency range.
    @returns Az Array of Azimuthal co-ordinates.
    @returns El Array of Elevation co-ordinates.
    """

    try:
        import MockSZ.Models as MModel
        import MockSZ.Constants as MConst
    except:
        clog.error("Could not import MockSZ. Check if it is installed.")
        exit(1)


    Te      = SZsourceDict.get("Te")
    ne0     = SZsourceDict.get("ne0")               # n / cm^-3 -> n / m^-3
    thetac  = SZsourceDict.get("thetac")      # arcsec
    beta    = SZsourceDict.get("beta")
    Da      = SZsourceDict.get("Da")      # Mpc -> pc -> m
    v_pec   = SZsourceDict.get("v_pec")               # km / s -> m / s

    lims_Az = SZsourceDict.get("Az")                 # arcsec
    lims_El = SZsourceDict.get("El")                 # arcsec

    nAz = SZsourceDict.get("nAz")
    nEl = SZsourceDict.get("nEl")

    freq_Hz = SZsourceDict.get("freqs_src") * 1e9               # GHz -> Hz 

    Az, El = np.mgrid[lims_Az[0]:lims_Az[1]:nAz * 1j, lims_El[0]:lims_El[1]:nEl * 1j]
    theta = np.sqrt(Az**2 + El**2)

    simObjIso = MModel.IsoBetaModel(param=Te, v_pec=v_pec, no_CMB=True)
    
    if trace_src is not None:
        isob = simObjIso.getIsoBeta(trace_src[0,:], trace_src[1,:], beta, ne0, thetac, Da, grid=False)
    else:
        isob = simObjIso.getIsoBeta(Az[:,0], El[0,:], beta, ne0, thetac, Da, grid=True)
    
    SZ = simObjIso.getIsoBetaCube(isob, freq_Hz)

    #if convolve_beam:
    #    clog.info("Convolving beam patterns with SZ maps...")
    #    SZ = _convolveMaps(SZ, Az, El, freq_Hz, telescopeDict.get("Dtel"))
    
    SZ *= (MConst.c / freq_Hz)**2 
    
    # Now, transpose Az and El s.t. az is semif-fast, el is slow
    #SZ = np.transpose(SZ, axes=(1,0,2))

    print(SZ.shape)
    if save:
        _saveSource(sourceName, path, SZ, Az[:,0], El[0,:], freq_Hz)
    return SZ, Az[:,0], El[0,:]

def loadSZMaps(SZsourceDict):
    sourceName = SZsourceDict.get("filename")
    path = SZsourceDict.get("path")

    totalPath = os.path.join(path, sourceName)
    
    SZ = np.load(os.path.join(totalPath, "SZ.npy"))
    Az = np.load(os.path.join(totalPath, "Az.npy"))
    El = np.load(os.path.join(totalPath, "El.npy"))
    freqs = np.load(os.path.join(totalPath, "freqs.npy"))

    return SZ, Az, El, freqs

def generateGalSpecMaps(GalSpecSourceDict, telescopeDict, convolve_beam=True, save=False, sourceName=None, path=None):
    """!
    Generate GalSpec maps on-sky over a range of frequencies.

    @param GalSpecSourceDict A GalSpec source dictionary.
    @param convolve_beam Convolve generated with diffraction-limited beams at specified frequencies. If True, need to specifiy a telescopeDict in order to get Dtel.
    @param telescopeDict A telescope dictionary. Need not be given if convolve_beam is False.

    @returns SZ Datacube containing SZ intensity sky maps across the frequency range.
    @returns CMB Cosmic Microwave Background. Can be added to SZ map.
    @returns Az Array of Azimuthal co-ordinates.
    @returns El Array of Elevation co-ordinates.
    """
    
    import galspec

    lims_Az = GalSpecSourceDict.get("Az")                 # degree
    lims_El = GalSpecSourceDict.get("El")                 # degree

    nAz = GalSpecSourceDict.get("nAz")
    nEl = GalSpecSourceDict.get("nEl")

    Az, El = np.mgrid[lims_Az[0]:lims_Az[1]:nAz * 1j, lims_El[0]:lims_El[1]:nEl * 1j]
    theta = np.sqrt(Az**2 + El**2)

    gal_freq, gal_flux = galspec.spectrum(GalSpecSourceDict.get("lum"),
                                          GalSpecSourceDict.get("z"),
                                          GalSpecSourceDict.get("f_lo"),
                                          GalSpecSourceDict.get("f_hi"),
                                          GalSpecSourceDict.get("nfreqs"),
                                          GalSpecSourceDict.get("lwidth"))
    gal_cube = np.zeros((nAz, nEl, gal_freq.size))
    
    Az0 = np.argwhere(Az[:,0] == 0)
    El0 = np.argwhere(El[0,:] == 0)

    gal_cube[Az0, El0, :] = gal_flux

    gal_cube *= (telescopeDict.get("Dtel")/2)**2 * np.pi / (3e8)**2 * (gal_freq*1e9)**2 * 1e-26

    pt.plot(gal_cube[Az0, El0, :].ravel())
    pt.show()

    if convolve_beam:
        gal_cube = _convolveMaps(gal_cube, Az, El, gal_freq*1e9, telescopeDict.get("Dtel"))# + CMB

    #if save:
    #    _saveSource(sourceName, path, SZ, CMB, Az[:,0], El[0,:], freq_Hz)


    return gal_cube, Az[:,0], El[0,:], gal_freq

def _saveSource(sourceName, path, CMB, Az, El, freqs):
    totalPath = os.path.join(path, sourceName)
    
    if os.path.exists(totalPath):
        shutil.rmtree(totalPath)
    os.makedirs(totalPath)
    
    np.save(os.path.join(totalPath, "SZ"), SZ)
    np.save(os.path.join(totalPath, "Az"), Az)
    np.save(os.path.join(totalPath, "El"), El)
    np.save(os.path.join(totalPath, "freqs"), freqs)

def _convolveMaps(SZ, Az, El, freqs, Dtel, numThreads=None):
    numThreads = os.cpu_count() if numThreads is None else numThreads
    k = 2 * np.pi * freqs / (3e8)
    chunks_k = np.array_split(k, numThreads)
    chunks_SZ = np.array_split(SZ, numThreads, axis=-1)

    args = zip(chunks_k, chunks_SZ, np.arange(0, numThreads))

    _parallelFuncPartial = partial(_parallelConvolve, 
                                   Rtel=Dtel/2,
                                   Az=Az,
                                   El=El)
    
    pool = Pool(numThreads)
    res = np.concatenate(pool.map(_parallelFuncPartial, args), axis=-1)
    return res

def _parallelConvolve(args, Rtel, Az, El):
    k, SZ, thread_id = args
    
    iterator = lambda x, idx : tqdm(enumerate(x), ncols=100, total=x.size, colour="GREEN") if idx == 0 else enumerate(x) 
    for i, ki in iterator(k, thread_id):
        lim = 3 * np.pi * 1.2 / ki / Rtel
        
        weights = _AiryDisk(Az, El, ki, Rtel)
        norm = np.sum(weights)
        weights /= norm
        SZ[:,:,i] = sci.convolve(SZ[:,:,i], weights)
    return SZ

def _AiryDisk(Az, El, k, R):
    np.seterr(divide='ignore', invalid='ignore')
    theta = np.radians(np.sqrt((Az/3600)**2 + (El/3600)**2))
    airy = np.nan_to_num((2*scp.j1(k * R * np.sin(theta)) / (k * R * np.sin(theta)))**2, nan=1)
    return airy 

if __name__ == "__main__":

    generateGalSpecMaps(GSDict)
