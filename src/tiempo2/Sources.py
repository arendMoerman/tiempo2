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
from multiprocessing import Pool
from functools import partial

def generateSZMaps(SZsourceDict, convolve_beam=False, telescopeDict=None, ret_unit="MJy", save=False, sourceName=None, path=None):
    """!
    Generate SZ maps on-sky over a range of frequencies.

    @param SZsourceDict An SZsource dictionary.
    @param convolve_beam Convolve generated with diffraction-limited beams at specified frequencies. If True, need to specifiy a telescopeDict in order to get Dtel.
    @param telescopeDict A telescope dictionary. Need not be given if convolve_beam is False.
    @param ret_unit Return maps in this unit. Choose between standard W m**-2 Hz**-1 sr**-1 ('SI'), MJy sr**-1 ('MJy') or mK ('mK').

    @returns SZ Datacube containing SZ intensity sky maps across the frequency range.
    @returns CMB Cosmic Microwave Background. Can be added to SZ map.
    @returns Az Array of Azimuthal co-ordinates.
    @returns El Array of Elevation co-ordinates.
    """
    
    import MockSZ.Models as MModel
    import MockSZ.Conversions as MConv
    import MockSZ.Backgrounds as MBack


    Te      = MConv.eV_Temp(SZsourceDict.get("Te") * 1e3)   # KeV -> eV -> K
    ne0     = SZsourceDict.get("ne0") * 1e2                 # n / cm^-2 -> n / m^-2
    rc      = MConv.pc_m(SZsourceDict.get("rc") * 1e3)      # Kpc -> pc -> m
    beta    = SZsourceDict.get("beta")
    Da      = MConv.pc_m(SZsourceDict.get("Da") * 1e6)      # Mpc -> pc -> m
    v_pec   = SZsourceDict.get("v_pec") * 1e3               # km / s -> m / s

    lims_Az = SZsourceDict.get("Az")                 # degree
    lims_El = SZsourceDict.get("El")                 # degree

    nAz = SZsourceDict.get("nAz")
    nEl = SZsourceDict.get("nEl")

    freq_Hz = SZsourceDict.get("freqs") * 1e9               # GHz -> Hz 

    isob = MModel.IsoBetaModel(Te, ne0, rc, beta, Da, v_pec)

    Az, El = np.mgrid[lims_Az[0]:lims_Az[1]:nAz * 1j, lims_El[0]:lims_El[1]:nEl * 1j]
    theta = np.sqrt(Az**2 + El**2)

    if ret_unit == "MJy":
        tSZ = MConv.SI_JySr(isob.tSZMap(theta, freq_Hz)) * 1e-6 # MJy / sr
        kSZ = MConv.SI_JySr(isob.kSZMap(theta, freq_Hz)) * 1e-6 # MJy / sr
        CMB = MConv.SI_JySr(MBack.getSpecificIntensityCMB(freq_Hz)) * 1e-6

    elif ret_unit == "SI":
        tSZ = isob.tSZMap(theta, freq_Hz)
        kSZ = isob.kSZMap(theta, freq_Hz)
        CMB = MBack.getSpecificIntensityCMB(freq_Hz)

    elif ret_unit == "mK":
        tSZ = MConv.SI_Temp(isob.tSZMap(theta, freq_Hz), freq_Hz) * 1e3 # mK
        kSZ = MConv.SI_Temp(isob.kSZMap(theta, freq_Hz), freq_Hz) * 1e3 # mK
        CMB = MConv.SI_Temp(MBack.getSpecificIntensityCMB(freq_Hz), freq_Hz) * 1e3

    SZ = tSZ + kSZ
    CMB = np.ones(SZ.shape) * CMB
    if convolve_beam:
        SZ = _convolveMaps(SZ, Az, El, freq_Hz, telescopeDict.get("Dtel"))# + CMB

    if save:
        _saveSource(sourceName, path, SZ, CMB, Az[:,0], El[0,:], freq_Hz)

    return SZ, CMB, Az[:,0], El[0,:]

def loadSZMaps(SZsourceDict):
    sourceName = SZsourceDict.get("filename")
    path = SZsourceDict.get("path")

    totalPath = os.path.join(path, sourceName)
    
    SZ = np.load(os.path.join(totalPath, "SZ.npy"))
    CMB = np.load(os.path.join(totalPath, "CMB.npy"))
    Az = np.load(os.path.join(totalPath, "Az.npy"))
    El = np.load(os.path.join(totalPath, "El.npy"))
    freqs = np.load(os.path.join(totalPath, "freqs.npy"))

    return SZ, CMB, Az, El, freqs

def _saveSource(sourceName, path, SZ, CMB, Az, El, freqs):
    totalPath = os.path.join(path, sourceName)
    
    if os.path.exists(totalPath):
        shutil.rmtree(totalPath)
    os.makedirs(totalPath)
    
    np.save(os.path.join(totalPath, "SZ"), SZ)
    np.save(os.path.join(totalPath, "CMB"), CMB)
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
    
    iterator = lambda x, idx : tqdm(enumerate(x), ncols=100, total=x.size) if idx == 0 else enumerate(x) 

    for i, ki in iterator(k, thread_id):
        weights = _AiryDisk(ki, Rtel, Az, El)
        SZ[:,:,i] = sci.convolve(SZ[:,:,i], weights)
    return SZ

def _AiryDisk(k, R, Az, El):
    np.seterr(divide='ignore', invalid='ignore')
    
    theta = np.radians(np.sqrt(Az**2 + El**2))
    airy = np.nan_to_num((2*scp.j1(k * R * np.sin(theta)) / (k * R * np.sin(theta)))**2, nan=1)
    airy /= np.sum(airy)
    return airy 

