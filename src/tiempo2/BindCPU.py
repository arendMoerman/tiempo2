"""!
@file
Bindings for the ctypes interface for TiEMPO2. 
"""

import ctypes
import numpy as np
import os
import pathlib

import tiempo2.Threadmgr as TManager
import tiempo2.Structs as TStructs
import tiempo2.BindUtils as TBUtils

def loadTiEMPO2lib():
    """!
    Load the TiEMPO2 shared library. Will detect the operating system and link the library accordingly.

    @returns lib The ctypes library containing the C/C++ functions.
    """

    path_cur = pathlib.Path(__file__).parent.resolve()
    try:
        lib = ctypes.CDLL(os.path.join(path_cur, "libtiempo2.dll"))
    except:
        try:
            lib = ctypes.CDLL(os.path.join(path_cur, "libtiempo2.so"))
        except:
            lib = ctypes.CDLL(os.path.join(path_cur, "libtiempo2.dylib"))

    lib.runTiEMPO2.argtypes = [ctypes.POINTER(TStructs.Instrument), 
                               ctypes.POINTER(TStructs.Telescope),
                               ctypes.POINTER(TStructs.Atmosphere), 
                               ctypes.POINTER(TStructs.Source),
                               ctypes.POINTER(TStructs.Output),
                               ctypes.c_int, ctypes.c_int]
    
    lib.calcW2K.argtypes = [ctypes.POINTER(TStructs.Instrument), 
                            ctypes.POINTER(TStructs.Telescope),
                            ctypes.POINTER(TStructs.Atmosphere), 
                            ctypes.POINTER(TStructs.CalOutput),
                            ctypes.c_int, ctypes.c_int]
    
    lib.getSourceSignal.argtypes = [ctypes.POINTER(TStructs.Instrument), 
                                    ctypes.POINTER(TStructs.Telescope),
                                    ctypes.POINTER(ctypes.c_double), 
                                    ctypes.POINTER(ctypes.c_double),
                                    ctypes.c_double, ctypes.c_bool]

    lib.getChopperCalibration.argtypes = [ctypes.POINTER(TStructs.Instrument),
                                          ctypes.POINTER(ctypes.c_double),
                                          ctypes.c_double]
    
    lib.getEtaAtm.argtypes = [TStructs.ArrSpec, 
                              ctypes.POINTER(ctypes.c_double), 
                              ctypes.c_double]
    
    lib.getNEP.argtypes = [ctypes.POINTER(TStructs.Instrument), 
                           ctypes.POINTER(TStructs.Telescope),
                           ctypes.POINTER(ctypes.c_double), 
                           ctypes.c_double, ctypes.c_double]
    
    lib.runTiEMPO2.restype = None
    lib.calcW2K.restype = None
    lib.getSourceSignal.restype = None
    lib.getChopperCalibration.restype = None
    lib.getEtaAtm.restype = None
    lib.getNEP.restype = None

    return lib

def loadTiEMPO2lib_CUDA():
    """!
    Load the TiEMPO2 shared library. Will detect the operating system and link the library accordingly.

    @returns The ctypes library containing the C/C++ functions.
    """

    path_cur = pathlib.Path(__file__).parent.resolve()
    try:
        lib = ctypes.CDLL(os.path.join(path_cur, "libcutiempo2.dll"))
    except:
        try:
            lib = ctypes.CDLL(os.path.join(path_cur, "libcutiempo2.so"))
        except:
            lib = ctypes.CDLL(os.path.join(path_cur, "libcutiempo2.dylib"))

    lib.runTiEMPO2_CUDA.argtypes = [ctypes.POINTER(TStructs.CuInstrument), 
                                    ctypes.POINTER(TStructs.CuTelescope),
                                    ctypes.POINTER(TStructs.CuAtmosphere), 
                                    ctypes.POINTER(TStructs.CuSource),
                                    ctypes.c_int, ctypes.c_char_p]
    
    lib.runTiEMPO2_CUDA.restype = None

    return lib

def runTiEMPO2(instrument, telescope, atmosphere, source, nTimes, nThreads):
    """!
    Binding for running the TiEMPO2 simulation on CPU.

    @param instrument Dictionary containing instrument parameters.
    @param telescope Dictionary containing telescope parameters.
    @param atmosphere Dictionary containing atmosphere parameters.
    @param source Dictionary containing astronomical source parameters.

    @returns 2D array containing timestreams of power in detector, for each channel frequency
    """

    lib = loadTiEMPO2lib()
    mgr = TManager.Manager()

    _instrument = TStructs.Instrument()
    _telescope = TStructs.Telescope()
    _atmosphere = TStructs.Atmosphere()
    _source = TStructs.Source()

    _output = TStructs.Output()

    TBUtils.allfillInstrument(instrument, _instrument)
    TBUtils.allfillTelescope(telescope, _telescope)
    TBUtils.allfillAtmosphere(atmosphere, _atmosphere)
    TBUtils.allfillSource(source, _source)

    cnTimes = ctypes.c_int(nTimes)
    cnThreads = ctypes.c_int(nThreads)

    TBUtils.allocateOutput(_output, nTimes, instrument["nf_ch"])

    args = [_instrument, _telescope, _atmosphere, _source, _output, cnTimes, cnThreads]

    mgr.new_thread(target=lib.runTiEMPO2, args=args)

    res = TBUtils.OutputStructToDict(_output, nTimes, instrument["nf_ch"], np_t=np.float64)

    return res

def calcW2K(instrument, telescope, atmosphere, nPWV, nThreads):
    """!
    Binding for calculating a power-temperature conversion table.

    @param instrument Dictionary containing instrument parameters.
    @param telescope Dictionary containing telescope parameters.
    @param atmosphere Dictionary containing atmosphere parameters.
    @param w2k Watt to Kelvin specification parameters.

    @returns A CalOutput dictionary, containing power-temperature relations.
    """

    lib = loadTiEMPO2lib()
    mgr = TManager.Manager()

    _instrument = TStructs.Instrument()
    _telescope = TStructs.Telescope()
    _atmosphere = TStructs.Atmosphere()

    _caloutput = TStructs.CalOutput()

    TBUtils.allfillInstrument(instrument, _instrument)
    TBUtils.allfillTelescope(telescope, _telescope)
    TBUtils.allfillAtmosphere(atmosphere, _atmosphere)

    TBUtils.allocateCalOutput(_caloutput, nPWV, instrument["nf_ch"])

    cnPWV = ctypes.c_int(nPWV)
    cnThreads = ctypes.c_int(nThreads)

    args = [_instrument, _telescope, _atmosphere, _caloutput, cnPWV, cnThreads]

    mgr.new_thread(target=lib.calcW2K, args=args)

    res = TBUtils.CalOutputStructToDict(_caloutput, nPWV, instrument["nf_ch"], np_t=np.float64)

    return res

def getSourceSignal(instrument, telescope, atmosphere, I_nu, PWV, ON):
    """!
    Binding for calculating the source signal, through the optical path, but without noise.

    @param instrument Dictionary containing instrument parameters.
    @param telescope Dictionary containing telescope parameters.
    @param atmosphere Dictionary containing atmosphere parameters.
    @param source Dictionary containing astronomical source parameters.
    @param Az_point Azimuth point on-sky at which to calculate source intensity.
    @param El_point Elevation point on-sky at which to calculate source intensity.
    @param PWV PWV value of atmosphere, in mm. If empty, no atmosphere used.
    @param ON Use ON-path. If False, uses OFF path.

    @returns 1D array containing power for each detector.
    """

    lib = loadTiEMPO2lib()
    mgr = TManager.Manager()

    _instrument = TStructs.Instrument()
    _telescope = TStructs.Telescope()

    coutput = (ctypes.c_double * instrument["nf_ch"]).from_buffer(np.zeros(instrument["nf_ch"]))

    cI_nu = (ctypes.c_double * I_nu.size)(*(I_nu.ravel().tolist()))

    cPWV = ctypes.c_double(PWV)
    cON = ctypes.c_bool(ON)

    TBUtils.allfillInstrument(instrument, _instrument)
    TBUtils.allfillTelescope(telescope, _telescope)
    
    args = [_instrument, _telescope, coutput, cI_nu, cPWV, cON]

    mgr.new_thread(target=lib.getSourceSignal, args=args)

    res = np.ctypeslib.as_array(coutput, shape=instrument["nf_ch"]).astype(np.float64)
    
    return res

def getChopperCalibration(instrument, Tcal):
    """!
    Binding for calculating the power obtained from a calibration blackbody.

    @param instrument Dictionary containing instrument parameters.
    @param Tcal Temperature of calibration blackbody in Kelvin.

    @returns 1D array containing power for each detector.
    """

    lib = loadTiEMPO2lib()
    mgr = TManager.Manager()

    _instrument = TStructs.Instrument()

    coutput = (ctypes.c_double * instrument["nf_ch"]).from_buffer(np.zeros(instrument["nf_ch"]))

    cTcal = ctypes.c_double(Tcal)

    TBUtils.allfillInstrument(instrument, _instrument)
    
    args = [_instrument, coutput, cTcal]

    mgr.new_thread(target=lib.getChopperCalibration, args=args)

    res = np.ctypeslib.as_array(coutput, shape=instrument["nf_ch"]).astype(np.float64)
    
    return res

def getEtaAtm(instrument, PWV):
    """!
    Binding for running the TiEMPO2 simulation.

    @param source Dictionary containing astronomical source parameters.
    @param atmosphere Dictionary containing atmosphere parameters.
    @param PWV PWV value of atmosphere, in mm.
    
    @returns 1D array containing atmospheric transmission for each detector.
    """


    lib = loadTiEMPO2lib()
    mgr = TManager.Manager()

    _source = TStructs.Source()
    coutput = (ctypes.c_double * instrument["nf_src"])(*(np.zeros(instrument["nf_src"]).tolist()))

    cPWV = ctypes.c_double(PWV)

    f_src_spec = TBUtils.arr2ArrSpec(instrument["f_src"])

    args = [f_src_spec, coutput, cPWV]

    mgr.new_thread(target=lib.getEtaAtm, args=args)

    res = np.ctypeslib.as_array(coutput, shape=instrument["nf_src"]).astype(np.float64)
    
    return res

def getNEP(instrument, telescope, atmosphere, PWV):
    """!
    Binding for running the TiEMPO2 simulation.

    @param instrument Dictionary containing instrument parameters.
    @param telescope Dictionary containing telescope parameters.
    @param atmosphere Dictionary containing atmosphere parameters.
    @param source Dictionary containing astronomical source parameters.
    @param PWV PWV value of atmosphere, in mm.
    
    @returns 1D array containing NEP for each detector.
    """

    lib = loadTiEMPO2lib()
    mgr = TManager.Manager()

    _instrument = TStructs.Instrument()
    _telescope = TStructs.Telescope()
    
    coutput = (ctypes.c_double * instrument["nf_ch"]).from_buffer(np.zeros(instrument["nf_ch"]))

    cPWV = ctypes.c_double(PWV)
    cTatm = ctypes.c_double(atmosphere["Tatm"])

    TBUtils.allfillInstrument(instrument, _instrument)
    TBUtils.allfillTelescope(telescope, _telescope)

    args = [_instrument, _telescope, coutput, cPWV, cTatm]

    mgr.new_thread(target=lib.getNEP, args=args)

    res = np.ctypeslib.as_array(coutput, shape=instrument["nf_ch"]).astype(np.float64)
    
    return res

def runTiEMPO2_CUDA(instrument, telescope, atmosphere, source, nTimes, outpath):
    """!
    Binding for running the TiEMPO2 simulation on GPU.

    @param instrument Dictionary containing instrument parameters.
    @param telescope Dictionary containing telescope parameters.
    @param atmosphere Dictionary containing atmosphere parameters.
    @param source Dictionary containing astronomical source parameters.
    @param simparams Dictionary containing simulation parameters.
    @param outpath Path to directory where TiEMPO2 output is stored.

    @returns 2D array containing timestreams of power in detector, for each channel frequency
    """
    import time


    lib = loadTiEMPO2lib_CUDA()
    mgr = TManager.Manager()

    _instrument = TStructs.CuInstrument()
    _telescope = TStructs.CuTelescope()
    _atmosphere = TStructs.CuAtmosphere()
    _source = TStructs.CuSource()

    _output = TStructs.CuOutput()

    ct_t = ctypes.c_float

    TBUtils.allfillInstrument(instrument, _instrument, ct_t)
    TBUtils.allfillTelescope(telescope, _telescope, ct_t)
    start = time.time()
    TBUtils.allfillAtmosphere(atmosphere, _atmosphere, ct_t, coalesce=True)
    end = time.time()
    TBUtils.allfillSource(source, _source, ct_t)

    cnTimes = ctypes.c_int(nTimes)
    coutpath = ctypes.c_char_p(outpath.encode())

    size_out = nTimes * instrument["nf_ch"]

    timed = end-start

    args = [_instrument, _telescope, _atmosphere, _source, cnTimes, coutpath]

    mgr.new_thread(target=lib.runTiEMPO2_CUDA, args=args)

