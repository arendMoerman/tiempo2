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

    lib.runTiEMPO2.argtypes = [ctypes.POINTER(TStructs.Instrument), ctypes.POINTER(TStructs.Telescope),
                            ctypes.POINTER(TStructs.Atmosphere), ctypes.POINTER(TStructs.Source),
                            ctypes.POINTER(TStructs.SimParams), ctypes.POINTER(TStructs.Output)]
    
    lib.getSourceSignal.argtypes = [ctypes.POINTER(TStructs.Instrument), ctypes.POINTER(TStructs.Telescope),
                            ctypes.POINTER(TStructs.Source), ctypes.POINTER(ctypes.c_double),
                            ctypes.c_double, ctypes.c_double]
    
    lib.getNEP.argtypes = [ctypes.POINTER(TStructs.Instrument), ctypes.POINTER(TStructs.Telescope),
                            ctypes.POINTER(TStructs.Atmosphere), ctypes.POINTER(TStructs.Source), 
                            ctypes.POINTER(ctypes.c_double), ctypes.c_double]
    
    lib.runTiEMPO2.restype = None
    lib.getSourceSignal.restype = None
    lib.getNEP.restype = None

    return lib

def loadTiEMPO2lib_CUDA():
    """!
    Load the TiEMPO2 shared library. Will detect the operating system and link the library accordingly.

    @returns lib The ctypes library containing the C/C++ functions.
    """

    path_cur = pathlib.Path(__file__).parent.resolve()
    try:
        lib = ctypes.CDLL(os.path.join(path_cur, "libcutiempo2.dll"))
    except:
        try:
            lib = ctypes.CDLL(os.path.join(path_cur, "libcutiempo2.so"))
        except:
            lib = ctypes.CDLL(os.path.join(path_cur, "libcutiempo2.dylib"))

    lib.runTiEMPO2_CUDA.argtypes = [ctypes.POINTER(TStructs.CuInstrument), ctypes.POINTER(TStructs.CuTelescope),
                            ctypes.POINTER(TStructs.CuAtmosphere), ctypes.POINTER(TStructs.CuSource),
                            ctypes.POINTER(TStructs.CuSimParams), ctypes.POINTER(TStructs.CuOutput)]
    
    lib.runTiEMPO2_CUDA.restype = None

    return lib

def runTiEMPO2(instrument, telescope, atmosphere, source, simparams):
    """!
    Binding for running the TiEMPO2 simulation.

    @param instrument Dictionary containing instrument parameters.
    @param telescope Dictionary containing telescope parameters.
    @param atmosphere Dictionary containing atmosphere parameters.
    @param source Dictionary containing astronomical source parameters.
    @param simparams Dictionary containing simulation parameters.
    """

    lib = loadTiEMPO2lib()
    mgr = TManager.Manager()

    _instrument = TStructs.Instrument()
    _telescope = TStructs.Telescope()
    _atmosphere = TStructs.Atmosphere()
    _source = TStructs.Source()
    _simparams = TStructs.SimParams()

    _output = TStructs.Output()

    TBUtils.allfillInstrument(instrument, _instrument)
    TBUtils.allfillTelescope(telescope, _telescope)
    TBUtils.allfillAtmosphere(atmosphere, _atmosphere)
    TBUtils.allfillSource(source, _source)
    TBUtils.allfillSimParams(simparams, _simparams)

    TBUtils.allocateOutput(_output, simparams["nTimes"], instrument["freqs_filt"].size)

    args = [_instrument, _telescope, _atmosphere, _source, _simparams, _output]

    mgr.new_thread(target=lib.runTiEMPO2, args=args)

    res = TBUtils.OutputStructToDict(_output, simparams["nTimes"], instrument["freqs_filt"].size, np_t=np.float64)

    return res

# BINDINGS FOR DEBUGGING / CHARACTERISATION
def getSourceSignal(instrument, telescope, source, Az_point, El_point):
    """!
    Binding for running the TiEMPO2 simulation.

    @param instrument Dictionary containing instrument parameters.
    @param telescope Dictionary containing telescope parameters.
    @param atmosphere Dictionary containing atmosphere parameters.
    @param source Dictionary containing astronomical source parameters.
    @param simparams Dictionary containing simulation parameters.
    """

    lib = loadTiEMPO2lib()
    mgr = TManager.Manager()

    _instrument = TStructs.Instrument()
    _telescope = TStructs.Telescope()
    _source = TStructs.Source()

    coutput = (ctypes.c_double * instrument["freqs_filt"].size)(*(np.zeros(instrument["freqs_filt"].size).tolist()))

    cAz_point = ctypes.c_double(Az_point)
    cEl_point = ctypes.c_double(El_point)

    TBUtils.allfillInstrument(instrument, _instrument)
    TBUtils.allfillTelescope(telescope, _telescope)
    TBUtils.allfillSource(source, _source)

    args = [_instrument, _telescope, _source, coutput, cAz_point, cEl_point]

    mgr.new_thread(target=lib.getSourceSignal, args=args)

    res = np.ctypeslib.as_array(coutput, shape=instrument["freqs_filt"]).astype(np.float64)
    
    return res

def getNEP(instrument, telescope, atmosphere, source, PWV_value):
    """!
    Binding for running the TiEMPO2 simulation.

    @param instrument Dictionary containing instrument parameters.
    @param telescope Dictionary containing telescope parameters.
    @param atmosphere Dictionary containing atmosphere parameters.
    @param source Dictionary containing astronomical source parameters.
    @param simparams Dictionary containing simulation parameters.
    """

    lib = loadTiEMPO2lib()
    mgr = TManager.Manager()

    _instrument = TStructs.Instrument()
    _telescope = TStructs.Telescope()
    _atmosphere = TStructs.Atmosphere()
    _source = TStructs.Source()

    coutput = (ctypes.c_double * instrument["freqs_filt"].size)(*(np.zeros(instrument["freqs_filt"].size).tolist()))

    cPWV = ctypes.c_double(PWV_value)

    TBUtils.allfillInstrument(instrument, _instrument)
    TBUtils.allfillTelescope(telescope, _telescope)
    TBUtils.allfillAtmosphere(atmosphere, _atmosphere)
    TBUtils.allfillSource(source, _source)

    args = [_instrument, _telescope, _atmosphere, _source, coutput, cPWV]

    mgr.new_thread(target=lib.getNEP, args=args)

    res = np.ctypeslib.as_array(coutput, shape=instrument["freqs_filt"]).astype(np.float64)
    
    return res

def runTiEMPO2_CUDA(instrument, telescope, atmosphere, source, simparams):
    """!
    Binding for running the TiEMPO2 simulation on GPU.

    @param instrument Dictionary containing instrument parameters.
    @param telescope Dictionary containing telescope parameters.
    @param atmosphere Dictionary containing atmosphere parameters.
    @param source Dictionary containing astronomical source parameters.
    @param simparams Dictionary containing simulation parameters.
    """

    lib = loadTiEMPO2lib_CUDA()
    mgr = TManager.Manager()

    _instrument = TStructs.CuInstrument()
    _telescope = TStructs.CuTelescope()
    _atmosphere = TStructs.CuAtmosphere()
    _source = TStructs.CuSource()
    _simparams = TStructs.CuSimParams()

    _output = TStructs.CuOutput()

    ct_t = ctypes.c_float

    TBUtils.allfillInstrument(instrument, _instrument, ct_t)
    TBUtils.allfillTelescope(telescope, _telescope, ct_t)
    TBUtils.allfillAtmosphere(atmosphere, _atmosphere, ct_t)
    TBUtils.allfillSource(source, _source, ct_t)
    TBUtils.allfillSimParams(simparams, _simparams, ct_t)

    TBUtils.allocateOutput(_output, simparams["nTimes"], instrument["freqs_filt"].size, ct_t)

    args = [_instrument, _telescope, _atmosphere, _source, _simparams, _output]

    mgr.new_thread(target=lib.runTiEMPO2_CUDA, args=args)

    res = TBUtils.OutputStructToDict(_output, simparams["nTimes"], instrument["freqs_filt"].size, np_t=np.float64)

    return res


