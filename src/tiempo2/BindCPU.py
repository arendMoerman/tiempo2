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
    
    lib.runTiEMPO2.restype = None

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

    TBUtils.allocateOutput(_output, instrument.get("freqs_filt").size, ctypes.c_double)

    args = [_instrument, _telescope, _atmosphere, _source, _simparams, _output]

    mgr.new_thread(target=lib.runTiEMPO2, args=args)

    res = TBUtils.OutputStructToDict(_output, instrument.get("freqs_filt").shape, np_t=np.float64)

    return res


