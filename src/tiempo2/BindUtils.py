"""!
@file
File containing utility functions for the ctypes bindings.

Most of these functions are concerned with allocating memory.
"""

import numpy as np
import tiempo2.Structs as TStructs
import ctypes
import os
import array as ar

def allfillInstrument(InstDict, InstStruct, ct_t=ctypes.c_double):
    """!
    Allocate and fill an instrument struct for ctypes.
    
    @param InstDict Dictionary containing instrument parameters.
    @param InstStruct Struct to be filled and passed to ctypes.
    @param ct_t Type of data. Use ctypes.c_double for CPU, ctypes.c_float for GPU.
    """

    arr_t = 'd' if ct_t == ctypes.c_double else 'f'
    arr_filterbank = ar.array(arr_t, InstDict["filterbank"].ravel())

    df_src = (InstDict["f1_src"] - InstDict["f0_src"]) / InstDict["nf_src"]

    InstStruct.nf_ch = ctypes.c_int(InstDict["nf_ch"])

    InstStruct.f_spec = arr2ArrSpec(InstDict["f_src"], ct_t)

    InstStruct.eta_inst = ct_t(InstDict["eta_inst"])
    InstStruct.eta_misc = ct_t(InstDict["eta_misc"])
    InstStruct.f_sample = ct_t(InstDict["f_sample"])
    InstStruct.filterbank = (ct_t * InstDict["filterbank"].size).from_buffer(arr_filterbank)
    InstStruct.delta = ct_t(InstDict["delta"])
    InstStruct.eta_pb = ct_t(InstDict["eta_pb"])

def allfillTelescope(TelDict, TelStruct, ct_t=ctypes.c_double):
    """!
    Allocate and fill a telescope struct for ctypes.
    
    @param TelDict Dictionary containing telescope parameters.
    @param TelStruct Struct to be filled and passed to ctypes.
    @param ct_t Type of data. Use ctypes.c_double for CPU, ctypes.c_float for GPU.
    """
    
    TelStruct.Ttel = ct_t(TelDict["Ttel"])
    TelStruct.Tgnd = ct_t(TelDict["Tgnd"])
    TelStruct.Dtel = ct_t(TelDict["Dtel"])
    TelStruct.chop_mode = ctypes.c_int(TelDict["chop_mode"])
    TelStruct.dAz_chop = ct_t(TelDict["dAz_chop"])
    TelStruct.freq_chop = ct_t(TelDict["freq_chop"])
    TelStruct.freq_nod = ct_t(TelDict["freq_nod"])
    TelStruct.eta_ap_ON = (ct_t * TelDict["eta_ap_ON"].size)(*(TelDict["eta_ap_ON"].ravel().tolist()))
    TelStruct.eta_ap_OFF = (ct_t * TelDict["eta_ap_OFF"].size)(*(TelDict["eta_ap_OFF"].ravel().tolist()))
    TelStruct.eta_mir = ct_t(TelDict["eta_mir"])
    TelStruct.eta_fwd = ct_t(TelDict["eta_fwd"])

    TelStruct.scantype = ctypes.c_int(TelDict["scantype"])
    TelStruct.El0 = ct_t(TelDict["El0"])
    TelStruct.Ax = ct_t(TelDict["Ax"])
    TelStruct.Axmin = ct_t(TelDict["Axmin"])
    TelStruct.Ay = ct_t(TelDict["Ay"])
    TelStruct.Aymin = ct_t(TelDict["Aymin"])
    TelStruct.wx = ct_t(TelDict["wx"])
    TelStruct.wxmin = ct_t(TelDict["wxmin"])
    TelStruct.wy = ct_t(TelDict["wy"])
    TelStruct.wymin = ct_t(TelDict["wymin"])
    TelStruct.phix = ct_t(TelDict["phix"])
    TelStruct.phiy = ct_t(TelDict["phiy"])


def allfillAtmosphere(AtmDict, AtmStruct, ct_t=ctypes.c_double, coalesce=False):
    """!
    Allocate and fill an atmosphere struct for ctypes.
    
    @param AtmDict Dictionary containing atmosphere parameters.
    @param AtmStruct Struct to be filled and passed to ctypes.
    @param ct_t Type of data. Use ctypes.c_double for CPU, ctypes.c_float for GPU.
    """

    AtmStruct.Tatm = ct_t(AtmDict["Tatm"])
    AtmStruct.v_wind = ct_t(AtmDict["v_wind"])
    AtmStruct.h_column = ct_t(AtmDict["h_column"])
    AtmStruct.dx = ct_t(AtmDict["dx"])
    AtmStruct.dy = ct_t(AtmDict["dy"])
    AtmStruct.path = ctypes.c_char_p(os.path.join(AtmDict["path"], "prepd").encode())

def allfillSource(SourceDict, SourceStruct, ct_t=ctypes.c_double):
    """!
    Allocate and fill a source object struct for ctypes.
    
    @param SourceDict Dictionary containing source angular extents and intensity maps.
    @param SourceStruct Struct to be filled and passed to ctypes.
    @param ct_t Type of data. Use ctypes.c_double for CPU, ctypes.c_float for GPU.
    """
       
    I_nu = SourceDict["I_nu"]

    nI_nu = I_nu.ravel().size

    SourceStruct.Az_spec = arr2ArrSpec(SourceDict["Az_src"], ct_t)
    SourceStruct.El_spec = arr2ArrSpec(SourceDict["El_src"], ct_t)
    
    SourceStruct.I_nu = (ct_t * nI_nu)(*(I_nu.ravel().tolist())) 
    SourceStruct.nI_nu = ctypes.c_int(nI_nu)

def allfillSimParams(SPDict, SPStruct, ct_t=ctypes.c_double):
    """!
    Allocate and fill parameters and settings that govern the simulation.

    @param SPDict Dictionary containing simulation parameters.
    @param SPStruct Struct to be filled and passed to ctypes.
    @param ct_t Type of data. Use ctypes.c_double for CPU, ctypes.c_float for GPU.
    """
    
    SPStruct.t_obs = ct_t(SPDict["t_obs"])
    SPStruct.nTimes = ctypes.c_int(SPDict["nTimes"])
    SPStruct.nThreads = ctypes.c_int(SPDict["nThreads"])

def allocateOutput(OutputStruct, size_t, size_f, ct_t=ctypes.c_double):
    """!
    Allocate memory for an output struct.

    @param OutputStruct Struct to be allocated and passed to ctypes.
    @param size_t Number of time evaluations.
    @param size_f Number of channels in filterbank.
    @param ct_t Type of data. Use ctypes.c_double for CPU, ctypes.c_float for GPU.
    """

    fill_sig = np.zeros(size_t * size_f)
    fill_t = np.zeros(size_t)

    OutputStruct.signal = (ct_t * (size_t * size_f)).from_buffer(fill_sig)
    OutputStruct.Az = (ct_t * size_t).from_buffer(fill_t)
    OutputStruct.El = (ct_t * size_t).from_buffer(fill_t)
    OutputStruct.flag = (ctypes.c_int * size_t).from_buffer(fill_t.astype(int))
    
    if isinstance(OutputStruct, TStructs.Output):
        OutputStruct.t_thread = ct_t(0.)

    else:
        OutputStruct.t_diag = (ct_t * 3)(0, 0, 0)

def allocateCalOutput(CalOutputStruct, size_t, size_f, ct_t=ctypes.c_double):
    """!
    Allocate memory for a power-temperature calibration output struct.

    @param CalOutputStruct Struct to be allocated and passed to ctypes.
    @param size_t Number of PWV evaluations.
    @param size_f Number of channels in filterbank.
    @param ct_t Type of data. Use ctypes.c_double for CPU, ctypes.c_float for GPU.
    """

    fill_power = np.zeros(size_t * size_f)
    fill_temperature = np.zeros(size_t * size_f)

    CalOutputStruct.power = (ct_t * (size_t * size_f)).from_buffer(fill_power)
    CalOutputStruct.temperature = (ct_t * (size_t * size_f)).from_buffer(fill_temperature)

def OutputStructToDict(OutputStruct, size_t, size_f, np_t, CPU=True):
    """!
    Convert an output struct to a dictionary.

    @param OutputStruct Struct filled with output.
    @param size_t Number of time evaluations.
    @param size_f Number of channels in filterbank.
    @param np_t Numpy type of array elements.
    """

    sig_shape = (size_f, size_t)
    t_shape = (size_t,)
    
    if CPU:
        times = float(OutputStruct.t_thread)

    else:
        times = np.ctypeslib.as_array(OutputStruct.t_diag, shape=(3,)).astype(np_t)

    OutputDict = {
    "signal" : np.ctypeslib.as_array(OutputStruct.signal, shape=sig_shape).astype(np_t).T,
    "Az" : np.ctypeslib.as_array(OutputStruct.Az, shape=t_shape).astype(np_t),
    "El" : np.ctypeslib.as_array(OutputStruct.El, shape=t_shape).astype(np_t),
    "flag" : np.ctypeslib.as_array(OutputStruct.flag, shape=t_shape).astype(ctypes.c_int),
    "t_diag" : times
    }

    return OutputDict

def CalOutputStructToDict(CalOutputStruct, size_t, size_f, np_t):
    """!
    Convert a calibration output struct to a dictionary.

    @param CalOutputStruct Struct filled with output.
    @param size_t Number of PWV evaluations.
    @param size_f Number of channels in filterbank.
    @param np_t Numpy type of array elements.
    """

    sig_shape = (size_f, size_t)
    t_shape = (size_t,)

    CalOutputDict = {
    "power" : np.ctypeslib.as_array(CalOutputStruct.power, shape=sig_shape).astype(np_t).T,
    "temperature" : np.ctypeslib.as_array(CalOutputStruct.temperature, shape=sig_shape).astype(np_t).T,
    }

    return CalOutputDict

def allfillArrSpec(arr, ArrSpecStruct, ct_t=ctypes.c_double):
    ArrSpecStruct.start = ct_t(arr[0])    
    ArrSpecStruct.step = ct_t(arr[1] - arr[0])    
    ArrSpecStruct.num = ctypes.c_int(arr.size)    

def arr2ArrSpec(arr, ct_t=ctypes.c_double):
    spec_t = TStructs.ArrSpec
    if ct_t == ctypes.c_float:
        spec_t = TStructs.CuArrSpec

    arrspec = spec_t()
    allfillArrSpec(arr, arrspec, ct_t)

    return arrspec





