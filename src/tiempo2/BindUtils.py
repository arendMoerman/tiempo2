"""!
@file
File containing utility functions for the ctypes bindings.

Most of these functions are concerned with allocating memory.
"""

import numpy as np
import ctypes

def allfillInstrument(InstDict, InstStruct):
    """!
    Allocate and fill an instrument struct for ctypes.
    
    @param InstDict Dictionary containing instrument parameters.
    @param InstStruct Struct to be filled and passed to ctypes.
    """
    
    InstStruct.freqs_filt = (ctypes.c_double * InstDict["freqs_filt"].size)(*(InstDict["freqs_filt"].ravel().tolist()))
    InstStruct.nfreqs_filt = ctypes.c_int(InstDict["freqs_filt"].size)
    InstStruct.R = ctypes.c_int(InstDict["R"])
    InstStruct.eta_inst = ctypes.c_double(InstDict["eta_inst"])
    InstStruct.freq_sample = ctypes.c_double(InstDict["freq_sample"])
    InstStruct.filterbank = (ctypes.c_double * InstDict["filterbank"].size)(*(InstDict["filterbank"].ravel().tolist()))

def allfillTelescope(TelDict, TelStruct):
    """!
    Allocate and fill a telescope struct for ctypes.
    
    @param TelDict Dictionary containing telescope parameters.
    @param TelStruct Struct to be filled and passed to ctypes.
    """
    
    TelStruct.Ttel = ctypes.c_double(TelDict["Ttel"])
    TelStruct.Tgnd = ctypes.c_double(TelDict["Tgnd"])
    TelStruct.Dtel = ctypes.c_double(TelDict["Dtel"])
    TelStruct.dAz_chop = ctypes.c_double(TelDict["dAz_chop"])
    TelStruct.freq_chop = ctypes.c_double(TelDict["freq_chop"])
    TelStruct.eta_ap = (ctypes.c_double * TelDict["eta_ap"].size)(*(TelDict["eta_ap"].ravel().tolist()))
    TelStruct.eta_mir = ctypes.c_double(TelDict["eta_mir"])
    TelStruct.eta_fwd = ctypes.c_double(TelDict["eta_fwd"])

def allfillAtmosphere(AtmDict, AtmStruct):
    """!
    Allocate and fill an atmosphere struct for ctypes.
    
    @param AtmDict Dictionary containing atmosphere parameters.
    @param AtmStruct Struct to be filled and passed to ctypes.
    """

    AtmStruct.Tatm = ctypes.c_double(AtmDict["Tatm"])
    AtmStruct.v_wind = ctypes.c_double(AtmDict["v_wind"])
    AtmStruct.h_column = ctypes.c_double(AtmDict["h_column"])
    AtmStruct.x_atm = (ctypes.c_double * AtmDict["x_atm"].size)(*(AtmDict["x_atm"].ravel().tolist()))
    AtmStruct.y_atm = (ctypes.c_double * AtmDict["y_atm"].size)(*(AtmDict["y_atm"].ravel().tolist())) 
    AtmStruct.nx = ctypes.c_int(AtmDict["x_atm"].size)
    AtmStruct.ny = ctypes.c_int(AtmDict["y_atm"].size)
    AtmStruct.PWV = (ctypes.c_double * AtmDict["PWV"].ravel().size)(*(AtmDict["PWV"].ravel().tolist()))
    AtmStruct.freqs_atm = (ctypes.c_double * AtmDict["freqs_atm"].size)(*(AtmDict["freqs_atm"].ravel().tolist()))
    AtmStruct.nfreqs_atm = ctypes.c_int(AtmDict["freqs_atm"].size)
    AtmStruct.PWV_atm = (ctypes.c_double * AtmDict["PWV_atm"].size)(*(AtmDict["PWV_atm"].ravel().tolist())) 
    AtmStruct.nPWV_atm = ctypes.c_int(AtmDict["PWV_atm"].size)
    AtmStruct.eta_atm = (ctypes.c_double * AtmDict["eta_atm"].ravel().size)(*(AtmDict["eta_atm"].ravel().tolist())) 

def allfillSource(SourceDict, SourceStruct):
    """!
    Allocate and fill a source object struct for ctypes.
    
    @param SourceDict Dictionary containing source angular extents and intensity maps.
    @param SourceStruct Struct to be filled and passed to ctypes.
    """
   
    nAzEl = SourceDict["I_nu"].shape[0] * SourceDict["I_nu"].shape[1]

    SourceStruct.Az = (ctypes.c_double * SourceDict["Az"].size)(*(SourceDict["Az"].ravel().tolist())) 
    SourceStruct.nAz = ctypes.c_int(SourceDict["Az"].size) 
    SourceStruct.El = (ctypes.c_double * SourceDict["El"].size)(*(SourceDict["El"].ravel().tolist())) 
    SourceStruct.nEl = ctypes.c_int(SourceDict["El"].size) 
    SourceStruct.I_nu = (ctypes.c_double * SourceDict["I_nu"].ravel().size)(*(SourceDict["I_nu"].ravel().tolist())) 
    SourceStruct.freqs_src = (ctypes.c_double * SourceDict["freqs_src"].ravel().size)(*(SourceDict["freqs_src"].ravel().tolist())) 
    SourceStruct.nfreqs_src = ctypes.c_int(SourceDict["freqs_src"].size) 

def allfillSimParams(SPDict, SPStruct):
    """!
    Allocate and fill parameters and settings that govern the simulation.

    @param SPDict Dictionary containing simulation parameters.
    @param SPStruct Struct to be filled and passed to ctypes.
    """
    SPStruct.t_obs = ctypes.c_double(SPDict["t_obs"])
    SPStruct.nThreads = ctypes.c_int(SPDict["nThreads"])

def allocateOutput(OutputStruct, size, ct_t):
    """!
    Allocate memory for an output struct.

    @param OutputStruct Struct to be allocated and passed to ctypes.
    @param size Size of arrays to be allocated.
    @param ct_t Ctypes type of array.
    """
    fill = np.zeros(size)
    OutputStruct.P_on = (ct_t * size)(*(fill.tolist()))
    OutputStruct.P_off = (ct_t * size)(*(fill.tolist()))

def OutputStructToDict(OutputStruct, shape, np_t):
    """!
    Convert an output struct to a dictionary.

    @param OutputStruct Struct filled with output.
    @param shape Shape of resulting arrays in dictionary.
    @param np_t Numpy type of array elements.
    """

    on = np.ctypeslib.as_array(OutputStruct.P_on, shape=shape).astype(np_t)
    off = np.ctypeslib.as_array(OutputStruct.P_off, shape=shape).astype(np_t)

    OutputDict = {
            "P_on"    : on,
            "P_off"   : off
            }

    return OutputDict
