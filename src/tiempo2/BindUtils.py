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
    InstStruct.nfreqs_filt = ctypes.c_int(InstDict["n_freqs"])
    InstStruct.R = ctypes.c_int(InstDict["R"])
    InstStruct.eta_inst = ctypes.c_double(InstDict["eta_inst"])
    InstStruct.freq_sample = ctypes.c_double(InstDict["freq_sample"])
    InstStruct.filterbank = (ctypes.c_double * InstDict["filterbank"].size)(*(InstDict["filterbank"].ravel().tolist()))
    InstStruct.delta = ctypes.c_double(InstDict["delta"])
    InstStruct.eta_pb = ctypes.c_double(InstDict["eta_pb"])

def allfillTelescope(TelDict, TelStruct):
    """!
    Allocate and fill a telescope struct for ctypes.
    
    @param TelDict Dictionary containing telescope parameters.
    @param TelStruct Struct to be filled and passed to ctypes.
    """
    
    TelStruct.Ttel = ctypes.c_double(TelDict["Ttel"])
    TelStruct.Tgnd = ctypes.c_double(TelDict["Tgnd"])
    TelStruct.Dtel = ctypes.c_double(TelDict["Dtel"])
    TelStruct.chop_mode = ctypes.c_int(TelDict["chop_mode"])
    TelStruct.dAz_chop = ctypes.c_double(TelDict["dAz_chop"])
    TelStruct.freq_chop = ctypes.c_double(TelDict["freq_chop"])
    TelStruct.freq_nod = ctypes.c_double(TelDict["freq_nod"])
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
       
    nAzEl = SourceDict["I_nu"].ravel().size
    
    present = 1
    if SourceDict.get("type") == "atmosphere":
        present = 0
    
    SourceStruct.present = ctypes.c_int(present)
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
    SPStruct.nTimes = ctypes.c_int(SPDict["nTimes"])
    SPStruct.nThreads = ctypes.c_int(SPDict["nThreads"])
    SPStruct.t0 = ctypes.c_double(SPDict["t0"])
    SPStruct.use_noise = ctypes.c_int(SPDict["use_noise"])

def allocateOutput(OutputStruct, size_t, size_f, ct_t):
    """!
    Allocate memory for an output struct.

    @param OutputStruct Struct to be allocated and passed to ctypes.
    @param size_t Number of time evaluations.
    @param size_f Number of channels in filterbank.
    @param ct_t Ctypes type of array.
    """
    
    fill_sig = np.zeros(size_t * size_f)
    fill_t = np.zeros(size_t)

    OutputStruct.signal = (ct_t * (size_t * size_f))(*(fill_sig.tolist()))
    OutputStruct.Az = (ct_t * size_t)(*(fill_t.tolist()))
    OutputStruct.El = (ct_t * size_t)(*(fill_t.tolist()))
    OutputStruct.flag = (ctypes.c_int * size_t)(*(fill_t.astype(int).tolist()))


def OutputStructToDict(OutputStruct, size_t, size_f, np_t):
    """!
    Convert an output struct to a dictionary.

    @param OutputStruct Struct filled with output.
    @param size_t Number of time evaluations.
    @param size_f Number of channels in filterbank.
    @param np_t Numpy type of array elements.
    """

    sig_shape = (size_t, size_f)
    t_shape = (size_t,)

    OutputDict = {
    "signal" : np.ctypeslib.as_array(OutputStruct.signal, shape=sig_shape).astype(np_t),
    "Az" : np.ctypeslib.as_array(OutputStruct.Az, shape=t_shape).astype(np_t),
    "El" : np.ctypeslib.as_array(OutputStruct.El, shape=t_shape).astype(np_t),
    "flag" : np.ctypeslib.as_array(OutputStruct.flag, shape=t_shape).astype(ctypes.c_int),
    }

    return OutputDict
