"""!
@file
Checker functions for input dictionaries.
"""

import numpy as np


def checkSourceDict(sourceDict):
    SZlist = ["type", "Az", "El", "nAz", "nEl", "Te", "ne0", "beta",
              "v_pec", "thetac", "Da", "freqs_src"]

    atmlist = ["type", "freqs_src"]
    
    loadlist = ["path", "filename"]

    errlist = []

    if sourceDict.get("type") == "SZ":
        checklist = SZlist

    elif sourceDict.get("type") == "atmosphere":
        checklist = atmlist
    
    elif sourceDict.get("type") == "load":
        checklist = loadlist
    
    elif sourceDict.get("type") == "GalSpec":
        return errlist
    else:
        errlist.append("type")
        return errlist

    for key in checklist:
        if sourceDict.get(key) is None:
            errlist.append(key)
    
    return errlist

def checkTelescopeDict(telescopeDict):
    checklist = ["Dtel", "Ttel", "Tgnd", "chop_mode", "eta_ap_ON", "eta_mir",
                     "eta_fwd", "freq_chop", "dAz_chop"]

    checklist_daisy = ["Ax", "Axmin", "Ay", "Aymin", "wx", "wxmin", "wy", "wymin"]

    errlist = []

    for key in checklist:
        if telescopeDict.get(key) is None:
            errlist.append(key)

    if telescopeDict.get("chop_mode") is None:
        telescopeDict["chop_mode"] = 0
    elif telescopeDict.get("chop_mode") == "direct":
        telescopeDict["chop_mode"] = 1
    elif telescopeDict.get("chop_mode") == "abba":
        telescopeDict["chop_mode"] = 2

    if telescopeDict.get("s_rms") is not None:
        telescopeDict["s_rms"] *= 1e-6 # Convert um to m

    if telescopeDict.get("scantype") is None:
        telescopeDict["scantype"] = "point"

    elif telescopeDict.get("scantype") == "daisy":
        if telescopeDict.get("phix") is None:
            telescopeDict["phix"] = 0
        if telescopeDict.get("phiy") is None:
            telescopeDict["phiy"] = 0
        
        for key in checklist_daisy:
            if telescopeDict.get(key) is None:
                errlist.append(key)

        telescopeDict["scantype"] = 1

    elif telescopeDict.get("scantype") == "point":
        for key in checklist_daisy:
            telescopeDict[key] = 0
        telescopeDict["scantype"] = 0


    return errlist

def checkInstrumentDict(instrumentDict):
    checklist = ["freq_0", "R", "n_freqs", "eta_inst", "freq_sample", "eta_filt", "delta", "eta_pb", "box_eq"]

    errlist = []

    if instrumentDict.get("delta") is None:
        instrumentDict["delta"] = 188 * 1.60218e-19 * 1e-6 # delta_Al, in micro eV

    else:
        instrumentDict["delta"] *= 1.60218e-19 * 1e-6

    if instrumentDict.get("eta_pb") is None:
        instrumentDict["eta_pb"] = 0.4

    if instrumentDict.get("box_eq") is None:
        instrumentDict["box_eq"] = True

    for key in checklist:
        if instrumentDict.get(key) is None:
            errlist.append(key)
    
    return errlist

def checkAtmosphereDict(atmosphereDict):
    checklist = ["Tatm", "filename", "path", "dx", "dy", "h_column", "v_wind", "PWV0"]

    errlist = []

    for key in checklist:
        if atmosphereDict.get(key) is None:
            errlist.append(key)
    
    return errlist

def checkObservationDict(observationDict):
    checklist = ["name_sim", "t_obs", "nThreads", "outDir", "OFF_empty", "get_t_diag"]

    errlist = []

    if observationDict.get("t0") is None:
        observationDict["t0"] = 0.

    if observationDict.get("use_noise") is None:
        observationDict["use_noise"] = 1
    
    if observationDict.get("OFF_empty") is None:
        observationDict["OFF_empty"] = 1
    
    if observationDict.get("get_t_diag") is None:
        observationDict["get_t_diag"] = 1

    for key in checklist:
        if observationDict.get(key) is None:
            errlist.append(key)
    
    return errlist
