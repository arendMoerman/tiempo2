"""!
@file
Checker functions for input dictionaries.
"""

import numpy as np


def checkSourceDict(sourceDict):
    SZlist = ["Az", "El", "nAz", "nEl", "Te", "ne0", "beta",
              "v_pec", "rc", "Da", "freqs"]

    loadlist = ["path", "filename"]

    errlist = []

    if sourceDict.get("type") == "SZ":
        checklist = SZlist

    elif sourceDict.get("type") == "load":
        checklist = loadlist

    else:
        errlist.append("type")
        return errlist

    for key in checklist:
        if sourceDict.get(key) is None:
            errlist.append(key)
    
    return errlist

def checkTelescopeDict(telescopeDict):
    checklist = ["Dtel", "Ttel", "Tgnd", "eta_ap", "eta_mir",
                     "eta_fwd", "f_chop", "dAz_chop"]

    errlist = []

    for key in checklist:
        if telescopeDict.get(key) is None:
            errlist.append(key)
    
    return errlist

def checkInstrumentDict(instrumentDict):
    checklist = ["freqs", "R", "eta_inst", "f_sample"]

    errlist = []

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
    checklist = ["name_sim", "t_obs", "nThreads", "outDir"]

    errlist = []

    for key in checklist:
        if observationDict.get(key) is None:
            errlist.append(key)
    
    return errlist
