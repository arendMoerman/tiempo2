"""!
@file
Checker functions for input dictionaries.

Apart from checking, the functions also set default parameters that are not given by the user in the input dictionary.
"""

import numpy as np
import tiempo2.Materials as TMaterials

def checkSourceDict(sourceDict):
    SZlist = ["type", "Az", "El", "nAz", "nEl", "Te", "ne0", "beta",
              "v_pec", "thetac", "Da"]

    backgroundlist = ["type"]
    
    loadlist = ["path", "filename"]

    errlist = []

    if sourceDict.get("type") == "SZ":
        checklist = SZlist

    elif sourceDict.get("type") == "background":
        checklist = backgroundlist
    
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


    if telescopeDict.get("chop_mode") is None:
        telescopeDict["chop_mode"] = 0
    elif telescopeDict.get("chop_mode") == "direct":
        telescopeDict["chop_mode"] = 1
    elif telescopeDict.get("chop_mode") == "abba":
        telescopeDict["chop_mode"] = 2


    if telescopeDict.get("scantype") is None:
        telescopeDict["scantype"] = "point"
        for item in checklist_daisy:
            telescopeDict[item] = 0
            telescopeDict["phix"] = 0
            telescopeDict["phiy"] = 0
        telescopeDict["scantype"] = 0

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

    for key in checklist:
        if telescopeDict.get(key) is None:
            errlist.append(key)


    return errlist

def checkInstrumentDict(instrumentDict):
    checklist = ["material", "f0_ch", "R", "nf_ch", "f0_src", "f1_src", "nf_src", "eta_inst", "f_sample", "eta_filt", "box_eq", "order"]

    errlist = []

    fail = False

    if isinstance(instrumentDict.get("f0_ch"), np.ndarray):
        instrumentDict["nf_ch"] = instrumentDict.get("f0_ch").size

    if instrumentDict.get("box_eq") is None:
        instrumentDict["box_eq"] = True

    if instrumentDict.get("order") is None:
        instrumentDict["order"] = 1

    for key in checklist:
        if instrumentDict.get(key) is None:
            errlist.append(key)
            fail = True
    
    if isinstance(instrumentDict.get("eta_filt"), float):
        instrumentDict["eta_filt"] *= np.ones(instrumentDict.get("nf_ch"))
    
    if instrumentDict.get("material") == "Al_NbTiN":
        instrumentDict["delta"] = TMaterials.Al_NbTiN["delta"]
        instrumentDict["eta_pb"] = TMaterials.Al_NbTiN["eta_pb"]

    else:
        errlist.append("material")

    return errlist

def checkAtmosphereDict(atmosphereDict):
    checklist = ["Tatm", "filename", "path", "dx", "dy", "h_column", "v_wind", "PWV0"]

    errlist = []

    for key in checklist:
        if atmosphereDict.get(key) is None:
            errlist.append(key)
    
    return errlist

def checkObservationDict(observationDict):
    checklist = ["name_sim", "t_obs", "nThreads", "outDir", "get_t_diag"]

    errlist = []

    if observationDict.get("t0") is None:
        observationDict["t0"] = 0.

    if observationDict.get("use_noise") is None:
        observationDict["use_noise"] = 1
    
    if observationDict.get("get_t_diag") is None:
        observationDict["get_t_diag"] = 1
    
    if observationDict.get("nThreads") is None:
        observationDict["nThreads"] = 1

    for key in checklist:
        if observationDict.get(key) is None:
            errlist.append(key)
    
    return errlist

def checkW2KDict(w2kDict):
    checklist = ["f0", "df", "nf", "nPWV", "nThreads"]

    errlist = []

    if w2kDict.get("nPWV") is None:
        w2kDict["nPWV"] = 100
    
    if w2kDict.get("nThreads") is None:
        w2kDict["nThreads"] = 1

    for key in checklist:
        if w2kDict.get(key) is None:
            errlist.append(key)
    
    return errlist
