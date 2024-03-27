"""!
@file
Generate or read filterbank matrix.
"""

import numpy as np

def generateFilterbankFromR(instrumentDict, sourceDict):
    """!
    Generate a Lorentzian filterbank matrix from resolving power R.

    @param instrumentDict Instrument dictionary.
    @param sourceDict Source dictionary.

    @returns filterbank The Lorentzian filterbank.
    """

    R = instrumentDict.get("R")
    freqs_filt = instrumentDict.get("freqs_filt")
    freqs_src = sourceDict.get("f_src")
    
    A = 1

    if instrumentDict["box_eq"]:
        A = 4 / np.pi

    filterbank = np.zeros((instrumentDict.get("n_freqs"), freqs_src.size))

    for j, f_j in enumerate(freqs_filt):
        gamma = f_j / (2 * R)
        filterbank[j,:] = instrumentDict.get("eta_filt")[j] * (A * gamma**2 / ((freqs_src - f_j)**2 + gamma**2))**instrumentDict.get("order")

    return filterbank
