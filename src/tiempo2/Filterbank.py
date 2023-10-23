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
    """

    R = instrumentDict.get("R")
    freqs_filt = instrumentDict.get("freqs_filt")
    freqs_src = sourceDict.get("freqs_src")
    

    filterbank = np.zeros((instrumentDict.get("n_freqs"), freqs_src.size))

    for j, f_j in enumerate(freqs_filt):
        gamma = f_j / (2 * R)
        filterbank[j,:] = 4 / np.pi * instrumentDict.get("eta_filt")[j] * gamma**2 / ((freqs_src - f_j)**2 + gamma**2)

    return filterbank
