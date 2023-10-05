"""!
@file
Generate or read filterbank matrix.
"""

import numpy as np

def generateFilterbankFromR(instrumentDict, sourceDict):
    """!
    Generate a filterbank matrix from resolving power R.

    Note that this method fails if the number of frequencies in the source are different from the instrument frequencies.
    Therefore, if freqs_src != freqs_filt, this method will automatically adjust freqs_src to freqs_filt.

    @param instrumentDict Instrument dictionary.
    @param sourceDict Source dictionary.
    """

    R = instrumentDict.get("R")
    freqs_filt = instrumentDict.get("freqs_filt")
    freqs_src = sourceDict.get("freqs_src")

    if freqs_filt.size != freqs_src.size:
        #Place a warning here that we change source frequencies.
        sourceDict["freqs_src"] = freqs_filt
        freqs_src = freqs_filt

    diag = freqs_filt / R * 1e9
    out = np.diag(diag)

    return out
