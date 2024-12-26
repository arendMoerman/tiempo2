"""!
@file File containing analysis functions that can be applied to simulation results.
These are not meant to be used immediately by users. These functions are used by interface functions and passed to the parallel job dispatcher. 
"""

import numpy as np
from scipy.signal import welch
from tqdm import tqdm

import tiempo2.FileIO as fio

def calcSignalPSD(chunk_idxs, result_path, xargs):
    """!
    Calculate signal PSD of a simulation output.

    @param chunk_idxs Array containing chunk indices to be processed by this particular instance of calcSignalPSD.
    @param result_path Path to simulation results.
    @param xargs List of extra arguments to function, in this case sampling frequency and nperseg parameter.
        If required, this could be updated to pass the full Scipy Welch argument list (maybe).

    @returns signal_psd Averaged PSDs over the chunks handled by this particular instance of calcSignalPSD.
    @returns freq_psd Frequencies at which the signal PSD is defined, in Hertz.
    """

    freq_sample, nperseg = xargs
    start_chunk = chunk_idxs[0]

    for idx in chunk_idxs:
        res_dict = fio.unpack_output(result_path, idx)
        freq_psd, Pxx = welch(res_dict["signal"], fs=freq_sample, nperseg=nperseg, axis=0)

        if idx == start_chunk:
            signal_psd = Pxx
        else:
            signal_psd += Pxx

    signal_psd /= chunk_idxs.size

    return signal_psd, freq_psd

def rebinSignal(self, output, freqs_old, nbins_add, final_bin=True):
    """!
    Rebin a simulation result into a coarser bin size.
    
    @param final_bin If number of old bins divided by nbins_add is not an integer, wether to rebin final new bin with less bins, or add extra bins to second-to-last bin.
    """

    shape_old  = output.get("signal").shape
    nbins_old = shape_old[1]
    
    if (not isinstance(nbins_add, int)) or (nbins_add == 1):
        self.clog.error(f"Rebin number must be an integer and larger than 1.")
        exit(1)

    if final_bin:
        nbins_new = math.ceil(nbins_old / nbins_add)
    else:
        nbins_new = math.floor(nbins_old / nbins_add)

    signal_new = np.zeros((shape_old[0], nbins_new))
    freqs_new = np.zeros(nbins_new)

    for nbin in range(nbins_new):
        start_bin = nbin * nbins_add
        if nbin == nbins_new - 1:
            signal_new[:,nbin] = np.mean(output.get("signal")[:,start_bin:], axis=1)
            freqs_new[nbin] = np.mean(freqs_old[start_bin:])

        else:
            signal_new[:,nbin] = np.mean(output.get("signal")[:,start_bin:start_bin+nbins_add], axis=1)
            freqs_new[nbin] = np.mean(freqs_old[start_bin:start_bin+nbins_add])

    output_binned = copy.deepcopy(output)
    output_binned["signal"] = signal_new

    return output_binned, freqs_new

def avgDirectSubtract(args, result_path, xargs):
    """!
    Apply full time-averaging and direct atmospheric subtraction.
    Only use for simulations of PSWSC (AB or ABBA) observations. If any scanning is involved, use directSubtract function.

    @param chunk_idxs Array containing chunk indices to be processed by this particular instance of avgDirectSubtract.
    @param result_path Path to simulation results.
    """

    select = lambda x,y, flags : np.squeeze(np.argwhere((flags == int(x)) | (flags == int(y))))
    
    chunk_idxs, thread_idx = args

    start_chunk = chunk_idxs[0]
    import matplotlib.pyplot as pt
    
    iterator = lambda x, idx : tqdm(x, ncols=100, total=x.size, colour="GREEN") if idx == 0 else x 
    for idx in iterator(chunk_idxs, thread_idx):
        res_dict = fio.unpack_output(result_path, idx)

        on = np.nanmean(res_dict["signal"][select(0, 2, res_dict["flag"]), :], axis=0)
        off = np.nanmean(res_dict["signal"][select(-1, 1, res_dict["flag"]), :], axis=0)

        if idx == start_chunk:
            signal_avg = on - off          
        elif idx == chunk_idxs[-1]:
            continue
        else:
            signal_avg += on - off
            #pt.plot(signal_avg, label=f"{idx}")

    #pt.legend()
    #pt.show()

    signal_avg /= chunk_idxs.size

    #if thread_idx == 0:

    return signal_avg
