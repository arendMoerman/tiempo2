"""!
@file File containing analysis functions that can be applied to simulation results.
These are not meant to be used immediately by users. These functions are used by interface functions and passed to the parallel job dispatcher. 
"""

import numpy as np
from scipy.signal import welch
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from tqdm import tqdm

import tiempo2.FileIO as fio

def avgDirectSubtract(args, result_path, xargs):
    """!
    Apply full time-averaging and direct atmospheric subtraction.
    Only use for simulations of PSWSC (AB or ABBA) observations. If any scanning is involved, use directSubtract function.

    The variance returned by this method is an estimate of the underlying distribution variance, i.e. the NEP.
    In the interface, this is converted to variance of the mean.

    @param chunk_idxs Array containing chunk indices to be processed by this particular instance of avgDirectSubtract.
    @param result_path Path to simulation results.
    """

    select = lambda x,y, flags : np.squeeze(np.argwhere((flags == int(x)) | (flags == int(y))))
    
    chunk_idxs, thread_idx = args
    f_sample, var_method = xargs

    for idx in _iterator(chunk_idxs, thread_idx):
        res_dict = fio.unpack_output(result_path, idx)

        idx_on = select(0, 2, res_dict["flag"])
        idx_off = select(-1, 1, res_dict["flag"])

        chunks_off = _consecutive(idx_off)

        mean_values_off = np.zeros((len(chunks_off), res_dict["signal"].shape[1]))
        mean_index_off = np.zeros(len(chunks_off))
        
        for chu_idx, chu_off in enumerate(chunks_off):
            mean_values_off[chu_idx] = np.nanmean(res_dict["signal"][chu_off,:], axis=0)
            mean_index_off[chu_idx] = np.nanmean(chu_off)
        
        bkg_on = interp1d(mean_index_off, mean_values_off, kind="linear", axis=0, fill_value="extrapolate")(idx_on)    

        on_sub = res_dict["signal"][idx_on,:] - bkg_on
        N_on = on_sub.shape[0]
        
        if var_method == "PSD":
            f_Pxx, Pxx = welch(on_sub, fs=f_sample, axis=0)
            var_term = (N_on - 1) * np.nanmean(Pxx, axis=0) * f_sample

        elif var_method == "TOD":
            var_term = 2 * (N_on - 1) * np.nanvar(on_sub, axis=0) 

        avg_term = (N_on - 1) * np.nanmean(on_sub, axis=0)

        if idx == chunk_idxs[0]:
            N_tot = N_on
            tot_avg = avg_term
            tot_var = var_term

        else:
            N_tot += N_on
            tot_avg += avg_term
            tot_var += var_term

    tot_avg /= (N_tot - len(chunk_idxs))    
    tot_var /= (N_tot - len(chunk_idxs))

    return tot_avg, tot_var, N_tot

def calcSignalPSD(args, result_path, xargs):
    """!
    Calculate signal PSD of a simulation output.

    @param chunk_idxs Array containing chunk indices to be processed by this particular instance of calcSignalPSD.
    @param result_path Path to simulation results.
    @param xargs List of extra arguments to function, in this case sampling frequency and nperseg parameter.
        If required, this could be updated to pass the full Scipy Welch argument list (maybe).

    @returns signal_psd Averaged PSDs over the chunks handled by this particular instance of calcSignalPSD.
    @returns freq_psd Frequencies at which the signal PSD is defined, in Hertz.
    """
    chunk_idxs, thread_idx = args
    f_sample = xargs[0]

    freq_sample, nperseg = xargs

    num_avg = 0
    for idx in _iterator(chunk_idxs, thread_idx):
        res_dict = fio.unpack_output(result_path, idx)
        f_Pxx, Pxx = welch(res_dict["signal"], fs=f_sample, axis=0)

        if idx == chunk_idxs[0]:
            tot_Pxx = Pxx
            out_f_Pxx = f_Pxx 
            num_avg += 1
        
        elif idx == chunk_idxs[-1]:
            # Because it is hard to add arrays of different size, I just discard the last chunk (if there are more than one) and only keep unifomrly sized chunks.
            break
        
        else:
            tot_Pxx += Pxx
            num_avg += 1

    tot_Pxx /= num_avg

    return tot_Pxx, out_f_Pxx

def _consecutive(data, stepsize=1):
    """
    Take numpy array and return list with arrays containing consecutive blocks
    
    @param data Array in which consecutive chunks are to be located.

    @returns List with arrays of consecutive chunks of data as elements.
    """

    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

def _iterator(x, idx_thread):
    return tqdm(x, ncols=100, total=x.size, colour="GREEN") if idx_thread == 0 else x 
