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

def avgDirectSubtract_spectrum(args, result_path, xargs):
    """!
    Calculate a spectrum by averaging over the time-domain signal.
    
    Atmosphere removal is done by direct subtraction, using the ON-OFF or ABBA scheme.
    Both schemes will work using this single method.

    @param output Output object obtained from a simulation.

    @returns spectrum Averaged and direct-subtracted spectrum.
    """
    
    chunk_idxs, thread_idx = args

    start_chunk = chunk_idxs[0]
    
    iterator = lambda x, idx : tqdm(x, ncols=100, total=x.size, colour="GREEN") if idx == 0 else x 
    for idx in iterator(chunk_idxs, thread_idx):
        res_dict = fio.unpack_output(result_path, idx)

        on = res_dict["signal"][select(0, 2, res_dict["flag"]), :]
        off = res_dict["signal"][select(-1, 1, res_dict["flag"]), :]
   
    # Determine total length of each timestream to not go out of index
    n_samps = output.get("signal").shape[0]

    n_eff = 0

    # We define two counters: 
    #   - pos_old
    #   - pos_now
    # The first keeps track of the chop position of the previous timestamp.
    # The last of the current timestamp. 
    # This is done to determine in which local average the spectrum goes.
    pos_old = 0
    pos_now = 0

    # Container for storing the local average.
    # Local, as in, for storing the average in a single chop position.
    loc_avg = np.zeros(output.get("signal").shape[1])
    
    # Two containers, each for the stored average in a chopping path during an A or B.
    # When pos_now changes w.r.t. pos_old, the local average container is averaged over the total timepoints.
    # The result is either stored in ON_container (0 and 2) or OFF_container (1 and -1).
    ON_container = np.zeros(output.get("signal").shape[1])
    OFF_container = np.zeros(output.get("signal").shape[1])

    # Container for storing the difference between two subsequent ON-OFF averages.
    # This is added to when a full ON-OFF cycle has been performed.
    diff_container = np.zeros(output.get("signal").shape[1])

    # Total average over observation
    tot_avg = np.zeros(output.get("signal").shape[1])

    # Counter for how many points to average over when averaging inside a single chop position.
    n_in_chop = 0

    # Two counters for how many ON and OFF positions have been recorded within a nod position.
    # If n_ON > n_OFF (n_OFF > n_ON) at the end of a nod position, the final ON (OFF) local average will be discarded.
    n_ON = 0
    n_OFF = 0

    # Trackers of old nod position and current one.
    # This is updated and checked by a lambda whenever the pos_now != pos_old
    check_nod = lambda x : "A" if ((x == 0) or (x == 1)) else "B"

    nod_old = "A"
    nod_now = "A"

    # Total number of chopping nods, as a check
    n_nod_tot_ch = 0

    for i in range(n_samps):
        # Start by updating current position
        pos_now = output.get("flag")[i]
       
        nod_now = check_nod(pos_now)
        
        # Just keep storing in loc_avg
        if pos_now == pos_old:
            loc_avg += output.get("signal")[i,:]
            
            n_in_chop += 1

            n_eff += 1

        # Encountered new chop position.
        else:
            
            # Calculate average over loc_avg in previous chop and store in appropriate container.
            # Also update approriate counter
            if ((pos_old == 0) or (pos_old == 2)):
                ON_container += loc_avg / n_in_chop
                n_ON += 1
            
            elif ((pos_old == 1) or (pos_old == -1)): 
                OFF_container += loc_avg / n_in_chop
                n_OFF += 1
            

            # Check if, by updating the numbers, a pair of ON-OFF averages has been obtained
            if n_ON == n_OFF:
                diff_container += ON_container - OFF_container
                ON_container.fill(0)
                OFF_container.fill(0)
                n_nod_tot_ch += 1

            # Check if also a new nod position has been entered by the new chop
            if nod_old != nod_now:
                ON_container.fill(0)
                OFF_container.fill(0)
                n_ON = 0
                n_OFF = 0
                

            # Reset loc_avg to new signal and set n_in_chop to one
            loc_avg = copy.deepcopy(output.get("signal")[i,:])
            #loc_avg = output.get("signal")[i,:]
            n_in_chop = 1
            n_eff += 1
        
            pos_old = pos_now
            nod_old = nod_now

    tot_avg = diff_container / n_nod_tot_ch

    return tot_avg

def avgDirectSubtract(args, result_path, xargs):
    """!
    Apply full time-averaging and direct atmospheric subtraction.
    Only use for simulations of PSWSC (AB or ABBA) observations. If any scanning is involved, use directSubtract function.

    @param chunk_idxs Array containing chunk indices to be processed by this particular instance of avgDirectSubtract.
    @param result_path Path to simulation results.
    """

    select = lambda x,y, flags : np.squeeze(np.argwhere((flags == int(x)) | (flags == int(y))))
    
    chunk_idxs, thread_idx = args
    print(args)

    start_chunk = chunk_idxs[0]
    
    from sys import getsizeof

    iterator = lambda x, idx : tqdm(x, ncols=100, total=x.size, colour="GREEN") if idx == 0 else x 
    for idx in iterator(chunk_idxs, thread_idx):
        res_dict = fio.unpack_output(result_path, idx)
        print(res_dict["signal"].shape, getsizeof(res_dict["signal"]) / 1e9)
        on = res_dict["signal"][select(0, 2, res_dict["flag"]), :]
        off = res_dict["signal"][select(-1, 1, res_dict["flag"]), :]

        N_on = on.shape[0]

        off_avg = np.nanmean(off, axis=0)

        on_sub = on - off_avg

        on_sum = np.nansum(on_sub, axis=0)
        on_var = np.nanvar(on_sub, axis=0)

        if idx == start_chunk:
            N_tot = N_on
            tot_sum = on_sum
            tot_var = on_var

        else:
            N_tot += N_on
            tot_sum += on_sum
            tot_var += on_var
        
    tot_avg = tot_sum / N_tot
    tot_var /= chunk_idxs.size**2

    return tot_sum, tot_var, N_tot
