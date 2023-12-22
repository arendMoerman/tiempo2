/*! \file InterfaceCPU.cpp
    \brief Implementations of library functions for simulations on CPU.

*/

#include "InterfaceCPU.h"

TIEMPO2_DLL void runTiEMPO2(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, 
        SimParams *simparams, Output *output) {
    
    // ALLOCATIONS
    // Doubles 
    double dt;      // Timestep used during observation.
    double freq;    // Frequency, used for initialising background sources.

    // Integers
    int num_AzEl;   // Number of source points in one Az-El slice, for one frequency.
    int step;       // Stepsize for each thread.

    // Double array types
    double* I_atm = new double[source->nfreqs_src];
    double* I_gnd = new double[source->nfreqs_src];
    double* I_tel = new double[source->nfreqs_src];

    double** slice_container = new double*[simparams->nThreads];    // Container for storing source slices.

    // Initialise constant efficiency struct
    Effs effs;
    effs.eta_tot_chain = instrument->eta_inst * telescope->eta_fwd * telescope->eta_mir * 0.5;
    effs.eta_tot_gnd = instrument->eta_inst * (1 - telescope->eta_fwd) * telescope->eta_mir * 0.5;
    effs.eta_tot_mir = instrument->eta_inst * (1 - telescope->eta_mir) * 0.5;

    // Make threadpool
    std::vector<std::thread> threadPool;
    threadPool.resize(simparams->nThreads);
    
    // PREAMBLE
    dt = 1. / instrument->freq_sample;
    step = ceil(simparams->nTimes / simparams->nThreads);
    
    printf("\033[1;32m\r");
    
    // Calculate I_atm, I_gnd, I_tel before entering time loop.
    // These stay constant during observation anyways.
    for(int j=0; j<source->nfreqs_src; j++)
    {
        freq = source->freqs_src[j];
        
        I_atm[j] = getPlanck(atmosphere->Tatm, freq); 
        I_gnd[j] = getPlanck(telescope->Tgnd, freq); 
        I_tel[j] = getPlanck(telescope->Ttel, freq);
    }
    
    // Allocate sub-arrays outside of thread loop - safer I guess
    num_AzEl = source->nAz * source->nEl;
    for(int n=0; n < simparams->nThreads; n++) {
        slice_container[n] = new double[num_AzEl]();
    }
    
    // Main thread spawning loop
    for(int n=0; n < simparams->nThreads; n++) {
        int final_step; // Final step for 
        
        if(n == (simparams->nThreads - 1)) {
            final_step = simparams->nTimes;
        }

        else {
            final_step = (n+1) * step;
        }
        
        threadPool[n] = std::thread(&parallelJobs, instrument, telescope, atmosphere, source, simparams, output, 
                &effs, n * step, final_step, dt, num_AzEl,
                slice_container[n], I_atm, I_gnd, I_tel, n);
    }

    // Wait with execution until all threads are done
    for (std::thread &t : threadPool) {
        if (t.joinable()) {
            t.join();
        }
    }
    printf("\033[0m\n");
    
    for(int n=0; n < simparams->nThreads; n++) {
        delete[] slice_container[n];
    }
    delete[] slice_container;

    delete[] I_atm;
    delete[] I_gnd;
    delete[] I_tel;
}

TIEMPO2_DLL void parallelJobs(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, SimParams* simparams, Output* output, 
        Effs* effs, int start, int stop, double dt, int num_AzEl, double* slice_container, 
        double* I_atm, double* I_gnd, double* I_tel, int threadIdx) {
    
    // Calculate total constant efficiencies
    
    // Get starting time and chop parameters
    double t_real; // Real time.
    double t_start; // Time from start of observation.
    double PWV_Gauss_interp; // Interpolated PWV of Gaussian smoothed screen.
    double eta_atm_interp; // Interpolated eta_atm, over frequency and PWV
    double freq; // Bin frequency
    double I_nu; // Specific intensity of source.
    double sigma_k; // Noise per channel.
    double eta_kj; // Filter efficiency for bin j, at channel k
    double sqrt_samp = sqrt(0.5 * instrument->freq_sample); // Constant term needed for noise calculation

    int n_chop, n_nod, start_slice, end_slice;

    double dfreq = source->freqs_src[1] - source->freqs_src[0];

    double* PSD_nu = new double[source->nfreqs_src];

    // Structs for storing sky and atmosphere co-ordinates
    AzEl center;
    center.Az = 0;
    center.El = 0;

    AzEl pointing;
    xy_atm point_atm;

    bool chop_flag;
    bool nod_flag;
    bool is_in_lower_half;

    // Debug utils
    bool debug = false;

    int position; // A (left) = 0, B (right) = 1
                  //
    int total = stop - start; // Total evaluations, for progress monitoring
    double prog_chunk = total / 100; // Divide total in chunks, s.t. about 100 chunks fit in total
    int chunk_count = 0;

    for(int i=start; i<stop; i++)
    { // Update time 
        
        if(threadIdx == 0 and chunk_count <= 100) {
            if(i >= chunk_count * prog_chunk) {
                printf("*** Progress: %d / 100 ***\r", chunk_count);
                fflush(stdout);
                chunk_count++;
            }
        }

        t_real = i * dt + simparams->t0;
        t_start = i * dt;

        n_chop = floor(t_start * telescope->freq_chop);
        n_nod = floor(t_start * telescope->freq_nod);
        
        chop_flag = (n_chop % 2 != 0); // If even (false), ON. Odd (true), OFF.
        nod_flag = (n_nod % 2 != 0); // If even (false), AB. Odd (true), BA.
        
        is_in_lower_half = (t_start - n_nod / telescope->freq_nod) < (1 / telescope->freq_nod / 2);
        if(nod_flag)
        {
            if(is_in_lower_half) {
                pointing = scanPoint(center, chop_flag, -1 * telescope->dAz_chop);
                position = 1;
            }
            else {
                pointing = scanPoint(center, chop_flag, telescope->dAz_chop);
                position = 0;
            }
        }

        else
        {
            if(is_in_lower_half) {
                pointing = scanPoint(center, chop_flag, telescope->dAz_chop);
                position = 0;
            }
            else {
                pointing = scanPoint(center, chop_flag, -1 * telescope->dAz_chop);
                position = 1;
            }
        }
        
        // STORAGE: Add current pointing to output array
        output->Az[i] = pointing.Az;
        output->El[i] = pointing.El;

        point_atm = convertAnglesToSpatialAtm(pointing, atmosphere->h_column);

        // Add wind to this - currently only along x-axis and pretty manual
        point_atm.xAz = point_atm.xAz + atmosphere->v_wind * t_real;

        // Interpolate on PWV_Gauss
        PWV_Gauss_interp = interpValue(point_atm.xAz, point_atm.yEl, 
                atmosphere->x_atm, atmosphere->y_atm, atmosphere->nx, 
                atmosphere->ny, atmosphere->PWV);
       
        // In this loop, calculate power in each bin
        for(int j=0; j<source->nfreqs_src; j++)
        {   
            freq = source->freqs_src[j];
            eta_atm_interp = interpValue(freq, PWV_Gauss_interp, 
                    atmosphere->freqs_atm, atmosphere->PWV_atm, atmosphere->nfreqs_atm, 
                    atmosphere->nPWV_atm, atmosphere->eta_atm);

            start_slice = num_AzEl * j;
            end_slice = num_AzEl * (j + 1);

            for(int k=start_slice; k<end_slice; k++)
            {
                slice_container[k - start_slice] = source->I_nu[k];
            }
            
            // Interpolate on source
            I_nu = interpValue(pointing.Az, pointing.El, 
                source->Az, source->El, source->nAz, 
                source->nEl, slice_container, debug);
            
            // Currently only uses optimal throughput
            // Store P_nu in array for later use
            PSD_nu[j] = (telescope->eta_ap[j] * eta_atm_interp * effs->eta_tot_chain * I_nu
                + effs->eta_tot_chain * (1 - eta_atm_interp) * I_atm[j] 
                + effs->eta_tot_gnd * I_gnd[j] 
                + effs->eta_tot_mir * I_tel[j]) * CL*CL / (freq*freq);
        }
        
        // In this loop, calculate P_k, NEP_k and noise
        for(int k=0; k<instrument->nfreqs_filt; k++) {
            double P_k = 0; // Initialise each channel to zero, for each timestep
            double NEP_accum = 0;

            // Can loop over bins again, cheap operations this time
            for(int j=0; j<source->nfreqs_src; j++) {   
                freq = source->freqs_src[j];
                eta_kj = instrument->filterbank[k*source->nfreqs_src + j];
                
                NEP_accum += PSD_nu[j] * eta_kj * (HP * freq + PSD_nu[j] * eta_kj + 2 * instrument->delta / instrument->eta_pb);
                P_k += PSD_nu[j] * eta_kj;
            }

            sigma_k = sqrt(2 * NEP_accum * dfreq) * sqrt_samp;
            P_k *= dfreq;

            if(simparams->use_noise) {
                P_k += drawGaussian(sigma_k);
            }
           
            // STORAGE: Add signal to signal array in output
            output->signal[i * instrument->nfreqs_filt + k] = P_k; 

            // STORAGE: Add correct flag to output array
            if (!chop_flag){
                output->flag[i] = 0;
            }
            
            else {
                if (!position) {
                    output->flag[i] = 2;
                }

                else {
                    output->flag[i] = 1;
                }
            }
        }
    }
}
