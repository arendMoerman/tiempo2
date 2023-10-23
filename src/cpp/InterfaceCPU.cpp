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
    int nTimes;     // Number of required time evaluations.

    // Double array types
    double* I_atm = new double[source->nfreqs_src];
    double* I_gnd = new double[source->nfreqs_src];
    double* I_tel = new double[source->nfreqs_src];

    double** ret_container_on = new double*[simparams->nThreads];   // Container for storing on-source power.
    double** ret_container_off_l = new double*[simparams->nThreads];   // Container for storing off-source power, nod A.
    double** ret_container_off_r = new double*[simparams->nThreads];   // Container for storing off-source power, nod B.
    double** slice_container = new double*[simparams->nThreads];    // Container for storing source slices.
    double** P_nu = new double*[simparams->nThreads];    // Container for storing power in each bin.

    int** n_times = new int*[simparams->nThreads]; // Timer for each thread to keep track of on, off (nod A) and off (nod B). 

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
    nTimes = ceil(simparams->t_obs / dt);
    step = ceil(nTimes / simparams->nThreads);

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

        ret_container_on[n] = new double[instrument->nfreqs_filt]();
        ret_container_off_l[n] = new double[instrument->nfreqs_filt]();
        ret_container_off_r[n] = new double[instrument->nfreqs_filt]();
        P_nu[n] = new double[source->nfreqs_src]();

        n_times[n] = new int[3]();
    }
    
    // Main thread spawning loop
    for(int n=0; n < simparams->nThreads; n++) {
        int final_step; // Final step for 
        
        if(n == (simparams->nThreads - 1)) {
            final_step = nTimes;
        }

        else {
            final_step = (n+1) * step;
        }
        
        threadPool[n] = std::thread(&parallelJobs, instrument, telescope, atmosphere, source, simparams, 
                &effs, n * step, final_step, dt, num_AzEl, P_nu[n],
                ret_container_on[n], ret_container_off_l[n], ret_container_off_r[n], slice_container[n],
                n_times[n], I_atm, I_gnd, I_tel);
    }

    // Wait with execution until all threads are done
    for (std::thread &t : threadPool) {
        if (t.joinable()) {
            t.join();
        }
    }
    for(int n=0; n < simparams->nThreads; n++)
    {
        for(int k=0; k<3; k++) {
            output->times[k] += n_times[n][k] * dt;
        }

        for(int j=0; j<instrument->nfreqs_filt; j++)
        {
            output->P_ON[j] += ret_container_on[n][j];
            output->P_OFF_L[j] += ret_container_off_l[n][j];
            output->P_OFF_R[j] += ret_container_off_r[n][j];
            
        }
    }
    
    for(int n=0; n < simparams->nThreads; n++) {
        delete[] ret_container_on[n]; 
        delete[] ret_container_off_l[n]; 
        delete[] ret_container_off_r[n]; 
        
        delete[] slice_container[n];
        delete[] P_nu[n];
        delete[] n_times[n];
    }

    delete[] ret_container_on; 
    delete[] ret_container_off_l; 
    delete[] ret_container_off_r; 
    delete[] slice_container;
    delete[] P_nu;

    delete[] I_atm;
    delete[] I_gnd;
    delete[] I_tel;

    delete[] n_times;
}

TIEMPO2_DLL void parallelJobs(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, SimParams* simparams, 
        Effs* effs, int start, int stop, double dt, int num_AzEl, double* P_nu, 
        double* ret_on, double* ret_off_l, double* ret_off_r, double* slice_container, 
        int* n_times, double* I_atm, double* I_gnd, double* I_tel) {
    
    // Calculate total constant efficiencies
    
    // Get starting time and chop parameters
    double t, PWV_Gauss_interp, eta_atm_interp, freq, I_nu, sigma_k, P_noise_k, Ap;
    double eta_kj; // Filter efficiency for bin j, at channel k
    double sqrt_samp = sqrt(0.5 * instrument->freq_sample); // Constant term needed for noise calculation

    int n_chop, n_nod, start_slice, end_slice;

    double dfreq = source->freqs_src[1] - source->freqs_src[0];
    Ap = M_PI * telescope->Dtel * telescope->Dtel / 4; 

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

    for(int i=start; i<stop; i++)
    { // Update time 
        t = i * dt + simparams->t0;
        n_chop = floor(t * telescope->freq_chop);
        n_nod = floor(t * telescope->freq_nod);
        // even -> chopper = open
        chop_flag = (n_chop % 2 != 0); // If even (false), ON. Odd (true), OFF.
        nod_flag = (n_nod % 2 != 0); // If even (false), AB. Odd (true), BA.
        //printf("%d\n", chop_flag);
        is_in_lower_half = (t - n_nod / telescope->freq_nod) < (1 / telescope->freq_nod / 2);
        //printf("%f, %f, %f, %f, %f\n", t, 1 / telescope->freq_nod / 2, ret_on[0], ret_off_l[0], ret_off_r[0]);        
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

        // Stamp time in correct container
        if (!chop_flag){
            n_times[0]++;
        }
        
        else {
            if (!position) {
                n_times[1]++;
            }

            else {
                n_times[2]++;
            }
        }

        point_atm = convertAnglesToSpatialAtm(pointing, atmosphere->h_column);

        // Add wind to this - currently only along x-axis and pretty manual
        point_atm.xAz = point_atm.xAz + atmosphere->v_wind * t;
        
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
            P_nu[j] = (telescope->eta_ap[j] * eta_atm_interp * effs->eta_tot_chain * I_nu
                + effs->eta_tot_chain * (1 - eta_atm_interp) * I_atm[j] 
                + effs->eta_tot_gnd * I_gnd[j] 
                + effs->eta_tot_mir * I_tel[j]) * CL*CL / (freq*freq) * dfreq;
        }
        
        // In this loop, calculate P_k, NEP_k and noise
        for(int k=0; k<instrument->nfreqs_filt; k++) {
            double P_k = 0; // Initialise each channel to zero, for each timestep
            double NEP_accum = 0;

            // Can loop over bins again, cheap operations this time
            for(int j=0; j<source->nfreqs_src; j++) {   
                freq = source->freqs_src[j];
                eta_kj = instrument->filterbank[k*source->nfreqs_src + j];
                
                NEP_accum += P_nu[j] * eta_kj * (HP * freq + P_nu[j] * eta_kj / dfreq + 2 * instrument->delta / instrument->eta_pb);
                P_k += P_nu[j] * eta_kj;
            }

            sigma_k = sqrt(2 * NEP_accum) * sqrt_samp;
            P_noise_k = drawGaussian(sigma_k);
            //printf("%15e, %15e\n", P_k, P_noise_k);
            //printf("%15e, %15e\n", P_noise_k, sigma_k);
            P_k += P_noise_k;
            //P_k = sqrt(2 * NEP_accum);//P_noise_k;
            
            if (!chop_flag){
                ret_on[k] = P_k;
            }
            
            else {
                if (!position) {
                    ret_off_l[k] = P_k;
                }

                else {
                    ret_off_r[k] = P_k;
                }
            }
        }
    }
}
