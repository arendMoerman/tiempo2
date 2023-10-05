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
    double** ret_container_off = new double*[simparams->nThreads];  // Container for storing off-source power.
    double** slice_container = new double*[simparams->nThreads];    // Container for storing source slices

    // Make threadpool
    std::vector<std::thread> threadPool;
    threadPool.resize(simparams->nThreads);
    
    // PREAMBLE
    dt = 1. / instrument->freq_sample;
    nTimes = ceil(simparams->t_obs / dt);
    step = ceil(nTimes / simparams->nThreads);
    num_AzEl = source->nAz * source->nEl;

    // Calculate I_atm, I_gnd, I_tel before entering time loop.
    // These stay constant during observation anyways.
    for(int j=0; j<source->nfreqs_src; j++)
    {
        freq = source->freqs_src[j];
        
        I_atm[j] = getPlanck(atmosphere->Tatm, freq) * SI_TO_MJY; 
        I_gnd[j] = getPlanck(telescope->Tgnd, freq) * SI_TO_MJY; 
        I_tel[j] = getPlanck(telescope->Ttel, freq) * SI_TO_MJY;
    }
    
    // Allocate sub-arrays outside of thread loop - safer I guess
    for(int n=0; n < simparams->nThreads; n++) {
        slice_container[n] = new double[num_AzEl]();
        ret_container_on[n] = new double[instrument->nfreqs_filt]();
        ret_container_off[n] = new double[instrument->nfreqs_filt]();
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

        threadPool[n] = std::thread(&parallelJobs, instrument, telescope, atmosphere, source, 
                n * step, final_step, dt, num_AzEl, 
                ret_container_on[n], ret_container_off[n], slice_container[n],
                I_atm, I_gnd, I_tel);
    }

    // Wait with execution until all threads are done
    for (std::thread &t : threadPool) {
        if (t.joinable()) {
            t.join();
        }
    }

    //for(int j=0; j<source->nfreqs_src; j++)
    //{
    //    printf("%8.7e : %8.7e\n", ret_container_on[0][j], ret_container_off[0][j]);
    //}
    
    for(int n=0; n < simparams->nThreads; n++)
    {
        for(int j=0; j<source->nfreqs_src; j++)
        {
            output->P_on[j] += ret_container_on[n][j];
            output->P_off[j] += ret_container_off[n][j];
        }
    }
    
    for(int n=0; n < simparams->nThreads; n++) {
        delete[] ret_container_on[n]; 
        delete[] ret_container_off[n];
        delete[] slice_container[n];
    }

    delete[] ret_container_on;
    delete[] ret_container_off;
    delete[] slice_container;

    delete[] I_atm;
    delete[] I_gnd;
    delete[] I_tel;
}

TIEMPO2_DLL void parallelJobs(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, 
        int start, int stop, double dt, int num_AzEl, 
        double* ret_on, double* ret_off, double* slice_container, 
        double* I_atm, double* I_gnd, double* I_tel) {
    
    // Get starting time and chop parameters
    double t, PWV_Gauss_interp, eta_atm_interp, freq, I_nu, I_tot, Ap;
    int n_switch, start_slice, end_slice;

    Ap = M_PI * telescope->Dtel * telescope->Dtel / 4; 

    // Structs for storing sky and atmosphere co-ordinates
    AzEl center;
    center.Az = 0;
    center.El = 0;

    AzEl pointing;
    xy_atm point_atm;

    bool chop_flag;

    double eta_tot_chain = instrument->eta_inst * telescope->eta_fwd * telescope->eta_mir * 0.5;
    double eta_tot_gnd = instrument->eta_inst * (1 - telescope->eta_fwd) * telescope->eta_mir * 0.5;
    double eta_tot_mir = instrument->eta_inst * (1 - telescope->eta_mir) * 0.5;

    for(int i=start; i<stop; i++)
    { // Update time 
        t = i * dt;
        n_switch = floor(t * telescope->freq_chop);
        
        // even -> chopper = open
        chop_flag = (n_switch % 2 != 0);
        //printf("%f \n", telescope->freq_chop);
        
        pointing = scanPoint(center, chop_flag, telescope->dAz_chop);
        point_atm = convertAnglesToSpatialAtm(pointing, atmosphere->h_column);

        // Add wind to this - currently only along x-axis and pretty manual
        point_atm.xAz = point_atm.xAz + atmosphere->v_wind * t;
        
        // Interpolate on PWV_Gauss
        PWV_Gauss_interp = interpValue(point_atm.xAz, point_atm.yEl, 
                atmosphere->x_atm, atmosphere->y_atm, atmosphere->nx, 
                atmosphere->ny, atmosphere->PWV);

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
                source->nEl, slice_container);
            
            I_tot = (telescope->eta_ap[j] * eta_atm_interp * eta_tot_chain * I_nu 
                + eta_tot_chain * (1 - eta_atm_interp) * I_atm[j] 
                + eta_tot_gnd * I_gnd[j] 
                + eta_tot_mir * I_tel[j]) * Ap * dt;
            
            for(int k=0; k<instrument->nfreqs_filt; k++)
            {
                if (chop_flag)
                {
                    ret_off[k] += (instrument->filterbank[k*source->nfreqs_src + j] * I_tot);
                }

                else 
                {
                    ret_on[k] += (instrument->filterbank[k*source->nfreqs_src + j] * I_tot);
                }
            }
        }
    }
}

