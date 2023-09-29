/*! \file InterfaceCPU.cpp
    \brief Implementations of library functions for simulations on CPU.

*/

#include "InterfaceCPU.h"

TIEMPO2_DLL void runTiEMPO2(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, SimParams *simparams) {
    std::vector<std::thread> threadPool;
    threadPool.resize(simparams->nThreads);
    // Timesteps
    double dt = 1. / instrument->freq_sample;
    
    // Number of time evaluations. Time is nTimes * dt
    int nTimes = ceil(simparams->t_obs / dt);
    
    // Number of steps per thread
    int step = ceil(nTimes / simparams->nThreads);
    
    int final_step;

    for(int n=0; n < simparams->nThreads; n++) {
        if(n == (simparams->nThreads - 1)) {
            final_step = nTimes;
        }

        else {
            final_step = (n+1) * step;
        }

        threadPool[n] = std::thread(parallelJobs, instrument, telescope, atmosphere, source, n * step, final_step, dt, nTimes, simparams->nThreads);
    }

    // Wait with execution until all threads are done
    for (std::thread &t : threadPool) {
        if (t.joinable()) {
            t.join();
        }
    }
}

TIEMPO2_DLL void parallelJobs(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, int start, int stop, double dt, int nTimes, int nThreads) {
    // Get starting time and chop parameters
    double t, PWV_Gauss_interp, eta_atm_interp, freq, I_nu;
    int n_switch, start_slice, end_slice, total_calcs, size_I_nu, plane_I_nu;

    total_calcs = floor(nTimes / nThreads);
    plane_I_nu = source->nAz * source->nEl;
    size_I_nu = plane_I_nu * source->nfreqs_src;

    // Structs for storing sky and atmosphere co-ordinates
    AzEl center;
    center.Az = 0;
    center.El = 0;

    AzEl pointing;
    xy_atm point_atm;

    bool chop_flag;
    double* slice_I_nu = (double*)malloc(plane_I_nu * sizeof(double));
    
    double* I_atm = (double*)malloc(source->nfreqs_src * sizeof(double));
    double* I_gnd = (double*)malloc(source->nfreqs_src * sizeof(double));
    double* I_tel = (double*)malloc(source->nfreqs_src * sizeof(double));

    // Calculate I_atm, I_gnd, I_tel before entering time loop.
    // These stay constant during observation anyways.
    for(int j=0; j<source->nfreqs_src; j++)
    {
        freq = source->freqs_src[j];
        
        I_atm[j] = getPlanck(atmosphere->Tatm, freq) * SI_TO_MJY; 
        I_gnd[j] = getPlanck(atmosphere->Tatm, freq) * SI_TO_MJY; 
        I_atm[j] = getPlanck(atmosphere->Tatm, freq) * SI_TO_MJY; 
    }


    for(int i=start; i<stop; i++)
    { // Update time t = i * dt;
        n_switch = floor(t * telescope->freq_chop);
        
        // even -> chopper = open
        chop_flag = (n_switch % 2 == 0);
        
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

            start_slice = plane_I_nu * j;
            end_slice = plane_I_nu * (j + 1);

            for(int k=start_slice; k<end_slice; k++)
            {
                slice_I_nu[k - start_slice] = source->I_nu[k];
            }
            
            // Interpolate on source
            I_nu = interpValue(pointing.Az, pointing.El, 
                source->Az, source->El, source->nAz, 
                source->nEl, slice_I_nu);
        }

        //if(i % 10 == 0 and start == 0)
        //{
        //    printf("%d / 100\r", i / total_calcs * 100);
        //}
    }
}


