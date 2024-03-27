/*! \file InterfaceCPU.cpp
    \brief Implementations of library functions for simulations on CPU.

*/

#include "InterfaceCPU.h"

double inline getPlanck(double T, double nu)
{
    double prefac = 2 * HP * nu*nu*nu / (CL*CL);
    double dist = 1 / (exp(HP*nu / (KB*T)) - 1); return prefac * dist;
}


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
    double* I_atm = new double[source->nf];
    double* I_gnd = new double[source->nf];
    double* I_tel = new double[source->nf];
    double* I_CMB = new double[source->nf];

    // Initialise constant efficiency struct
    Effs effs;
    effs.eta_tot_chain = instrument->eta_inst * instrument->eta_misc * telescope->eta_fwd * telescope->eta_mir * 0.5;
    effs.eta_tot_gnd = instrument->eta_inst  * instrument->eta_misc * (1 - telescope->eta_fwd) * telescope->eta_mir * 0.5;
    effs.eta_tot_mir = instrument->eta_inst  * instrument->eta_misc * (1 - telescope->eta_mir) * 0.5;

    // Make threadpool
    std::vector<std::thread> threadPool;
    threadPool.resize(simparams->nThreads);
    
    // PREAMBLE
    dt = 1. / instrument->freq_sample;
    step = ceil(simparams->nTimes / simparams->nThreads);
    
    printf("\033[1;32m\r");
    
    // Calculate I_atm, I_gnd, I_tel before entering time loop.
    // These stay constant during observation anyways.
    for(int j=0; j<source->nf; j++) { 
        freq = source->f0 + source->df * j;
        
        I_atm[j] = getPlanck(atmosphere->Tatm, freq); 
        I_gnd[j] = getPlanck(telescope->Tgnd, freq); 
        I_tel[j] = getPlanck(telescope->Ttel, freq);
        I_CMB[j] = getPlanck(2.725, freq);
    }
    
    // Allocate sub-arrays outside of thread loop - safer I guess
    num_AzEl = source->nAz * source->nEl;
   
    Timer timer;

    timer.start();
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
                I_atm, I_gnd, I_tel, I_CMB, n);
    }

    // Wait with execution until all threads are done
    for (std::thread &t : threadPool) {
        if (t.joinable()) {
            t.join();
        }
    }

    timer.stop();
    output->t_thread = timer.get();
    
    printf("\033[0m\n");

    delete[] I_atm;
    delete[] I_gnd;
    delete[] I_tel;
    delete[] I_CMB;
}

TIEMPO2_DLL void parallelJobs(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, SimParams* simparams, Output* output, 
        Effs* effs, int start, int stop, double dt, int num_AzEl, 
        double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx) {
    
    // Calculate total constant efficiencies
    
    // Get starting time and chop parameters
    double t_start; // Time from start of observation.
    double PWV_Gauss_interp; // Interpolated PWV of Gaussian smoothed screen.
    double eta_atm_interp; // Interpolated eta_atm, over frequency and PWV
    double freq; // Bin frequency
    double I_nu; // Specific intensity of source.
    double sigma_k; // Noise per channel.
    double eta_kj; // Filter efficiency for bin j, at channel k
    double sqrt_samp = sqrt(0.5 * instrument->freq_sample); // Constant term needed for noise calculation

    int n_chop, n_nod, start_slice, end_slice;
    
    double Az_src_max = source->Az0 + source->dAz*(source->nAz-1);
    double El_src_max = source->El0 + source->dEl*(source->nEl-1);

    double* PSD_nu = new double[source->nf];

    // Structs for storing sky and atmosphere co-ordinates
    AzEl center;

    center.Az = 0;
    center.El = 0;

    AzEl pointing;
    xy_atm point_atm;

    bool chop_flag;
    int nod_flag;
    
    double is_in_lower_half;

    // Debug utils
    bool debug = false;

    int position; // A (left) = 0, B (right) = 1
                  //
    int total = stop - start; // Total evaluations, for progress monitoring
    double prog_chunk = total / 100; // Divide total in chunks, s.t. about 100 chunks fit in total
    int chunk_count = 0;
    
    std::random_device rd{};
    std::mt19937 geno{rd()};
    
    int start_x0y0, start_x0y1, start_x1y0, start_x1y1;

    int _i_src, _j_src;

    double t_src, u_src;

    for(int i=start; i<stop; i++)
    { // Update time 
        
        if(threadIdx == 0 and chunk_count <= 100) {
            if(i >= chunk_count * prog_chunk) {
                printf("*** Progress: %d / 100 ***\r", chunk_count);
                fflush(stdout);
                chunk_count++;
            }
        }

        t_start = i * dt;

        if(telescope->chop_mode == 0) {getnochop_posflag(t_start, &center, &pointing, telescope, output->flag[i]);}
        else if(telescope->chop_mode == 1) {getONOFF_posflag(t_start, &center, &pointing, telescope, output->flag[i]);}
        else if(telescope->chop_mode == 2) {getABBA_posflag(t_start, &center, &pointing, telescope, output->flag[i]);}
        
        // STORAGE: Add current pointing to output array
        output->Az[i] = pointing.Az;
        output->El[i] = pointing.El;

        _i_src = floor((pointing.Az - source->Az0) / source->dAz);
        _j_src = floor((pointing.El - source->El0) / source->dEl);
        
        start_x0y0 = _i_src * source->nf + _j_src * source->nf * source->nAz;
        start_x1y0 = (_i_src + 1) * source->nf + _j_src * source->nf * source->nAz;
        start_x0y1 = _i_src * source->nf + (_j_src+1) * source->nf * source->nAz;
        start_x1y1 = (_i_src+1) * source->nf + (_j_src+1) * source->nf * source->nAz;

        t_src = (pointing.Az - (source->Az0 + source->dAz*_i_src)) / source->dAz;
        u_src = (pointing.El - (source->El0 + source->dEl*_j_src)) / source->dEl;
        
        convertAnglesToSpatialAtm(&pointing, &point_atm, atmosphere->h_column);
        
        bool offsource = ((pointing.Az < source->Az0) or (pointing.Az > Az_src_max)) or 
                        ((pointing.El < source->El0) or (pointing.El > El_src_max));

        // Add wind to this - currently only along x-axis and pretty manual
        point_atm.xAz = point_atm.xAz + atmosphere->v_wind * t_start;

        // Interpolate on PWV_Gauss
        PWV_Gauss_interp = interpValue(point_atm.xAz, point_atm.yEl, 
                atmosphere->x0, atmosphere->y0, atmosphere->nx, 
                atmosphere->ny, atmosphere->dx, atmosphere->dy, atmosphere->PWV);
      
        double eta_ap;
        // In this loop, calculate power in each bin
        for(int j=0; j<source->nf; j++)
        {   
            freq = source->f0 + source->df * j;
            eta_atm_interp = interpValue(freq, PWV_Gauss_interp, 
                    atmosphere->f0, atmosphere->PWV0, atmosphere->nf, 
                    atmosphere->nPWV, atmosphere->df, atmosphere->dPWV, atmosphere->eta_atm);

            start_slice = num_AzEl * j;
            
            if(offsource) {I_nu = I_CMB[j];}

            else {
                I_nu = (1-t_src)*(1-u_src) * source->I_nu[start_x0y0 + j];
                I_nu += t_src*(1-u_src) * source->I_nu[start_x1y0 + j];
                I_nu += (1-t_src)*u_src * source->I_nu[start_x0y1 + j];
                I_nu += t_src*u_src * source->I_nu[start_x1y1 + j];
                //I_nu = interpValue(pointing.Az, pointing.El, 
                //    source->Az0, source->El0, source->nAz, 
                //    source->nEl, source->dAz, source->dEl, source->I_nu, start_slice, debug);
            }
        
            if(output->flag[i] == 0 or output->flag[i] == -1) {
                eta_ap = telescope->eta_ap_ON[j];
            }

            else {
                eta_ap = telescope->eta_ap_OFF[j];
            }
            
            PSD_nu[j] = eta_ap * eta_atm_interp * effs->eta_tot_chain * I_nu
                + ( effs->eta_tot_chain * (1 - eta_atm_interp) * I_atm[j] 
                + effs->eta_tot_gnd * I_gnd[j] 
                + effs->eta_tot_mir * I_tel[j]) 
                * CL*CL / (freq*freq);
        }
        
        // In this loop, calculate P_k, NEP_k and noise
        for(int k=0; k<instrument->nfreqs_filt; k++) {
            double P_k = 0; // Initialise each channel to zero, for each timestep
            double NEP_accum = 0;

            // Can loop over bins again, cheap operations this time
            for(int j=0; j<source->nf; j++) {   
                freq = source->f0 + source->df * j;
                eta_kj = instrument->filterbank[k*source->nf + j];
                
                NEP_accum += PSD_nu[j] * eta_kj * (HP * freq + PSD_nu[j] * eta_kj + 2 * instrument->delta / instrument->eta_pb);
                P_k += PSD_nu[j] * eta_kj;
            }

            sigma_k = sqrt(2 * NEP_accum * source->df) * sqrt_samp;
            P_k *= source->df;

            std::normal_distribution<double> gg{0., sigma_k};
            if(i == 0){printf("%.12e\n", gg(geno));}
            P_k += gg(geno);
           
            // STORAGE: Add signal to signal array in output
            output->signal[k * simparams->nTimes + i] = P_k; 

        }
    }
    delete[] PSD_nu;
}

TIEMPO2_DLL void getSourceSignal(Instrument *instrument, Telescope *telescope, Source *source, double *output, double *eta_atm, double *freqs_atm, double *PWV_atm, int nfreqs_atm, int nPWV_atm, double Az, double El, double PWV, bool ON) {
    double freq; // Bin frequency
    double I_nu; // Specific intensity of source.
    double eta_kj; // Filter efficiency for bin j, at channel k
    double eta_atm_interp; // Interpolated eta_atm, over frequency and PWV

    double dfreqs_atm = freqs_atm[1] - freqs_atm[0];
    double dPWV_atm = PWV_atm[1] - PWV_atm[0];
    
    double PSD_nu;

    AzEl pointing;
    pointing.Az = Az;
    pointing.El = El;
    
    int _i_src = floor((Az - source->Az0) / source->dAz);
    int _j_src = floor((El - source->El0) / source->dEl);
        
    int start_x0y0 = _i_src * source->nf + _j_src * source->nf * source->nAz;
    int start_x1y0 = (_i_src + 1) * source->nf + _j_src * source->nf * source->nAz;
    int start_x0y1 = _i_src * source->nf + (_j_src+1) * source->nf * source->nAz;
    int start_x1y1 = (_i_src+1) * source->nf + (_j_src+1) * source->nf * source->nAz;

    double t_src = (pointing.Az - Az) / source->dAz;
    double u_src = (pointing.El - El) / source->dEl;
    
    double eta_tot_chain = instrument->eta_inst * instrument->eta_misc * telescope->eta_fwd * telescope->eta_mir * 0.5;
    int num_AzEl = source->nAz * source->nEl;
    double eta_ap;
    // In this loop, calculate power in each bin
    for(int j=0; j<source->nf; j++)
    {   
        freq = source->f0 + source->df * j;
        
        eta_atm_interp = 1.;

        if(PWV > 0) { 
            eta_atm_interp = interpValue(freq, PWV, 
                    freqs_atm[0], PWV_atm[0], nfreqs_atm, 
                    nPWV_atm, dfreqs_atm, dPWV_atm, eta_atm);
        }
        
        I_nu = (1-t_src)*(1-u_src) * source->I_nu[start_x0y0 + j];
        I_nu += t_src*(1-u_src) * source->I_nu[start_x1y0 + j];
        I_nu += (1-t_src)*u_src * source->I_nu[start_x0y1 + j];
        I_nu += t_src*u_src * source->I_nu[start_x1y1 + j];
        
        if(ON) {
            eta_ap = telescope->eta_ap_ON[j];
        }

        else {
            eta_ap = telescope->eta_ap_OFF[j];
        }


        PSD_nu = eta_ap * eta_atm_interp * eta_tot_chain * I_nu;
        for(int k=0; k<instrument->nfreqs_filt; k++) {
            eta_kj = instrument->filterbank[k*source->nf + j];
            output[k] += PSD_nu * eta_kj * source->df; 
        }
    }
}

TIEMPO2_DLL void getEtaAtm(Source *source, double *output, double *eta_atm, double *freqs_atm, double *PWV_atm, int nfreqs_atm, int nPWV_atm, double PWV) {
    double dfreqs_atm = freqs_atm[1] - freqs_atm[0];
    double dPWV_atm = PWV_atm[1] - PWV_atm[0];
    double freq;

    for(int j=0; j<source->nf; j++)
    {   
        freq = source->f0 + source->df * j;
        output[j] = interpValue(freq, PWV, 
                freqs_atm[0], PWV_atm[0], nfreqs_atm, 
                nPWV_atm, dfreqs_atm, dPWV_atm, eta_atm);
    }
}

TIEMPO2_DLL void getNEP(Instrument *instrument, Telescope *telescope, Source *source, double *eta_atm, double *freqs_atm, double *PWV_atm, int nfreqs_atm, int nPWV_atm, double *output, double PWV, double Tatm) {
    // Double array types
    double I_atm;
    double I_gnd;
    double I_tel;

    // Initialise constant efficiency struct
    Effs effs;
    effs.eta_tot_chain = instrument->eta_inst * instrument->eta_misc * telescope->eta_fwd * telescope->eta_mir * 0.5;
    effs.eta_tot_gnd = instrument->eta_inst  * instrument->eta_misc * (1 - telescope->eta_fwd) * telescope->eta_mir * 0.5;
    effs.eta_tot_mir = instrument->eta_inst  * instrument->eta_misc * (1 - telescope->eta_mir) * 0.5;
    
    double eta_atm_interp; // Interpolated eta_atm, over frequency and PWV
    double freq; // Bin frequency
    double eta_kj; // Filter efficiency for bin j, at channel k
    
    double dfreqs_atm = freqs_atm[1] - freqs_atm[0];
    double dPWV_atm = PWV_atm[1] - PWV_atm[0];

    double PSD_back;
    int num_AzEl = source->nAz * source->nEl;
    int start_slice;
    // In this loop, calculate power in each bin
    for(int j=0; j<source->nf; j++)
    {   
        freq = source->f0 + source->df * j;
        
        I_atm = getPlanck(Tatm, freq); 
        I_gnd = getPlanck(telescope->Tgnd, freq); 
        I_tel = getPlanck(telescope->Ttel, freq);
        
        eta_atm_interp = interpValue(freq, PWV, 
                freqs_atm[0], PWV_atm[0], nfreqs_atm, 
                nPWV_atm, dfreqs_atm, dPWV_atm, eta_atm);

        start_slice = num_AzEl * j;
        
        // Currently only uses optimal throughput
        // Store P_nu in array for later use
        PSD_back = (effs.eta_tot_chain * (1 - eta_atm_interp) * I_atm 
            + effs.eta_tot_gnd * I_gnd 
            + effs.eta_tot_mir * I_tel) * CL*CL / (freq*freq);
        
        for(int k=0; k<instrument->nfreqs_filt; k++) {   
            eta_kj = instrument->filterbank[k*source->nf + j];
            
            output[k] += 2 * source->df * PSD_back * eta_kj * (HP * freq + PSD_back * eta_kj + 2 * instrument->delta / instrument->eta_pb);
        }
    }
    
    for(int k=0; k<instrument->nfreqs_filt; k++) {
        output[k] = sqrt(output[k]);
    }
}

