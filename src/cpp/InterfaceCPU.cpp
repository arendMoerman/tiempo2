/*! \file InterfaceCPU.cpp
    \brief Implementations of library functions for simulations on CPU.

*/

#include "InterfaceCPU.h"

double inline getPlanck(double T, double nu)
{
    double prefac = 2 * HP * nu*nu*nu / (CL*CL);
    double dist = 1 / (exp(HP*nu / (KB*T)) - 1); 
    return prefac * dist;
}



TIEMPO2_DLL void runTiEMPO2(Instrument<double> *instrument, Telescope<double> *telescope, 
            Atmosphere<double> *atmosphere, Source<double> *source, 
            Output<double> *output, int nTimes, int nThreads) {
    
    // ALLOCATIONS
    // Doubles 
    double dt;      // Timestep used during observation.
    double freq;    // Frequency, used for initialising background sources.

    // Integers
    int step;       // Stepsize for each thread.

    int nf_src = instrument->f_spec.num;

    // Double array types
    double* I_atm = new double[nf_src];
    double* I_gnd = new double[nf_src];
    double* I_tel = new double[nf_src];
    double* I_CMB = new double[nf_src];

    // Initialise constant efficiency struct
    Effs<double> effs;
    effs.eta_tot_chain = instrument->eta_inst * instrument->eta_misc * telescope->eta_fwd * telescope->eta_mir * 0.5;
    effs.eta_tot_gnd = instrument->eta_inst  * instrument->eta_misc * (1 - telescope->eta_fwd) * telescope->eta_mir * 0.5;
    effs.eta_tot_mir = instrument->eta_inst  * instrument->eta_misc * (1 - telescope->eta_mir) * 0.5;

    // Make threadpool
    std::vector<std::thread> threadPool;
    threadPool.resize(nThreads);
    
    // PREAMBLE
    dt = 1. / instrument->f_sample;
    step = ceil(nTimes / nThreads);
    
    //printf("\033[1;32m\r");
    
    // Calculate I_atm, I_gnd, I_tel before entering time loop.
    // These stay constant during observation anyways.
    for(int j=0; j<nf_src; j++) { 
        freq = instrument->f_spec.start + instrument->f_spec.step * j;
        
        I_atm[j] = getPlanck(atmosphere->Tatm, freq); 
        I_gnd[j] = getPlanck(telescope->Tgnd, freq); 
        I_tel[j] = getPlanck(telescope->Ttel, freq);
        I_CMB[j] = getPlanck(2.725, freq);
    }
    
    // Allocate sub-arrays outside of thread loop - safer I guess
   
    Timer timer;

    timer.start();
    // Main thread spawning loop
    double *eta_atm;
    ArrSpec<double> PWV_atm;
    ArrSpec<double> f_atm;

    readEtaATM<double, ArrSpec<double>>(&eta_atm, &PWV_atm, &f_atm);

    for(int n=0; n < nThreads; n++) {
        int final_step; // Final step for 
        
        if(n == (nThreads - 1)) {
            final_step = nTimes;
        } else {
            final_step = (n+1) * step;
        }

        if(telescope->scantype == 0 && telescope->chop_mode == 0) {
            threadPool[n] = std::thread(&parallelJobs_1, instrument, 
                                        telescope, atmosphere, 
                                        source, output, 
                                        eta_atm, PWV_atm, f_atm,
                                        &effs, nTimes,
                                        n * step, final_step, dt,
                                        I_atm, I_gnd, I_tel, I_CMB, n);
        }
        
        else if(telescope->scantype == 0 && telescope->chop_mode == 1) {
            threadPool[n] = std::thread(&parallelJobs_2, instrument, 
                                        telescope, atmosphere, 
                                        source, output, 
                                        eta_atm, PWV_atm, f_atm,
                                        &effs, nTimes,
                                        n * step, final_step, dt,
                                        I_atm, I_gnd, I_tel, I_CMB, n);
        }
        
        else if(telescope->scantype == 0 && telescope->chop_mode == 2) {
            threadPool[n] = std::thread(&parallelJobs_3, instrument, 
                                        telescope, atmosphere, 
                                        source, output, 
                                        eta_atm, PWV_atm, f_atm,
                                        &effs, nTimes,
                                        n * step, final_step, dt,
                                        I_atm, I_gnd, I_tel, I_CMB, n);
        }
        
        else {
            threadPool[n] = std::thread(&parallelJobs, instrument, 
                                        telescope, atmosphere, 
                                        source, output, 
                                        eta_atm, PWV_atm, f_atm,
                                        &effs, nTimes,
                                        n * step, final_step, dt,
                                        I_atm, I_gnd, I_tel, I_CMB, n);
        }
    }

    // Wait with execution until all threads are done
    for (std::thread &t : threadPool) {
        if (t.joinable()) {
            t.join();
        }
    }

    timer.stop();
    output->t_thread = timer.get();
    
    //printf("\033[0m\n");

    delete[] I_atm;
    delete[] I_gnd;
    delete[] I_tel;
    delete[] I_CMB;
    delete[] eta_atm;
}


TIEMPO2_DLL void calcW2K(Instrument<double> *instrument, Telescope<double> *telescope, 
            Atmosphere<double> *atmosphere, CalOutput<double> *output, int nPWV, int nThreads) {
    // ALLOCATIONS
    // Doubles 
    double freq;    // Frequency, used for initialising background sources.

    // Integers
    int step;       // Stepsize for each thread.
    int nf_src = instrument->f_spec.num;

    // Double array types
    double* I_atm = new double[nf_src];
    double* I_gnd = new double[nf_src];
    double* I_tel = new double[nf_src];
    double* I_CMB = new double[nf_src];

    // Initialise constant efficiency struct
    Effs<double> effs;
    effs.eta_tot_chain = instrument->eta_inst * instrument->eta_misc * telescope->eta_fwd * telescope->eta_mir * 0.5;
    effs.eta_tot_gnd = instrument->eta_inst  * instrument->eta_misc * (1 - telescope->eta_fwd) * telescope->eta_mir * 0.5;
    effs.eta_tot_mir = instrument->eta_inst  * instrument->eta_misc * (1 - telescope->eta_mir) * 0.5;

    // Make threadpool
    std::vector<std::thread> threadPool;
    threadPool.resize(nThreads);
    
    // PREAMBLE
    step = ceil(nPWV / nThreads);
    
    //printf("\033[1;32m\r");
    
    // Calculate I_atm, I_gnd, I_tel before entering time loop.
    // These stay constant during observation anyways.
    for(int j=0; j<nf_src; j++) { 
        freq = instrument->f_spec.start + instrument->f_spec.step * j;
        
        I_atm[j] = getPlanck(atmosphere->Tatm, freq); 
        I_gnd[j] = getPlanck(telescope->Tgnd, freq); 
        I_tel[j] = getPlanck(telescope->Ttel, freq);
        I_CMB[j] = getPlanck(2.725, freq);
    }
    
    // Allocate sub-arrays outside of thread loop - safer I guess
    Timer timer;

    timer.start();
    
    double *eta_atm;
    ArrSpec<double> PWV_atm;
    ArrSpec<double> f_atm;

    readEtaATM<double, ArrSpec<double>>(&eta_atm, &PWV_atm, &f_atm);

    double dPWV_arr = (PWV_atm.num * PWV_atm.step - PWV_atm.start ) / nPWV;

    // Main thread spawning loop
    for(int n=0; n < nThreads; n++) {
        int final_step; // Final step for 
        
        if(n == (nThreads - 1)) {
            final_step = nPWV;
        } else {
            final_step = (n+1) * step;
        }
        
        threadPool[n] = std::thread(&parallelJobsW2K, instrument, atmosphere, 
                output, eta_atm, PWV_atm, f_atm, &effs, nPWV, n * step, final_step, dPWV_arr,
                I_atm, I_gnd, I_tel, I_CMB, n);
    }

    // Wait with execution until all threads are done
    for (std::thread &t : threadPool) {
        if (t.joinable()) {
            t.join();
        }
    }

    timer.stop();
    //output->t_thread = timer.get();
    
    //printf("\033[0m\n");

    delete[] I_atm;
    delete[] I_gnd;
    delete[] I_tel;
    delete[] I_CMB;
    delete[] eta_atm;
}


TIEMPO2_DLL void getSourceSignal(Instrument<double> *instrument, Telescope<double> *telescope, 
            double *output, double *I_nu, double PWV, bool ON) {
    double freq; // Bin frequency
    double eta_kj; // Filter efficiency for bin j, at channel k
    double eta_atm_interp; // Interpolated eta_atm, over frequency and PWV
    
    double PSD_nu;
    
    double eta_tot_chain = instrument->eta_inst * instrument->eta_misc * telescope->eta_fwd * telescope->eta_mir * 0.5;
    double eta_ap;
    
    double *eta_atm;
    ArrSpec<double> PWV_atm;
    ArrSpec<double> f_atm;

    readEtaATM<double, ArrSpec<double>>(&eta_atm, &PWV_atm, &f_atm);
    
    for(int j=0; j<instrument->f_spec.num; j++) { 
        freq = instrument->f_spec.start + instrument->f_spec.step * j;
        
        eta_atm_interp = 1.;

        if(PWV > 0) { 
            eta_atm_interp = interpValue(PWV, freq, 
                    PWV_atm, f_atm, eta_atm);
        }
        
        if(ON) {
            eta_ap = telescope->eta_ap_ON[j];
        }

        else {
            eta_ap = telescope->eta_ap_OFF[j];
        }


        PSD_nu = eta_ap * eta_atm_interp * eta_tot_chain * I_nu[j];
        for(int k=0; k<instrument->nf_ch; k++) {
            eta_kj = instrument->filterbank[k*instrument->f_spec.num + j];
            output[k] += PSD_nu * eta_kj * instrument->f_spec.step; 
        }
    }
    delete[] eta_atm;
}

TIEMPO2_DLL void getEtaAtm(ArrSpec<double> f_src, double *output, double PWV) {
    double freq;

    double *eta_atm;
    ArrSpec<double> PWV_atm;
    ArrSpec<double> f_atm;


    readEtaATM<double, ArrSpec<double>>(&eta_atm, &PWV_atm, &f_atm);

    for(int j=0; j<f_src.num; j++)
    {   
        freq = f_src.start + f_src.step * j;
        output[j] = interpValue(PWV, freq, 
                PWV_atm, f_atm, eta_atm);
    }

    delete[] eta_atm;
}

TIEMPO2_DLL void getNEP(Instrument<double> *instrument, Telescope<double> *telescope, 
            double *output, double PWV, double Tatm) {
    // Double array types
    double I_atm;
    double I_gnd;
    double I_tel;

    // Initialise constant efficiency struct
    Effs<double> effs;
    effs.eta_tot_chain = instrument->eta_inst * instrument->eta_misc * telescope->eta_fwd * telescope->eta_mir * 0.5;
    effs.eta_tot_gnd = instrument->eta_inst  * instrument->eta_misc * (1 - telescope->eta_fwd) * telescope->eta_mir * 0.5;
    effs.eta_tot_mir = instrument->eta_inst  * instrument->eta_misc * (1 - telescope->eta_mir) * 0.5;
    
    double eta_atm_interp; // Interpolated eta_atm, over frequency and PWV
    double freq; // Bin frequency
    double eta_kj; // Filter efficiency for bin j, at channel k

    double PSD_back;
    
    double *eta_atm;
    ArrSpec<double> PWV_atm;
    ArrSpec<double> f_atm;

    readEtaATM<double, ArrSpec<double>>(&eta_atm, &PWV_atm, &f_atm);

    for(int j=0; j<instrument->f_spec.num; j++) { 
        freq = instrument->f_spec.start + instrument->f_spec.step * j;
        
        I_atm = getPlanck(Tatm, freq); 
        I_gnd = getPlanck(telescope->Tgnd, freq); 
        I_tel = getPlanck(telescope->Ttel, freq);
        
        eta_atm_interp = interpValue(PWV, freq, 
                PWV_atm, f_atm, eta_atm);
        
        PSD_back = (effs.eta_tot_chain * (1 - eta_atm_interp) * I_atm 
            + effs.eta_tot_gnd * I_gnd 
            + effs.eta_tot_mir * I_tel) * CL*CL / (freq*freq);
        
        for(int k=0; k<instrument->nf_ch; k++) {   
            eta_kj = instrument->filterbank[k*instrument->f_spec.num + j];
            
            output[k] += 2 * instrument->f_spec.step * PSD_back * eta_kj * (HP * freq + PSD_back * eta_kj + 2 * instrument->delta / instrument->eta_pb);
        }
    }
    
    for(int k=0; k<instrument->nf_ch; k++) {
        output[k] = sqrt(output[k]);
    }
    delete[] eta_atm;
}

void parallelJobs_1(Instrument<double> *instrument, Telescope<double> *telescope, 
                  Atmosphere<double> *atmosphere, Source<double> *source, 
                  Output<double> *output, double *eta_atm,
                  ArrSpec<double> PWV_atm, ArrSpec<double> f_atm,
                  Effs<double> *effs, int nTimes,
                  int start, int stop, double dt, 
                  double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx) {
    
    // Calculate total constant efficiencies
    
    // Get starting time and chop parameters
    double t_start; // Time from start of observation.
    double PWV_Gauss_interp; // Interpolated PWV of Gaussian smoothed screen.
    double eta_atm_interp; // Interpolated eta_atm, over frequency and PWV
    double freq; // Bin frequency
    double I_nu; // Specific intensity of source.
    
    struct ArrSpec<double> _f_spec = instrument->f_spec;

    double* PSD_nu = new double[_f_spec.num];

    // Structs for storing sky and atmosphere co-ordinates
    AzEl center;

    center.Az = 0;
    center.El = 0;

    xy_atm point_atm;
    
    int total = stop - start; // Total evaluations, for progress monitoring
    double prog_chunk = total / 100; // Divide total in chunks, s.t. about 100 chunks fit in total
    int chunk_count = 0;
    
    std::random_device rd{};
    std::mt19937 geno{rd()};
    
    for(int i=start; i<stop; i++) { // Update time 
        //if(threadIdx == 0 and chunk_count <= 100) {
        //    if(i >= chunk_count * prog_chunk) {
        //        printf("*** Progress: %d / 100 ***\r", chunk_count);
        //        fflush(stdout);
        //        chunk_count++;
        //    }
        //}

        t_start = i * dt;
        
        // STORAGE: Add current pointing to output array
        output->Az[i] = center.Az;
        output->El[i] = center.El;
        
        convertAnglesToSpatialAtm(&center, &point_atm, atmosphere->h_column);

        // Add wind to this - currently only along x-axis and pretty manual
        point_atm.xAz = point_atm.xAz + atmosphere->v_wind * t_start;

        // Interpolate on PWV_Gauss
        PWV_Gauss_interp = interpValue(point_atm.xAz, point_atm.yEl, 
                atmosphere->x_spec, atmosphere->y_spec, atmosphere->PWV);
        
        output->flag[i] = 0;
      
        // In this loop, calculate power in each bin
        for(int j=0; j<_f_spec.num; j++)
        {   
            freq = _f_spec.start + _f_spec.step * j;
            eta_atm_interp = interpValue(PWV_Gauss_interp, freq, 
                    PWV_atm, f_atm, eta_atm);
        
            I_nu = source->I_nu[j];

            PSD_nu[j] = telescope->eta_ap_ON[j] * eta_atm_interp * effs->eta_tot_chain * I_nu
                + ( effs->eta_tot_chain * (1 - eta_atm_interp) * I_atm[j] 
                + effs->eta_tot_gnd * I_gnd[j] 
                + effs->eta_tot_mir * I_tel[j]) 
                * CL*CL / (freq*freq);
        }
        calcPhotonNoise(instrument, PSD_nu, geno, output, i, nTimes);
    }
    delete[] PSD_nu;
}

void parallelJobs_2(Instrument<double> *instrument, Telescope<double> *telescope, 
                  Atmosphere<double> *atmosphere, Source<double> *source, 
                  Output<double> *output, double *eta_atm,
                  ArrSpec<double> PWV_atm, ArrSpec<double> f_atm,
                  Effs<double> *effs, int nTimes,
                  int start, int stop, double dt, 
                  double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx) {
    
    // Calculate total constant efficiencies
    
    // Get starting time and chop parameters
    double t_start; // Time from start of observation.
    double PWV_Gauss_interp; // Interpolated PWV of Gaussian smoothed screen.
    double eta_atm_interp; // Interpolated eta_atm, over frequency and PWV
    double freq; // Bin frequency
    double I_nu; // Specific intensity of source.

    int n_chop;
    
    struct ArrSpec<double> _f_spec = instrument->f_spec;

    double* PSD_nu = new double[_f_spec.num];

    // Structs for storing sky and atmosphere co-ordinates
    AzEl center;

    center.Az = 0;
    center.El = 0;

    AzEl pointing;
    xy_atm point_atm;

    bool chop_flag;

    int total = stop - start; // Total evaluations, for progress monitoring
    double prog_chunk = total / 100; // Divide total in chunks, s.t. about 100 chunks fit in total
    int chunk_count = 0;
    
    std::random_device rd{};
    std::mt19937 geno{rd()};
    
    int _i_src, _j_src;

    double t_src, u_src;
    for(int i=start; i<stop; i++) { // Update time 
        //if(threadIdx == 0 and chunk_count <= 100) {
        //    if(i >= chunk_count * prog_chunk) {
        //        printf("*** Progress: %d / 100 ***\r", chunk_count);
        //        fflush(stdout);
        //        chunk_count++;
        //    }
        //}

        t_start = i * dt;

        getONOFF_posflag(t_start, &center, &pointing, telescope, output->flag[i]);
        
        // STORAGE: Add current pointing to output array
        output->Az[i] = pointing.Az;
        output->El[i] = pointing.El;
        
        convertAnglesToSpatialAtm(&pointing, &point_atm, atmosphere->h_column);
        
        // Add wind to this - currently only along x-axis and pretty manual
        point_atm.xAz = point_atm.xAz + atmosphere->v_wind * t_start;

        // Interpolate on PWV_Gauss
        PWV_Gauss_interp = interpValue(point_atm.xAz, point_atm.yEl, 
                atmosphere->x_spec, atmosphere->y_spec, atmosphere->PWV);
      
        double eta_ap;
        // In this loop, calculate power in each bin
        for(int j=0; j<_f_spec.num; j++)
        {   
            freq = _f_spec.start + _f_spec.step * j;
            eta_atm_interp = interpValue(PWV_Gauss_interp, freq, 
                    PWV_atm, f_atm, eta_atm);

            I_nu = source->I_nu[output->flag[i] * _f_spec.num + j];
        
            if(output->flag[i] == 0) {
                eta_ap = telescope->eta_ap_ON[j];
            } else {
                eta_ap = telescope->eta_ap_OFF[j];
            }
            
            PSD_nu[j] = eta_ap * eta_atm_interp * effs->eta_tot_chain * I_nu
                + ( effs->eta_tot_chain * (1 - eta_atm_interp) * I_atm[j] 
                + effs->eta_tot_gnd * I_gnd[j] 
                + effs->eta_tot_mir * I_tel[j]) 
                * CL*CL / (freq*freq);
        }
        calcPhotonNoise(instrument, PSD_nu, geno, output, i, nTimes);
    }
    delete[] PSD_nu;
}

void parallelJobs_3(Instrument<double> *instrument, Telescope<double> *telescope, 
                  Atmosphere<double> *atmosphere, Source<double> *source, 
                  Output<double> *output, double *eta_atm,
                  ArrSpec<double> PWV_atm, ArrSpec<double> f_atm,
                  Effs<double> *effs, int nTimes,
                  int start, int stop, double dt, 
                  double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx) {
    
    // Calculate total constant efficiencies
    
    // Get starting time and chop parameters
    double t_start; // Time from start of observation.
    double PWV_Gauss_interp; // Interpolated PWV of Gaussian smoothed screen.
    double eta_atm_interp; // Interpolated eta_atm, over frequency and PWV
    double freq; // Bin frequency
    double I_nu; // Specific intensity of source.

    int n_chop, n_nod;
    
    struct ArrSpec<double> _f_spec = instrument->f_spec;
    
    double* PSD_nu = new double[_f_spec.num];

    // Structs for storing sky and atmosphere co-ordinates
    AzEl center;

    center.Az = 0;
    center.El = 0;

    AzEl pointing;
    xy_atm point_atm;

    bool chop_flag;
    int nod_flag;
    
    double is_in_lower_half;

    int position; // A (left) = 0, B (right) = 1
                  //
    int total = stop - start; // Total evaluations, for progress monitoring
    double prog_chunk = total / 100; // Divide total in chunks, s.t. about 100 chunks fit in total
    int chunk_count = 0;
    
    std::random_device rd{};
    std::mt19937 geno{rd()};
    
    for(int i=start; i<stop; i++) { // Update time 
        //if(threadIdx == 0 and chunk_count <= 100) {
        //    if(i >= chunk_count * prog_chunk) {
        //        printf("*** Progress: %d / 100 ***\r", chunk_count);
        //        fflush(stdout);
        //        chunk_count++;
        //    }
        //}

        t_start = i * dt;

        getABBA_posflag(t_start, &center, &pointing, telescope, output->flag[i]);
        
        if (output->flag[i] == 0 or output->flag[i] == 2) {position = 1;}
        else if (output->flag[i] == 1) {position = 2;}
        else {position = 0;}

        // STORAGE: Add current pointing to output array
        output->Az[i] = pointing.Az;
        output->El[i] = pointing.El;

        convertAnglesToSpatialAtm(&pointing, &point_atm, atmosphere->h_column);
        
        // Add wind to this - currently only along x-axis and pretty manual
        point_atm.xAz = point_atm.xAz + atmosphere->v_wind * t_start;

        // Interpolate on PWV_Gauss
        PWV_Gauss_interp = interpValue(point_atm.xAz, point_atm.yEl, 
                atmosphere->x_spec, atmosphere->y_spec, atmosphere->PWV);
      
        double eta_ap;
        // In this loop, calculate power in each bin
        for(int j=0; j<_f_spec.num; j++)
        {   
            freq = _f_spec.start + _f_spec.step * j;
            eta_atm_interp = interpValue(PWV_Gauss_interp, freq, 
                    PWV_atm, f_atm, eta_atm);

            I_nu = source->I_nu[position * _f_spec.num + j];
        
            if(output->flag[i] == 0 or output->flag[i] == -1) {
                eta_ap = telescope->eta_ap_ON[j];
            } else {
                eta_ap = telescope->eta_ap_OFF[j];
            }
            
            PSD_nu[j] = eta_ap * eta_atm_interp * effs->eta_tot_chain * I_nu
                + ( effs->eta_tot_chain * (1 - eta_atm_interp) * I_atm[j] 
                + effs->eta_tot_gnd * I_gnd[j] 
                + effs->eta_tot_mir * I_tel[j]) 
                * CL*CL / (freq*freq);
        }
        calcPhotonNoise(instrument, PSD_nu, geno, output, i, nTimes);
    }
    delete[] PSD_nu;
}

void parallelJobs(Instrument<double> *instrument, Telescope<double> *telescope, 
                  Atmosphere<double> *atmosphere, Source<double> *source, 
                  Output<double> *output, double *eta_atm,
                  ArrSpec<double> PWV_atm, ArrSpec<double> f_atm,
                  Effs<double> *effs, int nTimes,
                  int start, int stop, double dt, 
                  double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx) {
    
    int num_AzEl;   // Number of source points in one Az-El slice, for one frequency.
    
    // Get starting time and chop parameters
    double t_start; // Time from start of observation.
    double PWV_Gauss_interp; // Interpolated PWV of Gaussian smoothed screen.
    double eta_atm_interp; // Interpolated eta_atm, over frequency and PWV
    double freq; // Bin frequency
    double I_nu; // Specific intensity of source.

    int n_chop, n_nod, start_slice, end_slice;
    
    struct ArrSpec<double> _f_spec = instrument->f_spec;
    struct ArrSpec<double> _Az_src = source->Az_spec;
    struct ArrSpec<double> _El_src = source->El_spec;
    
    double Az_src_max = _Az_src.start + _Az_src.step*(_Az_src.num - 1);
    double El_src_max = _El_src.start + _El_src.step*(_El_src.num - 1);

    double* PSD_nu = new double[_f_spec.num];
    
    num_AzEl = source->Az_spec.num * source->El_spec.num;

    // Structs for storing sky and atmosphere co-ordinates
    AzEl center;

    center.Az = 0;
    center.El = 0;

    AzEl pointing;
    xy_atm point_atm;

    bool chop_flag;
    int nod_flag;
    
    double is_in_lower_half;

    int total = stop - start; // Total evaluations, for progress monitoring
    double prog_chunk = total / 100; // Divide total in chunks, s.t. about 100 chunks fit in total
    int chunk_count = 0;
    
    std::random_device rd{};
    std::mt19937 geno{rd()};
    
    int start_x0y0, start_x0y1, start_x1y0, start_x1y1;

    int _i_src, _j_src;

    double t_src, u_src;
    for(int i=start; i<stop; i++) { // Update time 
        //if(threadIdx == 0 and chunk_count <= 100) {
        //   if(i >= chunk_count * prog_chunk) {
        //        printf("*** Progress: %d / 100 ***\r", chunk_count);
        //        fflush(stdout);
        //        chunk_count++;
        //    }
        //}

        t_start = i * dt;

        if(telescope->chop_mode == 0) {getnochop_posflag(t_start, &center, &pointing, telescope, output->flag[i]);}
        else if(telescope->chop_mode == 1) {getONOFF_posflag(t_start, &center, &pointing, telescope, output->flag[i]);}
        else if(telescope->chop_mode == 2) {getABBA_posflag(t_start, &center, &pointing, telescope, output->flag[i]);}
        
        // STORAGE: Add current pointing to output array
        output->Az[i] = pointing.Az;
        output->El[i] = pointing.El;

        _i_src = floor((pointing.Az - _Az_src.start) / _Az_src.step);
        _j_src = floor((pointing.El - _El_src.start) / _El_src.step);
        
        start_x0y0 = _f_spec.num * (_i_src + _j_src * _Az_src.num);
        start_x1y0 = _f_spec.num * (_i_src + 1 + _j_src * _Az_src.num);
        start_x0y1 = _f_spec.num * (_i_src + (_j_src+1) * _Az_src.num);
        start_x1y1 = _f_spec.num * (_i_src + 1 + (_j_src+1) * _Az_src.num);

        t_src = (pointing.Az - (_Az_src.start + _Az_src.step*_i_src)) / _Az_src.step;
        u_src = (pointing.El - (_El_src.start + _El_src.step*_j_src)) / _El_src.step;
        
        convertAnglesToSpatialAtm(&pointing, &point_atm, atmosphere->h_column);
        
        bool offsource = ((pointing.Az < _Az_src.start) or (pointing.Az > Az_src_max)) or 
                        ((pointing.El < _El_src.start) or (pointing.El > El_src_max));

        // Add wind to this - currently only along x-axis and pretty manual
        point_atm.xAz = point_atm.xAz + atmosphere->v_wind * t_start;

        // Interpolate on PWV_Gauss
        PWV_Gauss_interp = interpValue(point_atm.xAz, point_atm.yEl, 
                atmosphere->x_spec, atmosphere->y_spec, atmosphere->PWV);
      
        double eta_ap;
        // In this loop, calculate power in each bin
        for(int j=0; j<_f_spec.num; j++)
        {   
            freq = _f_spec.start + _f_spec.step * j;
            eta_atm_interp = interpValue(PWV_Gauss_interp, freq, 
                    PWV_atm, f_atm, eta_atm);

            start_slice = num_AzEl * j;
            
            if(offsource) {
                I_nu = I_CMB[j];
            } else {
                I_nu = (1-t_src)*(1-u_src) * source->I_nu[start_x0y0 + j];
                I_nu += t_src*(1-u_src) * source->I_nu[start_x1y0 + j];
                I_nu += (1-t_src)*u_src * source->I_nu[start_x0y1 + j];
                I_nu += t_src*u_src * source->I_nu[start_x1y1 + j];
            }
        
            if(output->flag[i] == 0 or output->flag[i] == -1) {
                eta_ap = telescope->eta_ap_ON[j];
            } else {
                eta_ap = telescope->eta_ap_OFF[j];
            }
            
            PSD_nu[j] = eta_ap * eta_atm_interp * effs->eta_tot_chain * I_nu
                + ( effs->eta_tot_chain * (1 - eta_atm_interp) * I_atm[j] 
                + effs->eta_tot_gnd * I_gnd[j] 
                + effs->eta_tot_mir * I_tel[j]) 
                * CL*CL / (freq*freq);
        }
        calcPhotonNoise(instrument, PSD_nu, geno, output, i, nTimes);
    }
    delete[] PSD_nu;
}

void calcPhotonNoise(Instrument<double> *instrument, double *PSD_nu, 
        std::mt19937 &geno, Output<double> *output, int idx, int nTimes) {
    
    double freq;
    double sigma_k; // Noise per channel.
    double eta_kj; // Filter efficiency for bin j, at channel k
    double sqrt_samp = sqrt(0.5 * instrument->f_sample); // Constant term needed for noise calculation

    for(int k=0; k<instrument->nf_ch; k++) {
        double P_k = 0; // Initialise each channel to zero, for each timestep
        double NEP_accum = 0;

        // Can loop over bins again, cheap operations this time
        for(int j=0; j<instrument->f_spec.num; j++)
        {   
            freq = instrument->f_spec.start + instrument->f_spec.step * j;
            eta_kj = instrument->filterbank[k*instrument->f_spec.num + j];
            
            NEP_accum += PSD_nu[j] * eta_kj * (HP * freq + PSD_nu[j] * eta_kj + 2 * instrument->delta / instrument->eta_pb);
            P_k += PSD_nu[j] * eta_kj;
        }

        sigma_k = sqrt(2 * NEP_accum * instrument->f_spec.step) * sqrt_samp;
        P_k *= instrument->f_spec.step;

        std::normal_distribution<double> gg{0., sigma_k};
        P_k += gg(geno);
       
        // STORAGE: Add signal to signal array in output
        output->signal[k * nTimes + idx] = P_k; 
    }
}

void parallelJobsW2K(Instrument<double> *instrument, Atmosphere<double> *atmosphere, 
        CalOutput<double> *output, double *eta_atm, ArrSpec<double> PWV_atm, ArrSpec<double> f_atm,
        Effs<double> *effs, int nPWV, int start, int stop, double dPWV,
        double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx) {
    
    // Get starting time and chop parameters
    double freq; // Bin frequency
    double eta_kj; // Filter efficiency for bin j, at channel k
    double _PWV;

    double* eta_atm_interp = new double[f_atm.num];
    double* PSD_nu = new double[f_atm.num];
    

    for(int i=start; i<stop; i++) {
        _PWV = PWV_atm.start + i * dPWV;
        for(int j=0; j<instrument->f_spec.num; j++) { 
            freq = instrument->f_spec.start + instrument->f_spec.step * j;

            eta_atm_interp[j] = interpValue(_PWV, freq, 
                    PWV_atm, f_atm, eta_atm);
            
            PSD_nu[j] = ( effs->eta_tot_chain * (1 - eta_atm_interp[j]) * I_atm[j] 
                + effs->eta_tot_gnd * I_gnd[j] 
                + effs->eta_tot_mir * I_tel[j]) 
                * CL*CL / (freq*freq);
        }
        
        // In this loop, calculate P_k, NEP_k and noise
        for(int k=0; k<instrument->nf_ch; k++) {
            double P_k = 0; // Initialise each channel to zero, for each timestep
            double eta_atm_avg = 0;
            double eta_kj_accum = 0;

            // Can loop over bins again, cheap operations this time
            for(int j=0; j<instrument->f_spec.num; j++) { 
                freq = instrument->f_spec.start + instrument->f_spec.step * j;
                eta_kj = instrument->filterbank[k * instrument->f_spec.num + j];
                
                eta_atm_avg += eta_atm_interp[j] * eta_kj;
                eta_kj_accum += eta_kj;
                P_k += PSD_nu[j] * eta_kj;
            }

            // STORAGE: Add signal to signal array in output
            //printf("%d\n", nPWV*k + i);
            output->power[k * nPWV + i] = P_k * instrument->f_spec.step; 
            output->temperature[k * nPWV + i] = atmosphere->Tatm * (1 - eta_atm_avg/eta_kj_accum); 
            //printf("%d\n", nPWV*k + i);
        }
    }
    delete[] eta_atm_interp;
    delete[] PSD_nu;
}
