/*! \file InterfaceCPU.cpp
    \brief Implementations of library functions for simulations on CPU.

*/

#include "InterfaceCPU.h"

TIEMPO2_DLL void runTiEMPO2(Instrument instrument, Telescope telescope, Atmosphere atmosphere, Source source, SimParams simparams) {
    std::vector<std::thread> threadPool;
    threadPool.resize(simparams.nThreads);
    
    // Timesteps
    double dt = 1. / instrument.freq_sample;
    
    // Number of time evaluations. Time is nTimes * dt
    double nTimes = ceil(simparams.t_obs / dt);
    
    // Number of steps per thread
    int step = ceil(nTimes / simparams.nThreads);
    
    int final_step;

    for(int n=0; n < simparams.nThreads; n++) {
        if(n == (simparams.nThreads - 1)) {
            final_step = nTimes;
        }

        else {
            final_step = (n+1) * step;
        }

        threadPool[n] = std::thread(parallelJobs, instrument, telescope, atmosphere, source, n * step, final_step);
    }

    // Wait with execution until all threads are done
    for (std::thread &t : threadPool) {
        if (t.joinable()) {
            t.join();
        }
    }
}

TIEMPO2_DLL void parallelJobs(Instrument instrument, Telescope telescope, Atmosphere atmosphere, Source source, int start, int stop) {
    printf("under construction\n"); 
}


