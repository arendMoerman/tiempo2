/*! \file InterfaceCPU.h
    \brief Declarations of library functions for simulations on CPU.
*/

#include <thread>
#include <vector>
#include <random>

#include "InterpUtils.h"
#include "Structs.h"
#include "Scan.h"
#include "Timer.h"

#ifdef _WIN32
#   define TIEMPO2_DLL __declspec(dllexport)
#else
#   define TIEMPO2_DLL
#endif

#ifndef __InterfaceCPU_h
#define __InterfaceCPU_h

#define SI_TO_MJY               1E20 /* SI to MJy*/

#define PI 3.14159265358979323846  /* pi */
#define CL 2.9979246E8 // m s^-1
#define HP 6.62607015e-34
#define KB 1.380649e-23

extern "C"
{
    TIEMPO2_DLL void runTiEMPO2(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, 
        SimParams *simparams, Output *output);

    TIEMPO2_DLL void parallelJobs(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, SimParams* simparams, Output* output, 
        Effs* effs, int start, int stop, double dt, int num_AzEl, 
        double* I_atm, double* I_gnd, double* I_tel, int threadIdx);
       
    TIEMPO2_DLL void getSourceSignal(Instrument *instrument, Telescope *telescope, Source *source, double *output, double Az, double El);

    TIEMPO2_DLL void getNEP(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, double *output, double PWV);
}
#endif
