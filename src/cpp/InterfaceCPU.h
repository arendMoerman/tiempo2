/*! \file InterfaceCPU.h
    \brief Declarations of library functions for simulations on CPU.
*/
#include <thread>
#include <vector>

#include "InterpUtils.h"
#include "Structs.h"
#include "Scan.h"
#include "Observe.h"
#include "DrawGaussian.h"

#ifdef _WIN32
#   define TIEMPO2_DLL __declspec(dllexport)
#else
#   define TIEMPO2_DLL
#endif

#ifndef __InterfaceCPU_h
#define __InterfaceCPU_h

#define SI_TO_MJY               1E20 /* SI to MJy*/

extern "C"
{
    TIEMPO2_DLL void runTiEMPO2(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, 
        SimParams *simparams, Output *output);

    TIEMPO2_DLL void parallelJobs(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, SimParams* simparams, Output* output, 
        Effs* effs, int start, int stop, double dt, int num_AzEl, 
        double* slice_container, double* I_atm, double* I_gnd, double* I_tel, int threadIdx);
}
#endif
