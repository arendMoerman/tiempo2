/*! \file InterfaceCPU.h
    \brief Declarations of library functions for simulations on CPU.

*/
#include <thread>
#include <vector>

#include "InterpUtils.h"
#include "Structs.h"
#include "Scan.h"

#ifdef _WIN32
#   define TIEMPO2_DLL __declspec(dllexport)
#else
#   define TIEMPO2_DLL
#endif

#ifndef __InterfaceCPU_h
#define __InterfaceCPU_h

extern "C"
{
    TIEMPO2_DLL void runTiEMPO2(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, SimParams *simparams);

    TIEMPO2_DLL void parallelJobs(Instrument *instrument, Telescope *telescope, Atmosphere *atmosphere, Source *source, int start, int stop);


}
#endif
