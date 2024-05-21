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
#define HP 6.62607015E-34
#define KB 1.380649E-23

extern "C"
{
    TIEMPO2_DLL void runTiEMPO2(Instrument *instrument, Telescope *telescope, 
                                Atmosphere *atmosphere, Source *source, 
                                Output *output, int nTimes, int nThreads);

    TIEMPO2_DLL void calcW2K(Instrument *instrument, Telescope *telescope, 
                             Atmosphere *atmosphere, CalOutput *output,
                             int nPWV, int nThreads);
       
    TIEMPO2_DLL void getSourceSignal(Instrument *instrument, Telescope *telescope, 
                                     double *output, double *I_nu, double *eta_atm, 
                                     ArrSpec f_atm, ArrSpec PWV_atm, 
                                     double PWV, bool ON);
    
    TIEMPO2_DLL void getEtaAtm(ArrSpec f_src, double *output, double *eta_atm, 
                               ArrSpec f_atm, ArrSpec PWV_atm, double PWV);

    TIEMPO2_DLL void getNEP(Instrument *instrument, Telescope *telescope, 
                            double *eta_atm, 
                            ArrSpec f_atm, ArrSpec PWV_atm, 
                            double *output, double PWV, double Tatm);

}

/**
 * Parallel job for no chop, no scan, observations.
 *
 * @param instrument Instrument structure.
 * @param telescope Telescope structure.
 * @param instrument Instrument structure.
 * @param telescope Telescope structure.
 */
void parallelJobs_1(Instrument *instrument, Telescope *telescope, 
                  Atmosphere *atmosphere, Source *source, 
                  Output* output, Effs* effs, int nTimes,
                  int start, int stop, double dt, 
                  double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx);

void parallelJobs_2(Instrument *instrument, Telescope *telescope, 
                  Atmosphere *atmosphere, Source *source, 
                  Output* output, Effs* effs, int nTimes,
                  int start, int stop, double dt, 
                  double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx);

void parallelJobs_3(Instrument *instrument, Telescope *telescope, 
                  Atmosphere *atmosphere, Source *source, 
                  Output* output, Effs* effs, int nTimes,
                  int start, int stop, double dt, 
                  double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx);

void parallelJobs(Instrument *instrument, Telescope *telescope, 
                  Atmosphere *atmosphere, Source *source, 
                  Output* output, Effs* effs, int nTimes,
                  int start, int stop, double dt,
                  double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx);

void calcPhotonNoise(Instrument *instrument, double *PSD_nu, 
                     std::mt19937 &geno, Output *output, int idx, int nTimes);

void parallelJobsW2K(Instrument *instrument, Atmosphere *atmosphere, 
                     CalOutput* output, Effs* effs, 
                     int nPWV, int start, int stop, double dPWV,
                     double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx);

#endif
