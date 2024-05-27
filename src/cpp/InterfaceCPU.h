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
#include "FileIO.h"

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
    TIEMPO2_DLL void runTiEMPO2(Instrument<double> *instrument, Telescope<double> *telescope, 
                                Atmosphere<double> *atmosphere, Source<double> *source, 
                                Output<double> *output, int nTimes, int nThreads);

    TIEMPO2_DLL void calcW2K(Instrument<double> *instrument, Telescope<double> *telescope, 
                             Atmosphere<double> *atmosphere, CalOutput<double> *output,
                             int nPWV, int nThreads);
       
    TIEMPO2_DLL void getSourceSignal(Instrument<double> *instrument, Telescope<double> *telescope, 
                                     double *output, double *I_nu, 
                                     double PWV, bool ON);
    
    TIEMPO2_DLL void getEtaAtm(ArrSpec<double> f_src, double *output, double PWV);

    TIEMPO2_DLL void getNEP(Instrument<double> *instrument, Telescope<double> *telescope, 
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
void parallelJobs_1(Instrument<double> *instrument, Telescope<double> *telescope, 
                  Atmosphere<double> *atmosphere, Source<double> *source, 
                  Output<double> *output, double *eta_atm,
                  ArrSpec<double> PWV_atm, ArrSpec<double> f_atm,
                  Effs<double> *effs, int nTimes,
                  int start, int stop, double dt, 
                  double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx);

void parallelJobs_2(Instrument<double> *instrument, Telescope<double> *telescope, 
                  Atmosphere<double> *atmosphere, Source<double> *source, 
                  Output<double> *output, double *eta_atm,
                  ArrSpec<double> PWV_atm, ArrSpec<double> f_atm,
                  Effs<double> *effs, int nTimes,
                  int start, int stop, double dt, 
                  double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx);

void parallelJobs_3(Instrument<double> *instrument, Telescope<double> *telescope, 
                  Atmosphere<double> *atmosphere, Source<double> *source, 
                  Output<double> *output, double *eta_atm,
                  ArrSpec<double> PWV_atm, ArrSpec<double> f_atm,
                  Effs<double> *effs, int nTimes,
                  int start, int stop, double dt, 
                  double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx);

void parallelJobs(Instrument<double> *instrument, Telescope<double> *telescope, 
                  Atmosphere<double> *atmosphere, Source<double> *source, 
                  Output<double> *output, double *eta_atm,
                  ArrSpec<double> PWV_atm, ArrSpec<double> f_atm,
                  Effs<double> *effs, int nTimes,
                  int start, int stop, double dt,
                  double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx);

void calcPhotonNoise(Instrument<double> *instrument, double *PSD_nu, 
                     std::mt19937 &geno, Output<double> *output, int idx, int nTimes);

void parallelJobsW2K(Instrument<double> *instrument, Atmosphere<double> *atmosphere, 
                     CalOutput<double> *output, double *eta_atm, 
                     ArrSpec<double> PWV_atm, ArrSpec<double> f_atm, Effs<double> *effs, 
                     int nPWV, int start, int stop, double dPWV,
                     double* I_atm, double* I_gnd, double* I_tel, double *I_CMB, int threadIdx);

#endif
