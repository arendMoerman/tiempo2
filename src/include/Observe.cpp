#include "Observe.h"

double getPlanck(double T, double nu)
{
    
    double prefac = 2 * HP * nu*nu*nu / (CL*CL);
    double dist = 1 / (exp(HP*nu / (KB*T)) - 1);
    return prefac * dist;
}

void procData(Telescope* telescope, Instrument* instrument, Source* source, 
        double* ret, int j, double PSD_nu, bool chop_flag) {
    if (telescope->chop_mode == 0) {
        for(int k=0; k<instrument->nfreqs_filt; k++) {
            ret[k] += (instrument->filterbank[k*source->nfreqs_src + j] * PSD_nu);
        }
    }
    
    else if (telescope->chop_mode == 1) {
        for(int k=0; k<instrument->nfreqs_filt; k++){
            if (chop_flag){
                ret[k] -= (instrument->filterbank[k*source->nfreqs_src + j] * PSD_nu);
            }

            else {
                ret[k] += (instrument->filterbank[k*source->nfreqs_src + j] * PSD_nu);
            }
        }
    }
    
    else if (telescope->chop_mode == 2) {
        for(int k=0; k<instrument->nfreqs_filt; k++){
            if (chop_flag){
                ret[k] -= (instrument->filterbank[k*source->nfreqs_src + j] * PSD_nu);
            }

            else {
                ret[k] += (instrument->filterbank[k*source->nfreqs_src + j] * PSD_nu);
            }
        }
    }
}
