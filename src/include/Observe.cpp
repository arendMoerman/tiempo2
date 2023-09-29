#include "Observe.h"

double getPlanck(double T, double nu)
{
    
    double prefac = 2 * M_HP * nu*nu*nu / (M_CL*M_CL);
    double dist = 1 / (exp(M_HP*nu / (M_KB*T)) - 1);

    return prefac * dist;
}
