#include "Observe.h"

double getPlanck(double T, double nu)
{
    
    double prefac = 2 * HP * nu*nu*nu / (CL*CL);
    double dist = 1 / (exp(HP*nu / (KB*T)) - 1);
    return prefac * dist;
}
