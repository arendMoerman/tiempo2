/*! \file InterpUtils.cpp
 * \brief Utilities for linear interpolation.
 **/

#include "InterpUtils.h"

int findIndexLow(double a, double *arr, int size_arr, bool debug) {
    double da = arr[1] - arr[0]; // Lowest value in arr. Note arr HAS TO BE monotonically increasing.
    if (da == 0) {
        printf("ERROR: input array not monotonically increasing.\n");
        exit(1);
    }

    int idx_start = floor((a - arr[0]) / da);

    return idx_start;
}

double interpValue(double x0, double y0, double *x, double *y, int size_x, int size_y, double *vals, int offset, bool debug) {
    int idx_x = findIndexLow(x0, x, size_x, debug); // Has to be size_x
    int idx_y = findIndexLow(y0, y, size_y, debug); // Has to be size_y
    double f00 = vals[idx_x * size_y + idx_y + offset];
    double f10 = vals[(idx_x + 1) * size_y + idx_y + offset];
    double f01 = vals[idx_x * size_y + idx_y + 1 + offset];
    double f11 = vals[(idx_x + 1) * size_y + idx_y + 1 + offset];
    
    //printf("x0 = %f, y0 = %f, f00 = %f, f10 = %f, f01 = %f, f11 = %f\n", x[idx_x], y[idx_y], f00, f10, f01, f11);
    
    double t = (x0 - x[idx_x]) / (x[idx_x + 1] - x[idx_x]);
    double u = (y0 - y[idx_y]) / (y[idx_y + 1] - y[idx_y]);
    
    double fxy = (1-t)*(1-u)*f00 + t*(1-u)*f10 + t*u*f11 + (1-t)*u*f01;

    return fxy;
}
