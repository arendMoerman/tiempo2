/*! \file InterpUtils.cpp
 * \brief Utilities for linear interpolation.
 **/

#include "InterpUtils.h"

double interpValue(double x, double y, 
            double x0, double y0, int size_x, int size_y, double dx, double dy, 
            double *vals, int offset, bool debug) {
    
    int idx_x = floorf((x - x0) / dx);
    int idx_y = floorf((y - y0) / dy);

    double f00 = vals[idx_x * size_y + idx_y + offset];
    double f10 = vals[(idx_x + 1) * size_y + idx_y + offset];
    double f01 = vals[idx_x * size_y + idx_y + 1 + offset];
    double f11 = vals[(idx_x + 1) * size_y + idx_y + 1 + offset];
    
    double t = (x - (x0 + dx*idx_x)) / dx;
    double u = (y - (y0 + dy*idx_y)) / dy;
    
    double fxy = (1-t)*(1-u)*f00 + t*(1-u)*f10 + t*u*f11 + (1-t)*u*f01;

    return fxy;
}
