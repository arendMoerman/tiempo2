/*! \file CuInterpUtils.h
 * \brief Utilities for linear interpolation for CUDA.
 **/

#include "cuda.h"
#include "math.h"

#ifndef __CUINTERPUTILS_H
#define __CUINTERPUTILS_H

/**
  Linearly interpolate bivariate function.
 
  @param x Point in x to interpolate on.
  @param y Point in y to interpolate on.
  @param x0 Start of x-coordinates.
  @param y0 Start of y-coordinates.
  @param size_x Size of array containing x-coordinates.
  @param size_y Size of array containing y-coordinates.
  @param dx Stepsize of x.
  @param dy Stepsize of y.
  @param vals Function values on grid spanning x and y.
  @param size_vals Size of array containing values.
  
  @returns val_interp Interpolated value of function on x0 and y0.
 */
__device__ float interpValue(float x, float y, 
            float x0, float y0, int size_x, int size_y, float dx, float dy,
            float *vals, int offset) {
    
    int idx_x = floorf((x - x0) / dx);
    int idx_y = floorf((y - y0) / dy);

    float f00 = vals[idx_x * size_y + idx_y + offset];
    float f10 = vals[(idx_x + 1) * size_y + idx_y + offset];
    float f01 = vals[idx_x * size_y + idx_y + 1 + offset];
    float f11 = vals[(idx_x + 1) * size_y + idx_y + 1 + offset];
    
    float t = (x - (x0 + dx*idx_x)) / dx;
    float u = (y - (y0 + dy*idx_y)) / dy;
    
    float fxy = (1-t)*(1-u)*f00 + t*(1-u)*f10 + t*u*f11 + (1-t)*u*f01;

    return fxy;
}

#endif
