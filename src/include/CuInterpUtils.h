/*! \file CuInterpUtils.h
 * \brief Utilities for linear interpolation for CUDA.
 **/

#include "cuda.h"
#include "math.h"

#ifndef __CUINTERPUTILS_H
#define __CUINTERPUTILS_H

/**
  Find index of lower bound of interval containing some specified number.
 
  @param a Value to find interval for.
  @param arr Array of doubles, specifying interval starts and ends. Note: should be complete interval, not a chunk.
  @param size_vals Size of array containing values.
 
  @returns idx Index of value in `arr', corresponding to lower bound of interval containing `a'.
 */
__device__ int findIndexLow(float a, float *arr, int size_arr) {
    float da = arr[1] - arr[0]; // Lowest value in arr. Note arr HAS TO BE monotonically increasing.
    if (da == 0) {
        printf("ERROR: input array not monotonically increasing.\n");
    }

    //int idx_start = floor(a / da);
    int idx_start = floor((a - arr[0]) / (arr[size_arr-1] - arr[0]) * size_arr);

    // Reduce index if larger than size_arr - 1
    while ((size_arr-1) <= idx_start) {
        idx_start--;
    }
     
jump:
    bool is_lower = arr[idx_start] <= a;
    bool is_higher = arr[idx_start + 1] > a;

    //if(debug) {printf("%d\n", idx_start);}

    if (is_lower && is_higher) {return idx_start;}

    else if (is_lower && !is_higher) {idx_start++; goto jump;}
    else if (is_higher && !is_lower) {idx_start--; goto jump;}
}

/**
  Linearly interpolate bivariate function.
 
  @param x0 Point in x to interpolate on.
  @param y0 Point in y to interpolate on.
  @param x Array of points representing x-coordinates.
  @param y Array of points representing y-coordinates.
  @param size_y Size of array containing y-coordinates.
  @param vals Function values on grid spanning x and y.
  @param size_vals Size of array containing values.
  
  @returns val_interp Interpolated value of function on x0 and y0.
 */
__device__ float interpValue(float x0, float y0, float *x, float *y, int size_x, int size_y, float *vals, int offset) {
    int idx_x = findIndexLow(x0, x, size_x); // Has to be size_x
    int idx_y = findIndexLow(y0, y, size_y); // Has to be size_y
    float f00 = vals[idx_x * size_y + idx_y + offset];
    float f10 = vals[(idx_x + 1) * size_y + idx_y + offset];
    float f01 = vals[idx_x * size_y + idx_y + 1 + offset];
    float f11 = vals[(idx_x + 1) * size_y + idx_y + 1 + offset];
    
    float t = (x0 - x[idx_x]) / (x[idx_x + 1] - x[idx_x]);
    float u = (y0 - y[idx_y]) / (y[idx_y + 1] - y[idx_y]);
    
    float fxy = (1-t)*(1-u)*f00 + t*(1-u)*f10 + t*u*f11 + (1-t)*u*f01;

    return fxy;
}

#endif
