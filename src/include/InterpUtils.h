/*! \file
 * \brief Utilities for linear interpolation.
 **/

#include <iostream>
#include <cmath>

#ifndef __INTERPUTILS_H
#define __INTERPUTILS_H

/*
 * Find index of lower bound of interval containing some specified number.
 *
 * @param a Value to find interval for.
 * @param arr Array of doubles, specifying interval start and ends. Note: should be complete interval, not a chunk.
 * @param size_vals Size of array containing values.
 *
 * @returns idx Index of value in `arr', corresponding to lower bound of interval containing `a'.
 **/
int findIndexLow(double a, double *arr, int size_vals);

/*
 * Linearly interpolate bivariate function.
 *
 * @param x0 Point in x to interpolate on.
 * @param y0 Point in y to interpolate on.
 * @param x Array of points representing x-coordinates.
 * @param y Array of points representing y-coordinates.
 * @param size_y Size of array containing y-coordinates.
 * @param vals Function values on grid spanning x and y.
 * @param size_vals Size of array containing values.
 * 
 * @returns val_interp Interpolated value of function on x0 and y0.
 **/
double interpValue(double x0, double y0, double *x, double *y, int size_x, int size_y, double *vals);

#endif
