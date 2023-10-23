/*! \file
 * \brief File containing functions used during simulations of observations.
 **/

#include <cmath>
#include <cstdio>
#include "InterpUtils.h"
#include "Structs.h"

#ifndef __OBSERVE_h
#define __OBSERVE_h

const double PI = 3.14159265358979323846;  /* pi */
const double CL = 299792458; /* Speed of light */
const double HP = 6.62607015e-34;
const double KB = 1.380649e-23;

/**
 * Get specific intensity of an object at temperature T and frequency nu.
 */
double getPlanck(double T, double nu);

#endif
