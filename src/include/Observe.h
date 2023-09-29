/*! \file
 * \brief File containing functions used during simulations of observations.
 **/

#include <cmath>
#include <cstdio>
#include "InterpUtils.h"

#ifndef __OBSERVE_h
#define __OBSERVE_h

#define M_PI           3.14159265358979323846  /* pi */
#define M_CL           299792458 /* Speed of light */
#define M_HP           6.62607015E-34
#define M_KB           1.380649E-23 

/**
 * Get specific intensity of an object at temperature T and frequency nu.
 */
double getPlanck(double T, double nu);

#endif
