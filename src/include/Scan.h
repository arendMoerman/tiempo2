/*! \file
 * \brief File containing scanning patterns and their implementations.
 **/

#include <cmath>
#include <cstdio>

#ifndef __SCAN_h
#define __SCAN_h

# define M_PI           3.14159265358979323846  /* pi */

struct AzEl {
    double Az;      /**< Azimuth angle on-sky, in degrees.*/
    double El;      /**< Elevation angle in degrees.*/
};

struct xy_atm {
    double xAz;     /** < X-coordinate on atmosphere at reference height corresponding to Az.*/
    double yEl;     /** < Y-coordinate on atmosphere at reference height corresponding to El.*/
};

/**
 * Obtain azimuth and elevation angle corresponding to point on sky.
 *
 * @param center Struct containing center azimuth and elevation angles. If chopping, corresponds to path containing source.
 * @param chop Whether to chop or not.
 * @param sep If chopping, angular separation between the two paths.
 *
 * @returns out Struct containing azimuth and elevation angles of point on sky.
 **/
AzEl scanPoint(AzEl center, bool chop, double sep = 0.);

/**
 * Convert an azimuth and elevation angle on-sky to a position on atmospheric screen.
 *
 * @param angles Struct containing azimuth and elevation angles.
 * @param h_column Reference height of atmospheric column.
 *
 * @returns out Struct containing x and y-coordinates on atmosphere corresponding to supplied angles.
 **/
xy_atm convertAnglesToSpatialAtm(AzEl angles, double h_column);

#endif
