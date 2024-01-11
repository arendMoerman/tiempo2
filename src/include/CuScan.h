/*! \file CuScan.h
 * \brief File containing scanning patterns and their implementations for CUDA.
 **/

#include <math.h>
#include <cuda.h>

#ifndef __CUSCAN_h
#define __CUSCAN_h

/**
  Structure for storing an Azimuth-Elevation co-ordinate.
 */
struct AzEl {
    float Az;      /**< Azimuth angle on-sky, in degrees.*/
    float El;      /**< Elevation angle in degrees.*/
};                                                                                             
                                                                                               
                                                                                            
struct xy_atm {
    float xAz;     /** < X-coordinate on atmosphere at reference height corresponding to Az.*/
    float yEl;     /** < Y-coordinate on atmosphere at reference height corresponding to El.*/
};

/**
  Calculate new Azimuth-Elevation co-ordinate, accoding to chop position.

  This function just "scans" a single point, so seems sort of pointless. 
  Still implemented for completeness.

  @param center Az-El co-ordinate of point to observe, w.r.t. source Az-El.
  @param out Container for storing output Az-El co-ordinate.
  @param chop Whether chopper is in A (false) or B (true).
  @param sep Angular throw between chop A and B, in degrees.
 */
__device__ __inline__ void scanPoint(AzEl* center, AzEl* out, bool chop, float sep = 0.) {
    float offset = 0.;
    
    if (chop) {
        offset = sep;
    }    
    
    out->Az = center->Az + offset;
    out->El = center->El;
}

/**
  Convert an Az-El co-ordinate to a projected x-y co-ordinate on the atmosphere.

  @param angles Az-El co-ordinate to convert.
  @param out Container for storing the calculated x-y point.
  @param h_column Reference column height of atmosphere.
 */
__device__ __inline__ void convertAnglesToSpatialAtm(AzEl* angles, xy_atm* out, float h_column) {
    
    float coord = tanf(M_PI * angles->Az / 180.) * h_column;
    
    out->xAz = coord;
    coord = tanf(M_PI * angles->El / 180.) * h_column;
    out->yEl = coord;
}

#endif
