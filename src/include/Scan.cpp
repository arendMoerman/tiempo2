/*! \file Scan.cpp
 * \brief Definition of scanning patterns and their implementations.
 **/

#include "Scan.h"

void scanPoint(AzEl* center, AzEl* out, bool chop, double sep) {
    double offset = 0.;
    
    if (chop) {
        offset = sep;
    }    
    
    out->Az = center->Az + offset;
    out->El = center->El;
}

void convertAnglesToSpatialAtm(AzEl* angles, xy_atm* out, double h_column) {
    
    double coord = tan(M_PI * angles->Az / 180.) * h_column;
    
    out->xAz = coord;
    coord = tan(M_PI * angles->El / 180.) * h_column;
    out->yEl = coord;
}
