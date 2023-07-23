#include "Scan.h"

AzEl scanPoint(AzEl center, bool chop, double sep) {
    AzEl out;
    
    double offset = 0.;
    
    if (chop) {
        offset = sep;
    }    
    
    out.Az = center.Az + offset;
    out.El = center.El;

    return out;
}

xy_atm convertAnglesToSpatialAtm(AzEl angles, double h_column) {
    xy_atm out;

    double coord = tan(M_PI * angles.Az / 180.) * h_column;
    out.xAz = coord;
    coord = tan(M_PI * angles.El / 180.) * h_column;
    out.yEl = coord;

    return out;
}
