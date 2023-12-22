#include "InterpUtils.h"

int findIndexLow(double a, double *arr, int size_arr, bool debug) {
    double da = arr[1] - arr[0]; // Lowest value in arr. Note arr HAS TO BE monotonically increasing.
    if (da == 0) {
        printf("ERROR: input array not monotonically increasing.\n");
        exit(1);
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

double interpValue(double x0, double y0, double *x, double *y, int size_x, int size_y, double *vals, bool debug) {
    int idx_x = findIndexLow(x0, x, size_x, debug); // Has to be size_x
    int idx_y = findIndexLow(y0, y, size_y, debug); // Has to be size_y
    double f00 = vals[idx_x * size_y + idx_y];
    double f10 = vals[(idx_x + 1) * size_y + idx_y];
    double f01 = vals[idx_x * size_y + idx_y + 1];
    double f11 = vals[(idx_x + 1) * size_y + idx_y + 1];
    
    //printf("x0 = %f, y0 = %f, f00 = %f, f10 = %f, f01 = %f, f11 = %f\n", x[idx_x], y[idx_y], f00, f10, f01, f11);
    
    double t = (x0 - x[idx_x]) / (x[idx_x + 1] - x[idx_x]);
    double u = (y0 - y[idx_y]) / (y[idx_y + 1] - y[idx_y]);
    
    double fxy = (1-t)*(1-u)*f00 + t*(1-u)*f10 + t*u*f11 + (1-t)*u*f01;

    return fxy;

    //double dy = y[idx_y + 1] - y[idx_y];
    //double dy0 = y0 - y[idx_y];


    //double sx0 = (val10 - val00) / dx;
    //double sx1 = (val11 - val01) / dx;
    //
    //double sy = (val01 - val00 + (sx1 - sx0)*dx0) / dy;
    ////double sy = (sx1 - sx0)*dx0 / dy;

    //double val_interp = val00 + sx0*dx0 + sy*dy0;
    ////double val_interp = sx0*dx0 + sy*dy0;
    ////printf("%f, %f, %f \n", val_interp, val00, val11);

    //return val_interp;
}
