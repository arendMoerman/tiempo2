#include "InterpUtils.h"

int findIndexLow(double a, double *arr, int size_arr) {
    double da = arr[1] - arr[0]; // Lowest value in arr. Note arr HAS TO BE monotonically increasing.

    if (da == 0) {
        printf("ERROR: input array not monotonically increasing.\n");
        exit(1);
    }

    int idx_start = floor(a / da);
    if (size_arr == idx_start) {
        printf("WARNING: value a = %f larger than maximum in arr. Reducing a.\n", a);
        a = a - da;
        idx_start--;
    }
    
jump:
    bool is_lower = arr[idx_start] <= a;
    bool is_higher = arr[idx_start + 1] > a;

    if (is_lower && is_higher) {return idx_start;}

    else if (is_lower) {idx_start++; goto jump;}
    else if (is_higher) {idx_start--; goto jump;}
}

double interpValue(double x0, double y0, double *x, double *y, int size_y, double *vals, int size_vals) {
    int idx_x = findIndexLow(x0, x, size_vals);
    int idx_y = findIndexLow(y0, y, size_vals);
    
    double val00 = vals[idx_x * size_y + idx_y];
    double val10 = vals[(idx_x + 1) * size_y + idx_y];
    double val01 = vals[idx_x * size_y + idx_y + 1];
    double val11 = vals[(idx_x + 1) * size_y + idx_y + 1];
    double dx = x[idx_x + 1] - x[idx_x];
    double dx0 = x0 - x[idx_x];
    
    double dy = y[idx_y + 1] - y[idx_y];
    double dy0 = y0 - y[idx_y];

    double sx0 = (val10 - val00) / dx;
    double sx1 = (val11 - val01) / dx;
    
    double sy = (val01 - val00 + (sx1 - sx0)*dx0) / dy;

    double val_interp = val00 + sx0*dx0 + sy*dy0;
    printf("%f, %f, %f \n", val_interp, val00, val11);

    return val_interp;
}
