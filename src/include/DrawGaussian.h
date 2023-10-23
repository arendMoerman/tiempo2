#include <iostream>
#include <array>
#include <cmath>
#include <complex>

#include "Random.h"
#define M_PI 3.14159265358979323846

#ifndef __DrawGaussian_h
#define __DrawGaussian_h

/*! \file DrawGaussian.h
    \brief Draw random number from Gaussian.
    
    Used to calculate random noise on top of MKID power.
*/

/** 
 * Draw random sample from Gaussian using rejection sampling.
 *
 * @param sigma standard deviation of Gaussian, obtained from NEP.
 */
double drawGaussian(double sigma);

/**
 * Probability density function for drawing sample from Gaussian.
 *
 * @param vars Vector of length 4, containing the xy positions and tilts of the ray to be checked.
 * @param scales Vector of length 4 containing the scale factors along xy and tilts
 */
double pdfGauss(double var, double sigma);

/** 
 * Initialize Gaussian beam from GPODict or GPODictf.
 *
 * Takes a GPODict or GPODictf and generates two c2Bundle or c2Bundlef objects, which contain the field and 
 *      associated currents and are allocated to passed pointer arguments.
 *
 * @param gdict GPODict or GPODictf object from which to generate a Gaussian beam.
 * @param refldict reflparams or reflparamsf object corresponding to surface on
 *      which to generate the Gaussian beam.
 * @param res_field Pointer to c2Bundle or c2Bundlef object.
 * @param res_current Pointer to c2Bundle or c2Bundlef object.
 *
 * @see GPODict
 * @see GPODictf
 * @see reflparams
 * @see reflparamsf
 * @see c2Bundle
 * @see c2Bundlef
 */

double drawGaussian(double sigma)
{
    double thres = 3.; // Choose 3-sigma point

    Random<double> rando;

    // Start rejection sampling. Use n_suc as succes counter
    bool suc = false;
    double yi;
    double xi;
    
    double lower = 0.;

    while (!suc)
    {
       // First, generate y-value between 0 and 1.
       yi = rando.generateUniform(lower);
       
       // Then, generate x value between -3*sigma and 3 *sigma.
       xi = rando.generateUniform() * thres * sigma;

       if (pdfGauss(xi, sigma) > yi) 
       {
           suc = true;
       }
    }
    return xi;
}

double pdfGauss(double var, double sigma)
{
    double norm = 1 / (sqrt(2*M_PI) * sigma);

    return norm * exp(-0.5 * var*var / (sigma*sigma));
}

#endif
