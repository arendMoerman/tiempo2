#include <random>
#include <iostream>
#include <string>
#include <vector>

#define _USE_MATH_DEFINES

#ifndef __Random_h
#define __Random_h

/*! \file Random.h
    \brief Class for generating random numbers.

    Generate object for creating random numbers for rejection sampling. 
    Currently only for Gaussian beams, but want to include Poisson disk sampling too.
    If called empty, generate a random seed to use, which can be retrieved.
    If called with seed, use that seed.
*/

/**
 * Class for generating random numbers.
 *
 * Note that no function returns. All values are stored inside a variable which is passed by reference to the function.
 */
template <typename T> class Random
{
private:
    unsigned int seed;
    std::mt19937 gen; 
    
public:
    Random();
    Random(unsigned int seed);

    T generateUniform(T lower = -1.0);
};
#endif

/**
 * Initialize RNG. This constructor generates a random seed for the random draws.
 */
template <typename T>
Random<T>::Random()
{
    std::random_device rd;
    std::mt19937 geno(rd());
    this->seed = rd();
    this->gen = geno;
}

/**
 * Initialize RNG. This constructor takes a pre-set seed for the random draws.
 *
 * @param seed Positive integer determining the RNG seed.
 */
template <typename T>
Random<T>::Random(unsigned int seed)
{
    std::mt19937 geno(seed);
    this->seed = seed;
    this->gen = geno;
}

/**
 * Generate a random sample.
 *
 * @param lower Lower value of range. Defaults to -1.0.
 * 
 * @returns out Number between lower and 1.0..
 */
template <typename T>
T Random<T>::generateUniform(T lower)
{
    std::uniform_real_distribution<T> dis(lower, 1.0);
    return dis(this->gen);
}
