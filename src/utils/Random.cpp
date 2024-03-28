/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Random.cpp
 *
 * Description: This file contains the implementation of random number
 * generation utilities.
 *****************************************************************************/

#include "Random.h"
#include <sys/time.h>
#include <Rmath.h>
#include <vector>
#undef length
#include <Rinternals.h>

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Methods for drawing random numbers
// ----------------------------------------------------------------------------

/**
 * Draws a uniformly distributed random double from the interval (0,1).
 */
double nextDouble()
{
    return unif_rand();
}


/**
 * Returns an exponential variate for the given rate parameter.
 * QAD for comparison with Delphi or for SienaProfile.exe
 */
double nextExponentialQAD(double lambda)
{
    return - log(nextDouble()) / lambda;
}


/**
 * Returns an exponential variate for the given rate parameter.
 * More sophisticated method from R for ongoing use.
 */
double nextExponential(double lambda)
{
#ifndef STANDALONE
	return rexp(1/lambda);
#endif
}

/**
 * Returns a gamma variate for the given shape and scale parameters.
 */

double nextGamma(double shape, double scale)
{
#ifndef STANDALONE
	return rgamma(shape, scale);
#endif
}

/**
 * Returns a normal variate for the given mean and standard deviation
 * parameters.
 */

double nextNormal(double mean, double standardDeviation)
{
#ifndef STANDALONE
	return rnorm(mean, standardDeviation);
#endif
}

/**
 * Returns the normal density for the given value, mean and standard deviation
 * parameters. The log will be returned if log is TRUE.
 */

double normalDensity(double value, double mean, double standardDeviation,
	 int log)
{
#ifndef STANDALONE
	return dnorm(value, mean, standardDeviation, log);
#endif
}
/**
 * Draws a uniformly distributed random integer from the interval [0,n).
 */
int nextInt(int n)
{
    return (int) (n * unif_rand());
}


/**
 * Draws a random integer from the interval [0,n) such that each value <i>i</i>
 * has the probability <i>p</i>[<i>i</i>] of occurrence. The probability vector
 * <i>p</i> should sum up to 1.
 */
int nextIntWithProbabilities(int n, const double * p)
{
    // Draw a uniform random number 'value' from (0,1) and iteratively
    // find the first i such that p[0] + p[1] + ... + p[i] >= value.

    double value = nextDouble();
    int i = 0;
    double sum = p[0];

    while (sum < value && i < n - 1)
    {
        i++;
        sum += p[i];
    }

    // Because of rounding errors, the sum of all probabilities might still
    // be lower than the random value. In this case return the last element
    // with a positive probability of occurrence.

    if (i == n - 1)
    {
        while (p[i] == 0 && i > 0)
        {
            i--;
        }
    }

    return i;
}


// Temporarily commented out as we don't have the Random class anymore.

///**
// * Stores the current state of the R random number seed
// * \return None
// */
//void Random::storeState()
//{
//    SEXP rs;
//    PROTECT(rs = install(".Random.seed"));
//    PutRNGstate();
//    SEXP ans;
//    PROTECT(ans = findVar(rs, R_GlobalEnv));
//    this->state = new std::vector<int>(Rf_length(ans));
//    /* copy ans to this->state */
//    int *start = INTEGER(ans);
//    for (int i = 0; i < Rf_length(ans); i++)
//    {
//        (*(this->state))[i] = *start++;
//    }
//    UNPROTECT(2);
//    return;
//}
//
//
///**
// * Restores the current state of the R random number seed from the value of the
// * state
// * \return None
// */
//void Random::restoreState()
//{
//    SEXP rs;
//    PROTECT(rs = install(".Random.seed"));
//    SEXP ans;
//    PROTECT(ans = allocVector(INTSXP, this->state->size()));
//    /* copy this->state to ans */
//    int *start = INTEGER(ans);
//    for (int i = 0; i < Rf_length(ans); i++)
//    {
//        *start++ = (*(this->state))[i];
//    }
//    /* set rs to be ans */
//    defineVar(rs, ans, R_GlobalEnv);
//    GetRNGstate();
//    UNPROTECT(2);
//    return;
//}

}
