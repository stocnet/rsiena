/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleEqualCovariateFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * DoubleEqualCovariateFunction.
 *****************************************************************************/

#include "DoubleEqualCovariateFunction.h"
#include <cmath>

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] covariateName1 the name of the first covariate this function is
 * associated with
 * @param[in] covariateName2 the name of the second covariate this function is
 * associated with
 * @param[in] excludeMissing whether to exclude missing values
 */

DoubleEqualCovariateFunction::DoubleEqualCovariateFunction(std::string covariateName1,
						std::string covariateName2, bool excludeMissing) :
	DoubleCovariateFunction(covariateName1, covariateName2)
{
	this->lexcludeMissing = excludeMissing;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DoubleEqualCovariateFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	DoubleCovariateFunction::initialize(pData, pState, period, pCache);
}

/**
 * Returns if the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double DoubleEqualCovariateFunction::value(int alter) const
{
	double statistic = 0;
	if (!(this->lexcludeMissing && 
				(this->firstMissing(this->ego()) || this->secondMissing(alter))))
	{
		if ((this->firstCovariateIntValue(this->ego()) ==
					this->secondCovariateValue(alter)) && 
							(this->secondCovariateIntValue(alter) >= 1))
		{
			statistic = 1;
		}
	}
	return statistic;
}
}
