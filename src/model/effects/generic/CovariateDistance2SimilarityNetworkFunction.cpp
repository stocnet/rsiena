/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2SimilarityNetworkFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * CovariateDistance2SimilarityNetworkFunction.
 *****************************************************************************/

#include "CovariateDistance2SimilarityNetworkFunction.h"
#include "network/Network.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 * @param[in] covariateName the name of the covariate this function is
 * associated with
 * @param[in] parameter the value of the internal effect parameter this
 * function is associated with
 * @param[in] excludeMissing: whether to exclude missing values
 */
CovariateDistance2SimilarityNetworkFunction::
CovariateDistance2SimilarityNetworkFunction(string networkName, string
	covariateName, bool excludeMissing) :
	CovariateDistance2NetworkFunction(networkName, covariateName, excludeMissing, true)
{
	this->lexcludeMissing = excludeMissing;
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double CovariateDistance2SimilarityNetworkFunction::value(int alter)
{
	double value = 0;
	if (!this->lexcludeMissing || (!this->missingDummy(alter) &&
				!this->missingDummy(this->ego())))
	{
		value = this->similarityAvAlt(this->ego(), alter);
  	}
	return value;
}


}
