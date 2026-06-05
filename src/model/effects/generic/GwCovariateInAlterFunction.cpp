/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwCovariateInAlterFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * GwCovariateInAlterFunction.
 *****************************************************************************/

#include <cmath>
#include <stdexcept>

#include "GwCovariateInAlterFunction.h"
#include "network/Network.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable
 * @param[in] covariateName the name of the behavior variable used as covariate
 * @param[in] excludeMissing whether to return 0 when all in-alters are missing
 * @param[in] parameter the GW decay parameter (must be >= 0)
 */
GwCovariateInAlterFunction::GwCovariateInAlterFunction(string networkName,
	string covariateName, bool excludeMissing, double parameter) :
	CovariateDistance2NetworkFunction(networkName, covariateName,
		excludeMissing, false, true) // loutgoing=false, lraw=true
{
	this->lexcludeMissing = excludeMissing;
	this->lweight = -0.01 * parameter;
	this->lexpmweight = exp(-this->lweight);
	this->lexpfactor = 1 - exp(this->lweight);

	if (parameter < 0)
	{
		throw runtime_error(
			"GwCovariateInAlterFunction must have nonnegative effect parameter");
	}
}


/**
 * Returns the GW-weighted raw in-alter behavior total for the given alter,
 * excluding ego's own contribution.
 */
double GwCovariateInAlterFunction::value(int alter) const
{
	if (this->lexcludeMissing && this->missingInDummy(alter))
	{
		return 0;
	}

	double total = this->totalInAlterValue(alter);

	// If ego -> alter tie exists, the cache includes ego's raw behavior;
	// subtract it to get T_alter^{(-ego)}.
	if (this->pNetwork()->tieValue(this->ego(), alter) == 1)
	{
		total -= CovariateNetworkAlterFunction::rawCovvalue(this->ego());
	}

	return this->gwWeight(total);
}


/**
 * Returns e^alpha * (1 - (1 - e^{-alpha})^total), or 0 if total <= 0.
 */
double GwCovariateInAlterFunction::gwWeight(double total) const
{
	if (total <= 0)
	{
		return 0;
	}
	return this->lexpmweight * (1 - pow(this->lexpfactor, total));
}

}
