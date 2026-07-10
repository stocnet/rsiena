/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwCovariateAlterFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * GwCovariateAlterFunction.
 *****************************************************************************/

#include <cmath>
#include <stdexcept>

#include "GwCovariateAlterFunction.h"
#include "network/Network.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable
 * @param[in] covariateName the name of the behavior variable used as covariate
 * @param[in] excludeMissing whether to return 0 when all out-alters are missing
 * @param[in] parameter the GW decay parameter (must be >= 0)
 */
GwCovariateAlterFunction::GwCovariateAlterFunction(string networkName,
	string covariateName, bool excludeMissing, double parameter) :
	CovariateDistance2NetworkFunction(networkName, covariateName,
		excludeMissing, true, true) // loutgoing=true, nc=true
{
	this->lexcludeMissing = excludeMissing;
	this->lweight = -0.01 * parameter;
	this->lexpmweight = exp(-this->lweight);
	this->lexpfactor = 1 - exp(this->lweight);

	if (parameter < 0)
	{
		throw runtime_error(
			"GwCovariateAlterFunction must have nonnegative effect parameter");
	}
}


/**
 * Returns the GW-weighted raw out-alter behavior total for the given alter,
 * excluding ego's own contribution.
 */
double GwCovariateAlterFunction::value(int alter) const
{
	if (this->lexcludeMissing && this->missingDummy(alter))
	{
		return 0;
	}

	double total = this->totalAlterValue(alter);

	// totalAlterValue(alter) sums over alter's OUT-neighbors; ego is one of
	// them only if a tie alter -> ego exists. Subtract ego's raw behavior in
	// that case to get T_alter^{(-ego)}.
	if (this->pNetwork()->tieValue(alter, this->ego()) == 1)
	{
		total -= CovariateNetworkAlterFunction::uncenteredCovvalue(this->ego());
	}

	return this->gwWeight(total);
}


/**
 * Returns e^alpha * (1 - (1 - e^{-alpha})^total), or 0 if total <= 0.
 */
double GwCovariateAlterFunction::gwWeight(double total) const
{
	if (total <= 0)
	{
		return 0;
	}
	return this->lexpmweight * (1 - pow(this->lexpfactor, total));
}

}
