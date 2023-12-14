/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2EgoAltSimNetworkFunction.h
 *
 * Description: This file contains the definition of the
 * CovariateDistance2EgoAltSimNetworkFunction class.
 *****************************************************************************/

#ifndef COVARIATEDISTANCE2EGOALTSIMNETWORKFUNCTION_H_
#define COVARIATEDISTANCE2EGOALTSIMNETWORKFUNCTION_H_

#include "CovariateDistance2NetworkFunction.h"

namespace siena
{

class CovariateDistance2EgoAltSimNetworkFunction: public
	CovariateDistance2NetworkFunction
{
public:
	CovariateDistance2EgoAltSimNetworkFunction(std::string networkName,
		std::string covariateName, bool excludeMissing, bool outgoing);
	virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
	bool loutgoing {};
};

}

#endif /* COVARIATEDISTANCE2EGOALTSIMNETWORKFUNCTION_H_ */
