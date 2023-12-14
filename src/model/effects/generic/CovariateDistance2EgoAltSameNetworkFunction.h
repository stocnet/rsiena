/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2EgoAltSameNetworkFunction.h
 *
 * Description: This file contains the definition of the
 * CovariateDistance2EgoAltSameNetworkFunction class.
 *****************************************************************************/

#ifndef COVARIATEDISTANCE2EGOALTSAMENETWORKFUNCTION_H_
#define COVARIATEDISTANCE2EGOALTSAMENETWORKFUNCTION_H_

#include "CovariateDistance2NetworkFunction.h"

namespace siena
{

class CovariateDistance2EgoAltSameNetworkFunction: public
	CovariateDistance2NetworkFunction
{
public:
	CovariateDistance2EgoAltSameNetworkFunction(std::string networkName,
		std::string covariateName, bool excludeMissing, bool outgoing, 
									double parameter);
	virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
	bool loutgoing {};
	bool ltrunc {};
};

}

#endif /* COVARIATEDISTANCE2EGOALTSAMENETWORKFUNCTION_H_ */
