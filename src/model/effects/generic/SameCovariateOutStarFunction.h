/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateOutStarFunction.h
 *
 * Description: This file contains the definition of the
 * SameCovariateOutStarFunction class.
 *****************************************************************************/

#ifndef SAMECOVARIATEOUTSTARFUNCTION_H_
#define SAMECOVARIATEOUTSTARFUNCTION_H_

#include "CovariateNetworkAlterFunction.h"

namespace siena
{

class SameCovariateOutStarFunction: public CovariateNetworkAlterFunction
{
public:
	SameCovariateOutStarFunction(std::string networkName,
		std::string covariateName, bool excludeMissing);
		
	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
};

}

#endif /* SAMECOVARIATEOUTSTARFUNCTION_H_ */
