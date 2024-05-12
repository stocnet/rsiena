/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateOutTiesFunction.h
 *
 * Description: This file contains the definition of the
 * SameCovariateOutTiesFunction class.
 *****************************************************************************/

#ifndef SAMECOVARIATEOUTTIESFUNCTION_H_
#define SAMECOVARIATEOUTTIESFUNCTION_H_

#include "CovariateNetworkAlterFunction.h"
#include "model/effects/NetworkEffect.h"

namespace siena
{

class SameCovariateOutTiesFunction: public CovariateNetworkAlterFunction
{
public:
	SameCovariateOutTiesFunction(std::string networkName,
		std::string covariateName, bool sameValue, bool excludeMissing);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);
	virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
	bool lsameValue {};
};

}

#endif /* SAMECOVARIATEOUTTIESFUNCTION_H_ */
