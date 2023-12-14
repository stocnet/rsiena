/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateInTiesFunction.h
 *
 * Description: This file contains the definition of the
 * SameCovariateInTiesFunction class.
 *****************************************************************************/

#ifndef SAMECOVARIATEINTIESFUNCTION_H_
#define SAMECOVARIATEINTIESFUNCTION_H_

#include "CovariateNetworkAlterFunction.h"
#include "model/effects/NetworkEffect.h"

namespace siena
{

class SameCovariateInTiesFunction: public CovariateNetworkAlterFunction
{
public:
	SameCovariateInTiesFunction(std::string networkName,
			std::string covariateName, bool sameValue, bool sameVariable,
			bool excludeMissing);

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual void preprocessEgo(int ego);
	virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
	bool lsameValue {};
	bool lsameVariable {};
};

}

#endif /* SAMECOVARIATEINTIESFUNCTION_H_ */
