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
			int parameter, bool excludeMissing);

	virtual ~SameCovariateInTiesFunction();
	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual void preprocessEgo(int ego);
	virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
	bool lsameValue {};
	bool lsameVariable {};
	bool lroot {}; // should the square root be taken?
	bool laverage {}; // should the average be used?
	int * lpCovariateNumbers {};
	double lCovEgo {};
	int lCovNumberEgo {};
};

}

#endif /* SAMECOVARIATEINTIESFUNCTION_H_ */
