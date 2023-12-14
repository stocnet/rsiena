/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: WeightedMixedTwoPathFunction.h
 *
 * Description: This file contains the definition of the
 * WeightedMixedTwoPathFunction class.
 *****************************************************************************/

#ifndef WEIGHTEDMIXEDTWOPATHFUNCTION_H_
#define WEIGHTEDMIXEDTWOPATHFUNCTION_H_

#include "DyadicCovariateMixedNetworkAlterFunction.h"

namespace siena
{

class WeightedMixedTwoPathFunction: public DyadicCovariateMixedNetworkAlterFunction
{
public:
	WeightedMixedTwoPathFunction(std::string firstNetworkName,
		std::string secondNetworkName, std::string dyadicCovariateName, bool excludeMissing);

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
};

}

#endif /* WEIGHTEDMIXEDTWOPATHFUNCTION_H_ */
