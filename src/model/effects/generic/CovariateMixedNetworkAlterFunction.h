/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateMixedNetworkAlterFunction.h
 *
 * Description: This file contains the definition of the
 * CovariateMixedNetworkAlterFunction class.
 *****************************************************************************/

// combines MixedNetworkAlterFunction and CovariateNetworkAlterFunction

#ifndef COVARIATEMIXEDNETWORKALTERFUNCTION_H_
#define COVARIATEMIXEDNETWORKALTERFUNCTION_H_

#include <string>
#include "AlterFunction.h"
#include "utils/NamedObject.h"
#include "MixedNetworkAlterFunction.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;
class TwoNetworkCache;
class NetworkCache;
class ConstantCovariate;
class ChangingCovariate;
class BehaviorVariable;
class BehaviorLongitudinalData;

class CovariateMixedNetworkAlterFunction: public MixedNetworkAlterFunction
{
public:
	CovariateMixedNetworkAlterFunction(std::string firstNetworkName,
		std::string secondNetworkName, std::string covariateName);
	virtual ~CovariateMixedNetworkAlterFunction();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

protected:
	double covvalue(int i) const;
	bool missing(int i) const;
	double actor_similarity(int i, int j) const;
	ConstantCovariate * pConstantCovariate() const;
	ChangingCovariate * pChangingCovariate() const;
	BehaviorLongitudinalData * pBehaviorData() const;

private:
	std::string lcovariateName {};
	int lperiod {};
	ConstantCovariate * lpConstantCovariate;
	ChangingCovariate * lpChangingCovariate;
	BehaviorLongitudinalData * lpBehaviorData;

	// The current value of a behavior variable per each actor.
	// This array is 0 for covariate-based effects.
	const int * lvalues {};
};


}

#endif /* COVARIATEMIXEDNETWORKALTERFUNCTION_H_ */
