/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateNetworkAlterFunction.h
 *
 * Description: This file contains the definition of the
 * CovariateNetworkAlterFunction class.
 *****************************************************************************/


#ifndef COVARIATENETWORKALTERFUNCTION_H_
#define COVARIATENETWORKALTERFUNCTION_H_

#include <string>
#include "NetworkAlterFunction.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;
class NetworkCache;
class ConstantCovariate;
class ChangingCovariate;
class BehaviorVariable;
class BehaviorLongitudinalData;


class CovariateNetworkAlterFunction: public NetworkAlterFunction
{
public:
	CovariateNetworkAlterFunction(std::string networkName, std::string covariateName);
	virtual ~CovariateNetworkAlterFunction();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

protected:

	// this is the overall observed mean;
	// except for centered actor covariates (not dependent behavior variables),
	// for which it is 0.
	double covmean() const;

// value of the covariate:
	double covvalue(int alter) const;
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

#endif /* COVARIATENETWORKALTERFUNCTION_H_ */
