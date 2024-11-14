/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TotalGwdspAlterEffect.h
 *
 * Description: This file contains the declaration of the class GwdspEffect.
 *
 *****************************************************************************/

#ifndef TOTALGWDSPALTEREFFECT_H_
#define TOTALGWDSPALTEREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"
#include <vector>

namespace siena
{

// should not be necessary
class ConfigurationTable;
class NetworkCache;

/**
 * This class defines the total geometrically weighted alter effect.
 */
class TotalGwdspAlterEffect : public NetworkDependentBehaviorEffect
{
public:
	TotalGwdspAlterEffect(const EffectInfo * pEffectInfo, bool forward);
	virtual void initialize(const Data * pData, State * pState,	int period,
			Cache * pCache);
	// do we need a deallocator?
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);


// protected:
// 	inline NetworkCache * pNetworkCache() const;


private:	
	// NetworkCache * lpNetworkCache;
	// ConfigurationTable * (NetworkCache::*lpTable)() const;
	double lparameter {};
	std::vector<double> lcumulativeWeight;
	double lforward {};
	double lweight {};
	double lexpmweight {};
	double lexpfactor {};
	ConfigurationTable *lpInitialisedTable;
};

// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

// should not be necessary


// NetworkCache * TotalGwdspAlterEffect::pNetworkCache() const
// {
// 	return this->lpNetworkCache;
// }

}

#endif /*TOTALGWDSPALTEREFFECT_H_*/
