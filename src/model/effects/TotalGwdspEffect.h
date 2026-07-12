/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TotalGwdspEffect.h
 *
 * Description: This file contains the declaration of the class TotalGwdspEffect.
 *
 *****************************************************************************/

#ifndef TOTALGWDSPEFFECT_H_
#define TOTALGWDSPEFFECT_H_

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
class TotalGwdspEffect : public NetworkDependentBehaviorEffect
{
public:
	TotalGwdspEffect(const EffectInfo * pEffectInfo, bool forward, bool nc);
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
	std::vector<int> lTwoPathCount;
	std::vector<int> lTouched;
	double lforward {};
	double lweight {};
	double lexpmweight {};
	double lexpfactor {};
	ConfigurationTable *lpInitialisedTable;
	bool lnc {};
};

// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

// should not be necessary


// NetworkCache * TotalGwdspffect::pNetworkCache() const
// {
// 	return this->lpNetworkCache;
// }

}

#endif /*TOTALGWDSPEFFECT_H_*/
