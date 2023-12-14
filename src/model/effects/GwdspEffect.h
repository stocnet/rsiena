/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwdspEffect.h
 *
 * Description: This file contains the declaration of the class GwdspEffect.
 *
 *****************************************************************************/

#ifndef GWDSPEFFECT_H_
#define GWDSPEFFECT_H_

#include "model/effects/NetworkEffect.h"
#include <vector>

namespace siena
{

class EgocentricConfigurationTable;
class NetworkCache;

/**
 * This class defines the geometrically weighted dyadwise shared partners effect.
 */
class GwdspEffect : public NetworkEffect
{
public:
	GwdspEffect(const EffectInfo * pEffectInfo,
			EgocentricConfigurationTable * (NetworkCache::*pTable)() const,
			bool forward);
	virtual void initialize(const Data * pData, State * pState,	int period,
			Cache * pCache);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double egoStatistic(int ego,
		const Network * pSummationTieNetwork);
	inline NetworkCache * pNetworkCache() const;


private:	
	NetworkCache * lpNetworkCache;
	EgocentricConfigurationTable * (NetworkCache::*lpTable)() const;
	double lparameter {};
	std::vector<double> lcumulativeWeight;
	double lforward {};
	double lweight {};
	double lexpmweight {};
	double lexpfactor {};
	EgocentricConfigurationTable *lpInitialisedTable;
};

// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

NetworkCache * GwdspEffect::pNetworkCache() const
{
	return this->lpNetworkCache;
}

}

#endif /*GWDSPEFFECT_H_*/
