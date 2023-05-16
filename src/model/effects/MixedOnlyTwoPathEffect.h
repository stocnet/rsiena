/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedOnlyTwoPathEffect.h
 *
 * Description: This file contains the definition of the
 * MixedOnlyTwoPathEffect class.
 *****************************************************************************/

#ifndef MIXEDONLYTWOPATHEFFECT_H_
#define MIXEDONLYTWOPATHEFFECT_H_

#include "MixedNetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Section: MixedOnlyTwoPathEffect class
// ----------------------------------------------------------------------------

/**
 * For toAny effect (see manual).
 */
class MixedOnlyTwoPathEffect : public MixedNetworkEffect
{
public:
	MixedOnlyTwoPathEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
			State * pState, int period, Cache * pCache);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double egoStatistic(int ego,
		const Network * pSummationTieNetwork);
};

}

#endif /* MIXEDONLYTWOPATHEFFECT_H_ */
