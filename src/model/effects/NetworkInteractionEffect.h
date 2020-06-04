/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkInteractionEffect.h
 *
 * Description: This file contains the definition of the
 * NetworkInteractionEffect class.
 *****************************************************************************/

#ifndef NETWORKINTERACTIONEFFECT_H_
#define NETWORKINTERACTIONEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

/**
 * A user-defined interaction effect between two or three other effects.
 */
class NetworkInteractionEffect: public NetworkEffect
{
public:
	NetworkInteractionEffect(const EffectInfo * pEffectInfo,
		NetworkEffect * pEffect1,
		NetworkEffect * pEffect2,
		NetworkEffect * pEffect3);
	virtual ~NetworkInteractionEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);

	virtual double calculateContribution(int alter) const;
	virtual bool egoEffect() const;

protected:
	virtual void initializeStatisticCalculation();
	virtual void onNextEgo(int ego);
	virtual double egoStatistic(int ego,
		const Network * pSummationTieNetwork);
	virtual double tieStatistic(int alter);
	virtual void cleanupStatisticCalculation();

private:
	// The interacting effects (lpEffect3 = 0 for two-way interactions)

	NetworkEffect * lpEffect1;
	NetworkEffect * lpEffect2;
	NetworkEffect * lpEffect3;
};

}

#endif /* NETWORKINTERACTIONEFFECT_H_ */
