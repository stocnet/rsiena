/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: PrimaryCompressionEffect.h
 *
 * Description: This file contains the definition of the
 * SettingSizeEffect class.
 *****************************************************************************/

#ifndef PRIMARYCOMPRESSIONEFFECT_H_
#define PRIMARYCOMPRESSIONEFFECT_H_

#include "NetworkWithPrimaryEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: PrimaryCompressionEffect class
// ----------------------------------------------------------------------------

/**
 * This class defines the setting size effects.
 */
 
class PrimaryCompressionEffect : public NetworkWithPrimaryEffect
{
public:
	PrimaryCompressionEffect(const EffectInfo * pEffectInfo,
		bool inside, bool useSize);
	virtual void initialize(const Data * pData, State * pState,	int period,
			Cache * pCache);

	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;
	virtual double egoStatistic(int ego, const Network * pNetwork);
	virtual bool egoEffect() const;
	
protected:
	virtual double tieStatistic(int alter);

private:
	double lparameter {};
	bool linside {};
	bool luseSize {};
	
	// log((n-1-lprimDegree)/lparameter):
	double llogNonPrimary;
	// log((lprimDegree-outdegree(ego))/lparameter):
	double llogPrimary;
};

}

#endif /*PRIMARYCOMPRESSIONEFFECT_H_*/
