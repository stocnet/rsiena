/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: PrimarySettingEffect.h
 *
 * Description: This file contains the definition of the
 * SettingSizeEffect class.
 *****************************************************************************/

#ifndef PRIMARYSETTINGEFFECT_H_
#define PRIMARYSETTINGEFFECT_H_

#include "NetworkWithPrimaryEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: PrimarySettingEffect class
// ----------------------------------------------------------------------------

/**
 * This class defines the setting size effects.
 */
 
class PrimarySettingEffect : public NetworkWithPrimaryEffect
{
public:
	PrimarySettingEffect(const EffectInfo * pEffectInfo,
		bool difference, bool logar, bool root, bool inv);
	virtual void initialize(const Data * pData, State * pState,	int period,
			Cache * pCache);

	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;
	virtual double egoStatistic(int ego, const Network * pNetwork);
	virtual bool egoEffect() const;	
protected:
	double transform(int value) const;
	virtual double tieStatistic(int alter);

private:
	double lparameter {};
	bool ldifference {};
	bool llogar {};
	bool lroot {};
	bool linv {};
	bool levalBoth {};
};

}

#endif /*PRIMARYSETTINGEFFECT_H_*/
