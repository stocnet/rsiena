/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SettingSizeEffect.h
 *
 * Description: This file contains the definition of the
 * SettingSizeEffect class.
 *****************************************************************************/

#ifndef SETTINGSIZEEFFECT_H_
#define SETTINGSIZEEFFECT_H_

#include "SettingsNetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: SettingSizeEffect class
// ----------------------------------------------------------------------------

/**
 * This class defines the setting size effects.
 */
 
class SettingSizeEffect : public SettingsNetworkEffect
{
public:
	SettingSizeEffect(const EffectInfo * pEffectInfo,
		bool difference,
		bool logar,
		bool root,
		bool inv,
		bool creation,
		bool evalDifference);

	virtual double calculateContribution(int alter) const;
	virtual double egoStatistic(int ego, const Network * pNetwork);

private:
	double lparameter {};
	bool ldifference {};
	bool llogar {};
	bool lroot {};
	bool linv {};
	bool lcreation {};
	bool levalDifference {};
	bool levalLog {};
	bool levalSqrt {};
};

}

#endif /*SETTINGSIZEEFFECT_H_*/
