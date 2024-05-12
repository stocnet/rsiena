/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EffectInfo.h
 *
 * Description: This file contains the definition of the
 * EffectInfo class.
 *****************************************************************************/

#ifndef EFFECTINFO_H_
#define EFFECTINFO_H_

#include <string>

namespace siena
{

/**
 * Encapsulates the information relevant to an effect, except for the
 * basic rate effects.
 */
class EffectInfo
{
public:
	EffectInfo(std::string variableName,
		std::string effectName,
		std::string effectType,
		double parameter,
		double internalEffectParameter,
		std::string interactionName1,
		std::string interactionName2,
		std::string rateType);
	EffectInfo(std::string variableName,
		std::string effectName,
		std::string effectType,
		double parameter,
		const EffectInfo * pEffect1,
		const EffectInfo * pEffect2,
		const EffectInfo * pEffect3);

	void parameter(double value);

	std::string variableName() const;
	std::string effectName() const;
	std::string effectType() const;
	double parameter() const;
	double internalEffectParameter() const;
	std::string interactionName1() const;
	std::string interactionName2() const;
	std::string rateType() const;
	const EffectInfo * pEffectInfo1() const;
	const EffectInfo * pEffectInfo2() const;
	const EffectInfo * pEffectInfo3() const;

private:
	// The name of the variable this effect is associated with
	std::string lvariableName {};

	// A short name of the effect used to identify the semantics
	// of this effect
	std::string leffectName {};

	// The type of the effect ("rate", "eval", or "endow")
	std::string leffectType {};

	// The multiplicative weight in the respective function
	double lparameter {};

	// The internal parameter, if applicable
	double linternalEffectParameter {};

	// The name of a variable or covariate this effect is interacting with,
	// if applicable
	std::string linteractionName1 {};

	// The name of the other interaction variable or covariate, if the
	// effect has two such interactions
	std::string linteractionName2 {};

	// Distinguishes between structural rate effects and covariate rate effects
	std::string lrateType {};

	// The interacting effect descriptions.
	// Undefined for non-interaction effects.
	// The third effect is undefined for two-way interactions.
	const EffectInfo * lpEffectInfo1;
	const EffectInfo * lpEffectInfo2;
	const EffectInfo * lpEffectInfo3;
};

}

#endif /*EFFECTINFO_H_*/
