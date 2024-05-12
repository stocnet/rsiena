/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DiffusionRateEffect.h
 *
 * Description: This file contains the definition of the
 * DiffusionRateEffect class.
 *****************************************************************************/

#ifndef DIFFUSIONRATEEFFECT_H_
#define DIFFUSIONRATEEFFECT_H_

#include <string>

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class NetworkVariable;
class BehaviorVariable;
class DiffusionEffectValueTable;
class ConstantCovariate;
class ChangingCovariate;
class Network;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Encapsulates the information necessary for calculating the contributions
 * of a diffusion rate effect. This includes the effect type,
 * the network variable and the behavior variable
 * the effect depends on, and the statistical parameter of the effect.
 */
class DiffusionRateEffect
{
public:
	DiffusionRateEffect(const NetworkVariable * pVariable,
		const BehaviorVariable * pBehaviorVariable,
		std::string effectName,
		double parameter,
		double internalEffectParameter);
	DiffusionRateEffect(const NetworkVariable * pVariable,
		const BehaviorVariable * pBehaviorVariable,
		const ConstantCovariate * pCovariate,
		const ChangingCovariate * pChangingCovariate,
		std::string effectName,
		double parameter,
		double internalEffectParameter);

	virtual ~DiffusionRateEffect();
	double proximityValue(Network * pNetwork, int i, int egoNumer,
			int egoDenom) const;
	double value(int i, int period) const;
	void parameter(double parameterValue) const;
	double parameter() const;
	void internalEffectParameter(int parValue);
	int internalEffectParameter() const;

private:
	// The network variable this effect depends on
	const NetworkVariable * lpVariable;

	// The behavior variable this effect depends on
	const BehaviorVariable * lpBehaviorVariable;

	// The covariates some effects depend on
	const ConstantCovariate * lpConstantCovariate;
	const ChangingCovariate * lpChangingCovariate;

	// A table for efficient calculation of contributions. If two actors have
	// the same effect value, then the table ensures that we don't
	// calculate the same contribution twice.

	DiffusionEffectValueTable * lpTable;
	std::string leffectName {};
	int linternalEffectParameter {};
	int labsInternalEffectParameter {};
	bool linternalNonZero {};
};

}

#endif /* DIFFUSIONRATEEFFECT_H_ */
