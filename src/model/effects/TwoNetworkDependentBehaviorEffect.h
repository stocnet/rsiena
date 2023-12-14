/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TwoNetworkDependentBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * TwoNetworkDependentBehaviorEffect class.
 *****************************************************************************/

#ifndef TWONETWORKDEPENDENTBEHAVIOREFFECT_H_
#define TWONETWORKDEPENDENTBEHAVIOREFFECT_H_

#include "BehaviorEffect.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;

// ----------------------------------------------------------------------------
// Section: TwoNetworkDependentBehaviorEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all behavior effects depending on two network variables.
 */
class TwoNetworkDependentBehaviorEffect : public BehaviorEffect
{
public:
	TwoNetworkDependentBehaviorEffect(const EffectInfo * pEffectInfo);
	virtual ~TwoNetworkDependentBehaviorEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

protected:
	inline const Network * pFirstNetwork() const;
	inline const Network * pSecondNetwork() const;
	double firstTotalAlterValue(int i) const;
	double firstTotalInAlterValue(int i) const;
	virtual void preprocessEgo(int ego);

private:
	// The network this effect is interacting with
	const Network * lpFirstNetwork;
	const Network * lpSecondNetwork;
	// total out- and in-alter values
	double * lfirstTotalAlterValues {};
	double * lfirstTotalInAlterValues {};
};

// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the networks this effect is interacting with.
 */
const Network * TwoNetworkDependentBehaviorEffect::pFirstNetwork()
	const
{
	return this->lpFirstNetwork;
}

const Network * TwoNetworkDependentBehaviorEffect::pSecondNetwork()
	const
{
	return this->lpSecondNetwork;
}

}

#endif /*NETWORKDEPENDENTBEHAVIOREFFECT_H_*/
