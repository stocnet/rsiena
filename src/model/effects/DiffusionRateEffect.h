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
#include "NetworkEffect.h"
#include "BehaviorRateEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;
// class BehaviorVariable;
class ConstantCovariate;
class ChangingCovariate;
// class NetworkVariable;

// Remove: class DiffusionEffectValueTable;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * The base class for all diffusion rate effects depending on some network variable.
 */
class DiffusionRateEffect : public BehaviorRateEffect
{
public:
    DiffusionRateEffect(const EffectInfo* pEffectInfo);

    virtual ~DiffusionRateEffect();
    
    virtual void initialize(const Data* pData, 
        State* pState, int period, Cache* pCache);
    
    double calculateContribution(int i) const;

protected:
    inline const Network * pNetwork() const;
    
    // Helper method that calculates raw proximity/exposure statistic
    // Contains all effect-specific conditional logic
    // double proximityValue(const Network* pNetwork, int i, int period) const;
    // Only declare as pure virtual, do not implement here:
    virtual double proximityValue(const Network* pNetwork, int i) const;
    // Pure stateless utility used in most diffusionRate effects
    static double applyThreshold(double value, int numInfectedAlter, int absThreshold, bool capAtThreshold);

    // The covariates some effects depend on
    const ConstantCovariate * lpConstantCovariate;
    const ChangingCovariate * lpChangingCovariate;  

private:
     // The network variable this effect depends on
    const Network * lpNetwork;

};

// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the network this effect is interacting with.
 */
const Network * DiffusionRateEffect::pNetwork()
	const
{
	return this->lpNetwork;
}

}

#endif /* DIFFUSIONRATEEFFECT_H_ */