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


namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;
class BehaviorVariable;
class ConstantCovariate;
class ChangingCovariate;
class NetworkVariable;

// Remove: class DiffusionEffectValueTable;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Encapsulates the information necessary for calculating the contributions
 * of a diffusion rate effect. This includes the effect type,
 * the network variable and the behavior variable
 * the effect depends on, and the statistical parameter of the effect.
 */
class DiffusionRateEffect : public Effect
{
public:
    DiffusionRateEffect(const EffectInfo* pEffectInfo,
                        const Network* pNetwork,
                        const BehaviorVariable* pBehaviorVariable,
                        std::string effectName,
                        double parameter,
                        double internalEffectParameter);
    DiffusionRateEffect(const EffectInfo* pEffectInfo,
        const Network * pNetwork,
        const BehaviorVariable * pBehaviorVariable,
        const ConstantCovariate * pCovariate,
        const ChangingCovariate * pChangingCovariate,
        std::string effectName,
        double parameter,
        double internalEffectParameter);

    virtual ~DiffusionRateEffect();
    
    // Returns the raw statistic (for scores and statistics)
    double value(int i, int period) const;
    
    // Returns the exponentiated rate contribution (for rate calculations)
    double logRate(int i, int period) const;
    
    void parameter(double parameterValue);
    double parameter() const;

    void setInternalEffectParameter(int parValue);
    int getInternalEffectParameter() const;
protected:
    // Helper method that calculates raw proximity/exposure statistic
    // Contains all effect-specific conditional logic
    // double proximityValue(const Network* pNetwork, int i, int period) const;
    // Only declare as pure virtual, do not implement here:
    virtual double proximityValue(const Network* pNetwork, int i, int period) const;

	double applyInternalEffectParameter(double value, int numInfectedAlter) const;

    // The network variable this effect depends on
    const Network * lpNetwork;

    // The behavior variable this effect depends on
    const BehaviorVariable * lpBehaviorVariable;

    // The covariates some effects depend on
    const ConstantCovariate * lpConstantCovariate;
    const ChangingCovariate * lpChangingCovariate;

    std::string leffectName {};
private:
    double lparameter {};
    int linternalEffectParameter {};
    int labsInternalEffectParameter {};
    bool linternalNonZero {};
};

}

#endif /* DIFFUSIONRATEEFFECT_H_ */