/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InfectEffect.h
 *
 * Description: This file contains the definition of the class
 * InfectEffect, which handles infection-based diffusion rate effects.
 *****************************************************************************/

#ifndef INFECTEFFECT_H_
#define INFECTEFFECT_H_

#include "DiffusionRateEffect.h"

namespace siena {

class InfectEffect : public DiffusionRateEffect
{
public:
    InfectEffect(const EffectInfo* pEffectInfo);

    void initialize(const Data* pData, State* pState, int period, Cache* pCache) override;

protected:
    double proximityValue(const Network* pNetwork, int i) override;

private:
    // Threshold state variables
    int labsThreshold{0};
    bool lcapAtThreshold{false};
    bool lhasThreshold{false};

    // Optimization flags: pre-parsed routing switches to avoid string matching in the loop
    bool lIsInfectIn{false};
    bool lIsInfectDegOrOut{false};
    bool lIsInfectCovar{false};

    // Fast threshold application wrapper
    double applyThreshold(double value, int numInfectedAlter) const;
};

}

#endif /* INFECTEFFECT_H_ */