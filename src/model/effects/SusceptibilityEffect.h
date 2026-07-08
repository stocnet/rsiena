/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SusceptibilityEffect.h
 *
 * Description: This file contains the definition of the class
 * SusceptibilityEffect, which handles susceptibility diffusion rate effects.
 *****************************************************************************/

#ifndef SUSCEPTIBILITYEFFECT_H_
#define SUSCEPTIBILITYEFFECT_H_

#include "DiffusionRateEffect.h"

namespace siena {

class SusceptibilityEffect : public DiffusionRateEffect 
{
public:
    SusceptibilityEffect(const EffectInfo* pEffectInfo);

    void initialize(const Data* pData, State* pState, int period, Cache* pCache) override;

protected:
    double proximityValue(const Network* pNetwork, int i) const override;

private:
    int labsThreshold{0};
    bool lcapAtThreshold{false};
    bool lhasThreshold{false};

    bool lIsSusceptAvIn{false};
    bool lIsSusceptAvCovar{false};

    double applyThreshold(double value, int numInfectedAlter) const;
};

}

#endif /* SUSCEPTIBILITYEFFECT_H_ */