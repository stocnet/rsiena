/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ExposureEffect.h
 *
 * Description: This file contains the definition of the class
 * ExposureEffect, which handles direct exposure-based diffusion rate effects.
 *****************************************************************************/

#ifndef EXPOSUREEFFECT_H_
#define EXPOSUREEFFECT_H_

#include "DiffusionRateEffect.h"

namespace siena {

class ExposureEffect : public DiffusionRateEffect
{
public:
    ExposureEffect(const EffectInfo* pEffectInfo);

protected:
    double proximityValue(const Network* pNetwork, int i) const override;
private:
    int labsThreshold{0};
    bool lcapAtThreshold{false};
    bool lhasThreshold{false};
    bool lIsAvExposure{false};

    double applyThreshold(double value, int numInfectedAlter) const;
};

}

#endif /* EXPOSUREEFFECT_H_ */