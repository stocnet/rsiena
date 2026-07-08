/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Distance2ExposureEffect.h
 *
 * Description: This file contains the definition of the class
 * Distance2ExposureEffect, which handles distance-2 exposure rate effects.
 *****************************************************************************/

#ifndef DISTANCE2EXPOSUREEFFECT_H_
#define DISTANCE2EXPOSUREEFFECT_H_

#include "DiffusionRateEffect.h"

namespace siena
{

/**
 * Handles all distance-2 diffusion rate effects (e.g., anyInExposureDist2, 
 * totInExposureDist2, avTinExposureDist2, totAInExposureDist2).
 */
class Distance2ExposureEffect : public DiffusionRateEffect
{
public:
    Distance2ExposureEffect(const EffectInfo* pEffectInfo);
    void initialize(const Data* pData, State* pState, int period, Cache* pCache) override;
    
protected:
    double proximityValue(const Network* pNetwork, int i) const override;
private:
    double applyThreshold(double value, int numInfectedAlter) const;

    int labsThreshold{0};
    bool lcapAtThreshold{false};
    bool lhasThreshold{false};

    bool lIsAvTinExposureDist2{false};
    bool lIsTotAInExposureDist2{false};
    bool lIsAnyInExposureDist2{false};
};

}

#endif /* DISTANCE2EXPOSUREEFFECT_H_ */