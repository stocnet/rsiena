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
    using DiffusionRateEffect::DiffusionRateEffect;
protected:
    double proximityValue(const Network* pNetwork, int i) const;
};

}

#endif /* EXPOSUREEFFECT_H_ */